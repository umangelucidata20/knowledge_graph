"""
This code gets disease-gene associations using the Harmonizome API.
"""

import sys
import math
import json
import requests
import pandas as pd
from bs4 import BeautifulSoup
from multiprocessing import Process, Manager
from dotenv import load_dotenv

sys.path.append(
    "/Users/pawan/Documents/Elucidata/Knowledge_Graphs/biomedical_ontology/biomedical_ontologies_kg/"
)

from scripts.harmonizomeapi import Harmonizome, Entity
from scripts.create_db import save_to_db, add_table_name, read_from_db

load_dotenv()


def return_json(addr):
    """Get JSON data from a URL."""
    base_url = "https://maayanlab.cloud/Harmonizome"
    url = base_url + addr
    try:
        response = requests.get(url)
    except Exception as e:
        print("Error fetching:", addr)
        return " "

    if response.status_code == 200:
        try:
            return json.loads(response.text)
        except Exception:
            return " "
    else:
        print("Error fetching:", addr)
        return " "


def convert_name_to_url(name):
    """
    Convert a disease name into two possible URLs for fetching gene-disease associations.
    """
    url_variant1 = name.replace(",", "%").replace(" ", "+")
    url_variant2 = name.replace(",", "%2C").replace(" ", "+")
    url1 = f"https://maayanlab.cloud/Harmonizome/gene_set/{url_variant1}/CTD+Gene-Disease+Associations"
    url2 = f"https://maayanlab.cloud/Harmonizome/gene_set/{url_variant2}/CTD+Gene-Disease+Associations"
    return url1, url2


def fetch_mesh_id_from_url(url):
    """Fetch and return the MeSH ID from the provided URL, if available."""
    try:
        response = requests.get(url)
        if response.status_code == 200:
            soup = BeautifulSoup(response.text, "html.parser")
            for link in soup.find_all("a"):
                link_url = link.get("href")
                if link_url and "ctdbase.org" in link_url:
                    return link_url.split("=")[-1]
    except Exception:
        pass
    return " "


def find_mesh_id(url_tuple):
    """
    Try both URL variants to extract the MeSH ID.
    Returns the first valid MeSH ID found or a blank string if none is found.
    """
    for url in url_tuple:
        mesh_id = fetch_mesh_id_from_url(url)
        if mesh_id != " ":
            return mesh_id
    return " "


def disease_gene_associations(genesets, num, return_dict):
    """
    Process a subset of gene sets, retrieving disease MeSH IDs and associated gene symbols.
    The results are stored in a shared dictionary.
    """
    associations = {"disease_id": [], "Associated_Gene_Symbols": []}

    for geneset in genesets:
        disease_name = geneset["name"].split("/")[0]

        url_tuple = convert_name_to_url(disease_name)
        mesh_id = find_mesh_id(url_tuple)

        gene_info = return_json(geneset["href"])
        if mesh_id != " " and gene_info != " ":
            associations["disease_id"].append(mesh_id)
            genes = [
                gene["gene"]["symbol"] for gene in gene_info.get("associations", [])
            ]
            associations["Associated_Gene_Symbols"].append(genes)

    return_dict[num] = associations


def retrive_dis_gene_associations(dataset):
    """
    Retrieve disease-gene associations using multiple processes and return the consolidated DataFrame.
    """
    num_processes = 20
    manager = Manager()
    return_dict = manager.dict()
    jobs = []
    dataset_length = len(dataset)
    chunk_size = math.ceil(dataset_length / num_processes)

    for i in range(num_processes):
        lb = i * chunk_size
        up = min((i + 1) * chunk_size, dataset_length)
        p = Process(
            target=disease_gene_associations, args=(dataset[lb:up], i, return_dict)
        )
        jobs.append(p)
        p.start()

    for process in jobs:
        process.join()

    # Combine all DataFrames from the individual processes
    df_list = [pd.DataFrame(return_dict[i]) for i in range(len(return_dict))]
    return pd.concat(df_list, ignore_index=True)


def main_disease_gene():
    """
    Main function to retrieve, filter, and store disease-gene associations.
    """
    # Get the list of datasets from the Harmonizome API and select the CTD dataset (index 24)
    dataset_lst = Harmonizome.get(Entity.DATASET)
    dataset = return_json(dataset_lst["entities"][24]["href"])["geneSets"]

    dis_gene = retrive_dis_gene_associations(dataset)

    # Filter: Only include diseases present in the disease ontology database
    mesh_dis = read_from_db("disease__nodes")
    valid_diseases = mesh_dis[["disease_id"]]
    dis_gene = pd.merge(dis_gene, valid_diseases, on="disease_id", how="inner")

    # Filter: Only include genes present in the gene ontology database with HGNC symbols
    hgnc_gene = read_from_db("gene__nodes")
    gene_df = hgnc_gene[["gene_id", "gene_symbol"]]

    # Explode the gene symbol list and drop any empty entries
    dis_gene = dis_gene.explode("Associated_Gene_Symbols")
    dis_gene = dis_gene.dropna(subset=["Associated_Gene_Symbols"])

    # Merge to keep only valid gene associations
    final_df = pd.merge(
        dis_gene,
        gene_df,
        left_on="Associated_Gene_Symbols",
        right_on="gene_symbol",
        how="inner",
    )
    final_df = final_df[["disease_id", "gene_id"]]
    final_df.columns = ["disease_source", "gene_target"]
    final_df["relation"] = "associated_with"

    # Add metadata and save the final DataFrame to the database
    add_table_name(
        [
            "disease__associated_with__relation",
            "disease",
            "gene",
            "associated_with",
            "not_mapped",
        ]
    )
    save_to_db(final_df, "disease__associated_with__relation")


# main
if __name__ == "__main__":
    main_disease_gene()
