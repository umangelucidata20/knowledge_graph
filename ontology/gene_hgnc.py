"""
This module collects gene information using the HGNC API.
"""

import math
import json
from multiprocessing import Process, Manager

import pandas as pd
import httplib2 as http
from dotenv import load_dotenv
from scripts.create_db import save_to_db, read_from_db, add_node_table_name

# sys.path.append("/Users/pawan/Documents/Elucidata/Knowledge_Graphs/biomedical_ontology/biomedical_ontologies_kg/")

try:
    from urlparse import urlparse
except ImportError:
    from urllib.parse import urlparse


def format_dict(dct):
    """
    Converts any list values in the dictionary to a semicolon-separated string.
    """
    return {
        key: (
            ";".join(str(ele) for ele in dct[key])
            if isinstance(dct[key], list)
            else dct[key]
        )
        for key in dct
    }


def retrive_gene_data(sym):
    headers = {
        "Accept": "application/json",
    }
    uri = "http://rest.genenames.org"
    path = "/fetch/symbol/" + sym

    target = urlparse(uri + path)
    method = "GET"
    body = ""
    h = http.Http()
    response, content = h.request(target.geturl(), method, body, headers)

    if response["status"] == "200":
        data = json.loads(content)
        if len(data["response"]["docs"]) > 0:
            return format_dict(data["response"]["docs"][0])
        else:
            return " "

    else:
        print("Error detected: " + response["status"])


def get_gene_nodes(gene_lst, num, return_dict):
    lst = []
    for gene in gene_lst:
        node_data = retrive_gene_data(gene)
        if node_data != " ":
            lst.append(node_data)
    return_dict[num] = lst


def main_genes_hgnc():
    """
    Main function to retrieve gene nodes concurrently, update the gene node table,
    and update the disease-to-gene relation table in the database.
    """
    dis_gene = read_from_db("disease__associated_with__relation")
    gene_lst = list(set(dis_gene["gene_target"]))

    num_processes = 40
    return_dict = Manager().dict()
    jobs = []

    chunk_size = math.ceil(len(gene_lst) / num_processes)
    
    for i in range(num_processes):
        lb = i * chunk_size
        up = min((i + 1) * chunk_size, len(gene_lst))
        p = Process(target=get_gene_nodes, args=(gene_lst[lb:up], i, return_dict))
        jobs.append(p)
        p.start()

    for process in jobs:
        process.join()

    dfs = [pd.DataFrame(lst) for lst in return_dict.values()]
    final_df = pd.concat(dfs, ignore_index=True)

    # final_df = pd.read_csv("gene__nodes.csv")

    required_cols = [
        "hgnc_id",
        "symbol",
        "name",
        "prev_name",
        "alias_symbol",
        "location",
        "date_modified",
        "ena",
        "entrez_id",
        "mgd_id",
        "pubmed_id",
        "refseq_accession",
        "vega_id",
        "ensembl_gene_id",
        "ccds_id",
        "locus_group",
        "omim_id",
        "uniprot_ids",
        "ucsc_id",
        "rgd_id",
        "location_sortable",
        "agr",
        "mane_select",
        "enzyme_id",
        "gene_group",
        "prev_symbol",
        "alias_name",
        "pseudogene.org",
        "imgt",
    ]

    final_df = final_df[required_cols]

    new_names = [
        "gene_id",
        "gene_symbol",
        "name",
        "prev_name",
        "alias_symbol",
        "location",
        "date_modified",
        "ena",
        "entrez_id",
        "mgd_id",
        "pubmed_id",
        "refseq_accession",
        "vega_id",
        "ensembl_gene_id",
        "ccds_id",
        "locus_group",
        "omim_id",
        "uniprot_ids",
        "ucsc_id",
        "rgd_id",
        "location_sortable",
        "agr",
        "mane_select",
        "enzyme_id",
        "gene_group",
        "prev_symbol",
        "alias_name",
        "pseudogene_org",
        "imgt",
    ]
    final_df.columns = new_names
    add_node_table_name(
        ["gene__nodes", "gene", ";".join(final_df.columns), "not_mapped"]
    )
    save_to_db(final_df, "gene__nodes")

    genes_present = final_df[["gene_symbol"]]
    genes_present.columns = ["gene_target"]

    dis_gene = pd.merge(dis_gene, genes_present, on="gene_target", how="inner")
    save_to_db(dis_gene, "disease__associated_with__relation")


if __name__ == "__main__":
    main_genes_hgnc()
