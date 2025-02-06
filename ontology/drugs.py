import sys

sys.path.append(
    "/Users/pawan/Documents/Elucidata/Knowledge_Graphs/biomedical_ontology/biomedical_ontologies_kg/"
)
import pandas as pd
import tqdm
import requests
from dotenv import load_dotenv
import json
from scripts.create_db import (
    save_to_db,
    read_from_db,
    add_node_table_name,
    add_table_name,
)

load_dotenv()


def fetch_json(url):
    """Helper function to fetch JSON data from a URL"""
    response = requests.get(url)
    if response.status_code >= 202:
        return None
    return response.json()


def process_drug_gene_regulation(file_path, relation_value, table_meta):
    """Helper function to process drug-gene regulation data"""
    df = pd.read_csv(file_path)
    df["relation"] = relation_value

    genes = read_from_db("gene__nodes")
    df["gene_target"] = df["gene_target"].str.upper()
    df = pd.merge(
        df,
        genes[["gene_id", "gene_symbol"]],
        left_on="gene_target",
        right_on="gene_symbol",
    )[["gene_id", "drug_id"]]

    print(df.head())
    add_table_name(table_meta)
    save_to_db(df, table_meta[0])


def drug_node():
    with open("ontologies/pubchem_onto_upper.json") as f:
        data = json.load(f)

    dct = {"drug_id": [], "name": [], "synonyms": []}

    for drug in data:
        dct["drug_id"].append(drug["id"])
        dct["name"].append(drug["name"])
        syn = ";".join(drug["synonym"])
        dct["synonyms"].append(syn)

    drug__nodes = pd.DataFrame(dct)
    chebi_id = pd.read_csv(
        "ontologies/drug_chebi_ids.tsv", sep="\t", names=["drug_id", "chebi_id"]
    )
    chembl_id = pd.read_csv(
        "ontologies/drug_chembl_ids.tsv", sep="\t", names=["drug_id", "chembl_id"]
    )
    drug__nodes = pd.merge(drug__nodes, chebi_id, on="drug_id", how="left")
    drug__nodes = pd.merge(drug__nodes, chembl_id, on="drug_id", how="left")
    drug__nodes = drug__nodes.fillna("0")

    polly_drugs = pd.read_csv("ontologies/drugs_on_polly.csv")
    polly_drugs = list(polly_drugs["curated_drugs"])

    def filter_rows(row):
        if row["name"].lower() in polly_drugs:
            return 1
        else:
            syns = list(map(lambda x: x.lower(), row["synonyms"].split(";")))
            inter = set(polly_drugs).intersection(set(syns))
            if len(inter) == 0:
                return 0
            else:
                return 1

    drug__nodes["matched"] = drug__nodes.apply(filter_rows, axis=1)
    drug__nodes.to_csv("ontologies/drug__nodes.csv", index=False)
    df = pd.read_csv("ontologies/drug__nodes.csv")
    final_df = df[df["matched"] == 1]
    final_df.drop(columns=["matched"], inplace=True)

    save_to_db(final_df, "drug__nodes")
    add_node_table_name(
        ["drug__nodes", "drug", ";".join(final_df.columns), "not_mapped"]
    )


def up_drug_gene():
    process_drug_gene_regulation(
        file_path="ontologies/drug_upregulates_gene.csv",
        relation_value="drug_upregulates",
        table_meta=[
            "drug__upregulates__gene",
            "drug",
            "gene",
            "drug_upregulates",
            "not_mapped",
        ],
    )


def down_drug_gene():
    process_drug_gene_regulation(
        file_path="ontologies/drug_downregulates_gene.csv",
        relation_value="drug_downregulates",
        table_meta=[
            "drug__downregulates__gene",
            "drug",
            "gene",
            "drug_downregulates",
            "not_mapped",
        ],
    )


def get_drug_gene_interaction():
    """Retrieve drug-gene interaction data from PubChem and save to the database"""
    drugs_df = read_from_db("drug__nodes")
    compound_ids = drugs_df["drug_id"].tolist()
    c = ["1117", "24524", "10461508"]
    gene_ids_list = []

    for cid in tqdm.tqdm(compound_ids):
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/GeneID,PubMedID/JSON"
        response_json = fetch_json(url)
        if response_json is None:
            continue

        info = response_json.get("InformationList", {}).get("Information", [{}])[0]
        if "GeneID" in info:
            gids = info["GeneID"]
            gene_ids_list.append(gids if gids else ["None"])

    interactions_df = pd.DataFrame(
        list(zip(compound_ids, gene_ids_list)), columns=["drug_id", "gene_id"]
    )
    interactions_df = (
        interactions_df.explode("gene_id").dropna(subset=["gene_id"]).drop_duplicates()
    )

    # Map Entrez IDs to HGNC IDs
    hgnc_gene = read_from_db("gene__nodes")
    mapping_df = hgnc_gene[["gene_id", "entrez_id"]].copy()
    mapping_df.columns = ["hgnc_id", "entrez_id"]
    interactions_df = pd.merge(
        interactions_df,
        mapping_df,
        left_on="gene_id",
        right_on="entrez_id",
        how="inner",
    )
    interactions_df = interactions_df[["drug_id", "hgnc_id"]]
    interactions_df.columns = ["drug_source", "gene_target"]
    interactions_df["relation"] = "interacts_with"

    add_table_name(
        [
            "drug__interacts_with__relation",
            "drug",
            "gene",
            "interacts_with",
            "not_mapped",
        ]
    )
    save_to_db(interactions_df, "drug__interacts_with__relation")


def get_3d_similar_drugs():
    """Retrieve 3D structurally similar drugs from PubChem and save to the database."""
    drugs_df = read_from_db("drug__nodes")
    compound_ids = drugs_df["drug_id"].tolist()
    c = ["2244", "5282452", "10461508"]
    analog_ids_list = []

    for cid in tqdm.tqdm(compound_ids):
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_3d/cid/{cid}/cids/JSON"
        response_json = fetch_json(url)
        if response_json is None:
            continue

        identifier_list = response_json.get("IdentifierList", {})
        if "CID" in identifier_list:
            cids = identifier_list["CID"]
            cid_int = int(cid)
            if cid_int in cids:
                cids.remove(cid_int)
            analog_ids_list.append(cids)
        else:
            analog_ids_list.append(None)

    interactions_df = pd.DataFrame(
        list(zip(compound_ids, analog_ids_list)), columns=["drug_id", "similar_drug_id"]
    )
    interactions_df = interactions_df.dropna(subset=["similar_drug_id"]).explode(
        "similar_drug_id"
    )
    interactions_df.columns = ["drug_source", "drug_target"]
    interactions_df = pd.merge(
        interactions_df,
        drugs_df[["drug_id"]],
        left_on="drug_target",
        right_on="drug_id",
        how="inner",
    ).drop(columns=["drug_id"])
    interactions_df = interactions_df.drop_duplicates()
    interactions_df["relation"] = "has_similar_structure"

    add_table_name(
        [
            "drug__has_similar_structure__relation",
            "drug",
            "drug",
            "has_similar_structure",
            "not_mapped",
        ]
    )
    save_to_db(interactions_df, "drug__has_similar_structure__relation")


def main_drug():
    # Uncomment the functions you wish to run:
    # drug_node()
    # down_drug_gene()
    # up_drug_gene()
    # get_drug_gene_interaction()
    get_3d_similar_drugs()


if __name__ == "__main__":
    main_drug()
