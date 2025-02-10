# This code fetches properties of genes from Harmonizome api
import math
import json
import requests
import pandas as pd
from multiprocessing import Process, Manager
from scripts.create_db import save_to_db, read_from_db, add_node_table_name

DEFAULT_GENE_RESPONSE = {
    "symbol": "",
    "synonyms": [],
    "name": "",
    "description": "",
    "ncbiEntrezGeneId": -1,
    "ncbiEntrezGeneUrl": "",
    "proteins": [],
    "hgncRootFamilies": [],
}


def return_json(url):
    """
    Fetches JSON data from the provided URL.
    If the request fails or returns a non-200 status code, returns a default gene response.
    """
    try:
        response = requests.get(url)
        if response.status_code == 200:
            return json.loads(response.text)
    except Exception:
        print("Timed out")
        print(url)
    return DEFAULT_GENE_RESPONSE


def join_items(items, sep, key=None):
    """
    Joins a list of items into a string separated by 'sep'.
    If 'key' is provided, it extracts that key's value from each dictionary in the list.
    """
    if not items:
        return ""
    if key:
        return sep.join(item.get(key, "") for item in items)
    return sep.join(items)


def gene_properties(gene_lst, num, return_dict):
    base_url = "https://maayanlab.cloud/Harmonizome/api/1.0/gene/"

    gene_prop_dict = {
        "gene_id": [],
        "synonyms": [],
        "name": [],
        "definition": [],
        "NcbiEntrezGeneId": [],
        "NcbiEntrezGeneUrl": [],
        "proteins": [],
        "HgncRootFamilies": [],
    }

    for gene in gene_lst:
        url = base_url + gene
        gene_props = return_json(url)
        gene_prop_dict["gene_id"].append(gene)
        gene_prop_dict["synonyms"].append(
            join_items(gene_props.get("synonyms", []), ";")
        )
        gene_prop_dict["name"].append(gene_props.get("name", ""))
        gene_prop_dict["definition"].append(gene_props.get("description", ""))
        gene_prop_dict["NcbiEntrezGeneId"].append(
            gene_props.get("ncbiEntrezGeneId", -1)
        )
        gene_prop_dict["NcbiEntrezGeneUrl"].append(
            gene_props.get("ncbiEntrezGeneUrl", "")
        )
        gene_prop_dict["proteins"].append(
            join_items(gene_props.get("proteins", []), ",", key="symbol")
        )
        gene_prop_dict["HgncRootFamilies"].append(
            join_items(gene_props.get("hgncRootFamilies", []), ",", key="name")
        )

    return_dict[num] = gene_prop_dict


def retrive_gene_properties(gene_lst):
    num_processes = 20

    return_dict = Manager().dict()
    jobs = []
    chunk_size = math.ceil(len(gene_lst) / num_processes)

    for i in range(num_processes):
        lb = i * chunk_size
        up = min((i + 1) * chunk_size, len(gene_lst))
        if lb >= len(gene_lst):
            break

        p = Process(target=gene_properties, args=(gene_lst[lb:up], i, return_dict))
        jobs.append(p)
        p.start()

    for process in jobs:
        process.join()

    dfs = [pd.DataFrame(return_dict[i]) for i in sorted(return_dict.keys())]
    df = pd.concat(dfs, ignore_index=True)
    return df


def main_genes():
    # present in disease gene associations
    df = read_from_db("disease__associated_with__relation")
    gene_lst = list(set(df["gene_target"]))

    gene_nodes = retrive_gene_properties(gene_lst)  # remove the list slice
    gene_nodes["type"] = "gene"
    add_node_table_name(
        ["gene__nodes", "gene", ";".join(gene_nodes.columns), "not_mapped"]
    )
    save_to_db(gene_nodes, "gene__nodes")


# main
if __name__ == "__main__":
    main_genes()
