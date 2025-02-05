import os
import sys
import pronto
import pandas as pd
from dotenv import load_dotenv

from scripts.create_db import read_from_db, save_to_db, add_table_name, add_node_table_name
from get_all_terms import *

load_dotenv()

def get_polly_cl(indexes):
    """
    Retrieve cell line terms from Polly based on the repository indexes.
    """
    term_repos = {}  # term: {repos}
    field = 'cell_line'
    polly_env = 'prod'
    
    for i, oa in enumerate(list(indexes.values())):
        if i % 4 == 0:
            allz = get_all_cell_lines_polly(field, oa, polly_env)
        
        for w in {d.strip() for d in allz}:
            term_repos[w] = term_repos.get(w, set()) | {oa.replace("_files", "")}
    
    term_repos.pop("", None)
    return term_repos.keys()

def cell_line_node(clo_v3):
    """
    Create nodes table for cell lines from the cellosaurus ontology.
    """
    gender_terms = [
        "Male",
        "Female",
        "Mixed_sex",
        "Sex_ambiguous",
        "Sex_unspecified",
    ] 
    
    df = pd.DataFrame(
        [
            (
                term.id, 
                term.name, 
                list(term.subsets.intersection(gender_terms)),
                list(term.subsets.difference(gender_terms)), 
                term.synonyms, 
                'cell-line'
            )
            for term in clo_v3.terms()
        ],
        columns=['cell_line_id', 'name', 'gender', 'category', 'synonyms', 'type']
    )
    
    df['synonyms'] = df["synonyms"].apply(lambda x: ';'.join([str(i).split("'")[1] for i in x]))
    df['gender'] = df['gender'].apply(lambda x: x[0] if len(x) > 0 else " ")
    df['category'] = df["category"].apply(lambda x: x[0] if len(x) > 0 else " ")
    
    df = df[df['name'].isin(polly_cell_lines)]
    
    add_node_table_name(['cell_line__nodes', 'cell_line', ';'.join(df.columns), 'not_mapped'])
    save_to_db(df, "cell_line__nodes")


def cell_line_rel(clo_v3):
    """
    Create relationship files for cell lines from the cellosaurus ontology.
    """
    rels_clo = [r.id for r in clo_v3.relationships()]
    clo_id, clo_rel, clo_tar = [], [], []
    
    for term in clo_v3.terms():
        for rel in rels_clo:
            if rel in str(list(term.relationships)):
                ts = term.relationships[clo_v3.get_relationship(rel)]
                clo_id.append(str(term.id))
                clo_rel.append(rel)
                clo_tar.append(str(list(ts)[0].id))
    
    rel_df = pd.DataFrame(list(zip(clo_id, clo_rel, clo_tar)), 
                          columns=['cell_line_source', 'relation', 'cell_line_target'])
    rel_df = rel_df.explode('cell_line_target').dropna(subset=['cell_line_target'])
    
    node_df = read_from_db("cell_line__nodes")[["cell_line_id"]]
    node_df.columns = ['cell_line_target']
    rel_df = pd.merge(node_df, rel_df, on='cell_line_target', how='inner')
    
    for rel in rels_clo:
        df = rel_df[rel_df['relation'] == rel]
        add_table_name([f"cell_line__{rel}__relation", "cell_line", "cell_line", rel, 'not_mapped'])
        save_to_db(df, f"cell_line__{rel}__relation")
        

def disease_cell_line(clo_2):
    """
    Create a table mapping cell lines to diseases based on the ontology.
    """
    clo_all = []
    for term in clo_2.terms():
        for xref in term.xrefs:
            if xref.id.startswith("MESH"):
                clo_all.append((term.id, xref.id, term.name))
    
    dis = pd.DataFrame(clo_all, columns=['cell_line_source', 'disease_target', "cell_line_name"])
    dis = dis[dis['cell_line_name'].isin(polly_cell_lines)]
    dis['relation'] = 'obtained_from_sample_with_disease'
    
    add_table_name([
        'cell_line__obtained_from_sample_with_disease__relation',
        'cell_line',
        'disease',
        'obtained_from_sample_with_disease',
        'not_mapped'
    ])
    save_to_db(dis, "cell_line__obtained_from_sample_with_disease__relation")


def main_cellosaurus(clo_v3):
    """
    Process cell line nodes and relationships from cellosaurus ontology.
    """
    ontology, new = clo_v3
    if new:
        cell_line_node(ontology)
        cell_line_rel(ontology)


def main_cell_line(clo_2):
    """
    Process disease and cell line relationships.
    """
    ontology, new = clo_2
    if new:
        disease_cell_line(ontology)


if __name__ == "__main__":

    refresh_token = os.environ.get('PROD_REFRESH_TOKEN')
    library_client = OmixAtlas(refresh_token, polly_env='prod')
    
    repo_id_index = get_all_oas(library_client)
    polly_cell_lines = get_polly_cl(repo_id_index)
    print(polly_cell_lines)
    
    main_cell_line((pronto.Ontology("ontologies/cell_line.obo"), True))
    main_cellosaurus((pronto.Ontology("ontologies/cell_line_cellosaurus.obo"), True))