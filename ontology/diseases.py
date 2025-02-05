# This code reads the disease ontology file and saves it to db
import pronto
import pandas as pd
from scripts.load import load_file
from scripts.create_db import (
    add_table_name,
    read_from_db,
    save_to_db,
    add_node_table_name,
)


def get_synonyms(synonyms):
    """
    Converts a collection of synonyms into a semicolon-separated string.
    """
    syn = list(synonyms)
    temp = [str(str(s).split("'")[1]) for s in syn]
    return ";".join(temp)


def get_subclass(term):
    """
    Retrieves direct subclass IDs (distance=1) for the given term.
    """
    mesh = term.subclasses(with_self=False, distance=1)
    return [str(m.id) for m in mesh]


def disease_nodes(mesh):
    """
    Reads the disease ontology terms, extracts relevant node information,
    and saves the nodes to the database.
    """
    data = []
    for term in mesh.terms():
        term_id = str(term.id)
        term_name = str(term.name)
        term_synonyms = get_synonyms(term.synonyms)
        term_type = "disease"
        data.append((term_id, term_name, term_synonyms, term_type))

    df = pd.DataFrame(data, columns=["disease_id", "name", "synonyms", "type"])

    table_columns = ";".join(df.columns)
    add_node_table_name(["disease__nodes", "disease", table_columns, "not_mapped"])
    save_to_db(df, "disease__nodes")


def disease_subclass(mesh):
    """
    Reading the disease ontology terms, extracts subclass relationships,
    and saves the relations to the database.
    """
    data = []
    for term in mesh.terms():
        target_id = term.id
        subclass_list = get_subclass(term)
        data.append((target_id, subclass_list))

    sub_df = pd.DataFrame(data, columns=["disease_target", "disease_source"])
    sub_df = sub_df.explode("disease_source")
    sub_df = sub_df.dropna(subset=["disease_source"])
    sub_df["relation"] = "is_a"

    node_df = read_from_db("disease__nodes")[["disease_id"]]
    node_df = node_df.rename(columns={"disease_id": "disease_source"})

    merged_df = pd.merge(node_df, sub_df, on="disease_source", how="inner")

    add_table_name(
        ["disease__is_a__relation", "disease", "disease", "is_a", "not_mapped"]
    )
    save_to_db(merged_df, "disease__is_a__relation")


def main_mesh_disease(mesh):
    """
    Processes the disease [new]ontology.
    """
    ontology, new = mesh
    if new:
        disease_nodes(ontology)
        disease_subclass(ontology)


if __name__ == "__main__":
    ontology = pronto.Ontology("ontologies/disease.obo")
    main_mesh_disease((ontology, True))
