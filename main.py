from dotenv import load_dotenv
import argparse

from ontology import *
from ontology.drugs import main_drug
from ontology.gene_hgnc import main_genes_hgnc
from ontology.pathways import main_pathway

from scripts.load import load_file, update_mapping
from scripts.ontop_cli import materialize
from scripts.neo_export import ontology_export
from scripts.update_owl import update_owl


class OntologyGroup:
    def load_ontologies(self):
        self.cell_line_onto = load_file("cell_line")
        self.cellosaurus_onto = load_file("cell_line_cellosaurus")
        self.mesh_onto = load_file("disease")
        self.tissue_onto = load_file("tissue")
        
        # self.cell_line_onto = (pronto.Ontology('ontologies/cell_line.obo'),True)
        # self.cellosaurus_onto = (pronto.Ontology('ontologies/cell_line_cellosaurus.obo'),True)
        # self.mesh_onto = (pronto.Ontology('ontologies/disease.obo'),True)
        # self.tissue_onto = (pronto.Ontology('ontologies/tissue.obo'),True)

    def start(self):
        load_dotenv()

        main_cell_line(self.cell_line_onto)
        main_cellosaurus(self.cellosaurus_onto)
        main_mesh_disease(self.mesh_onto)
        main_disease_gene()
        # main_genes()
        main_genes_hgnc()
        main_load_tissue_cell_type(self.tissue_onto)
        main_tissue_cell_line_rel(self.tissue_onto)
        main_drug()
        main_pathway()

def process_ontologies():
    """
    Note: uncomment this function call in main() if ontology processing is required.
    """
    ontology_group = OntologyGroup()
    ontology_group.load_ontologies()
    ontology_group.start()


def update_semantics(mapping_file: str, owl_file: str):
    update_mapping(mapping_file)
    update_owl(owl_file)


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output_file")
    return parser.parse_args()

def main():
    args = parse_arguments()
    args = vars(args)

    # Uncomment the following line if ontology loading and processing is required.
    # process_ontologies()

    mapping_file_name = "semantics/BioMedOnto.obda"
    owl_file_name = "semantics/BioMedOnto.owl"

    update_semantics(mapping_file_name, owl_file_name)

    # Materialize the ontop using provided files and arguments
    materialize(
        mapping=mapping_file_name,
        ontology=owl_file_name,
        output=args["output_file"],
        properties="semantics/BioMedOnto.properties"
    )

    # Uncomment if export is needed:
    # ontology_export(args.output_file)


if __name__ == '__main__':
    main()