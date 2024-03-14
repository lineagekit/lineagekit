import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), "../../.."))
from src.basic.GenGraph import *


def get_file_path(input_request: str):
    while True:
        file_path = input(input_request)
        if not os.path.exists(file_path):
            print("The specified file does not exist, try again")
        elif not os.path.isfile(file_path):
            print("The specified path is not a file, try again")
        else:
            return file_path


filepath = get_file_path("Specify the path to the pedigree:")
pedigree = GenGraph.get_haploid_graph_from_file(filepath=filepath, separation_symbol=";",
                                                missing_parent_notation=[""], skip_first_line=True)
# Get the ascending genealogy grouped by levels
ascending_genealogy_by_levels = pedigree.get_ascending_genealogy_from_vertices_by_levels(vertices=[2, 3])
# Alternatively, just get all the vertices in the ascending genealogy
ascending_genealogy = pedigree.get_ascending_graph_from_vertices(vertices=[2, 3])
# Perturb the edges
pedigree.remove_edge(parent=6, child=2)
pedigree.add_edge(parent=8, child=2)
# Recalculate the ascending genealogy
modified_ascending_genealogy = pedigree.get_ascending_genealogy_from_vertices_by_levels(vertices=[2, 3])
is_founder = pedigree.is_founder(10)
verify_pedigree_correctness = pedigree.verify_max_parents_number(2)
every_vertex_has_only_one_parent = pedigree.verify_max_parents_number(1)