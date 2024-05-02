import sys

from src.basic.Pedigree import Pedigree
from src.utility.utility import *

sys.path.append(os.path.join(os.path.dirname(__file__), "../../.."))

filepath = get_file_path("Specify the path to the pedigree:")
pedigree = Pedigree.get_pedigree_graph_from_file(filepath=filepath, separation_symbol=" ",
                                                 missing_parent_notation=["-1"], skip_first_line=False)

#(vertex_to_index, kinship_matrix) = pedigree.calculate_kinship()
#print(vertex_to_index)
#print(kinship_matrix)

sparse_kinship_matrix = pedigree.calculate_kinship_for_probands_sparse()
print(sparse_kinship_matrix)

# for proband in pedigree.get_sink_vertices():
#     for other_proband in pedigree.get_sink_vertices():
#         first_proband_index = vertex_to_index[proband]
#         second_proband_index = vertex_to_index[other_proband]
#         simple_algorithm_kinship = kinship_matrix[first_proband_index, second_proband_index]
#         improved_algorithm_kinship = sparse_kinship_matrix[proband][other_proband]
#         assert simple_algorithm_kinship == improved_algorithm_kinship


def get_vertex_id():
    while True:
        try:
            vertex_id = int(input("Enter the vertex id: "))
            if vertex_id not in pedigree:
                print("Invalid vertex id")
            return vertex_id
        except Exception:
            pass


# while True:
#     first_vertex = get_vertex_id()
#     second_vertex = get_vertex_id()
#     first_vertex_index = vertex_to_index[first_vertex]
#     second_vertex_index = vertex_to_index[second_vertex]
#     print(kinship_matrix[first_vertex_index, second_vertex_index])
