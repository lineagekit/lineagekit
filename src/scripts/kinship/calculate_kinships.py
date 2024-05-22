import sys

from src.basic.Pedigree import Pedigree
from src.utility.utility import *
import time

sys.path.append(os.path.join(os.path.dirname(__file__), "../../.."))

filepath = get_file_path("Specify the path to the pedigree:")
pedigree = Pedigree.get_pedigree_graph_from_file(filepath=filepath, separation_symbol=" ",
                                                 missing_parent_notation=["-1"], skip_first_line=False)

# pedigree = Pedigree.get_pedigree_graph_from_file(filepath=filepath, separation_symbol=";",
#                                                  missing_parent_notation=[""], skip_first_line=True)

print("Calculating kinships")
start_time = time.time()
# kinship_matrix = pedigree.calculate_kinship_for_probands_sparse()
# (vertex_to_index, direct_kinship_matrix) = pedigree.calculate_kinship()
# half_kinship_matrix = pedigree.calculate_kinship_for_probands_half()
# kinship_matrix = pedigree.calculate_kinship_by_levels()
kinship_matrix, vertex_to_representative = pedigree.calculate_sparse_kinship_queue()
# iterative_kinships = pedigree.calculate_sparse_kinship_queue()
end_time = time.time()
print(f"Time taken after compilation: {end_time - start_time}")
# print(kinship_matrix)


# for proband in pedigree.get_sink_vertices():
#     for other_proband in pedigree.get_sink_vertices():
#         simple_algorithm_kinship = direct_kinship_matrix[vertex_to_index[proband], vertex_to_index[other_proband]]
#         # old_kinship = kinship_matrix[proband][other_proband]
#         # half_algorithm_kinship = iterative_kinships[proband][other_proband]
#         equivalence_kinship = kinship_matrix[vertex_to_representative[proband]][vertex_to_representative[other_proband]]
#         assert simple_algorithm_kinship == equivalence_kinship

print("Kinships have been calculated")


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
