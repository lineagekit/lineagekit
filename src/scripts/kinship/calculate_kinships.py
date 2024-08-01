import itertools

from basic.Pedigree import Pedigree
from utility.utility import *
import time

filepath = get_file_path("Specify the path to the pedigree:")

pedigree = Pedigree.get_pedigree_graph_from_file(filepath=filepath, separation_symbol=" ",
                                                 missing_parent_notation=["-1"], skip_first_line=True)

print("Calculating kinships")
start_time = time.time()
print("All pair-wise:")
(vertex_to_index, direct_kinship_matrix) = pedigree.calculate_kinship()
end_time = time.time()
print(f"Time taken: {end_time - start_time}")
print("Probands:")
kinship_matrix_c = pedigree.calculate_probands_kinship()
end_time = time.time()
print(f"Time taken: {end_time - start_time}")

for proband in pedigree.get_sink_vertices():
    for other_proband in pedigree.get_sink_vertices():
        simple_algorithm_kinship = direct_kinship_matrix[vertex_to_index[proband], vertex_to_index[other_proband]]
        kinship_c_implementation = kinship_matrix_c.get_kinship(proband, other_proband)

for vertex_1, vertex_2 in itertools.combinations(pedigree, 2):
    kinship = direct_kinship_matrix[vertex_to_index[vertex_1], vertex_to_index[vertex_2]]

print("Kinships have been calculated")
