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

# Use both matrices to get the kinship between every pair of probands
for proband in pedigree.get_sink_vertices():
    for other_proband in pedigree.get_sink_vertices():
        simple_algorithm_kinship = direct_kinship_matrix[vertex_to_index[proband], vertex_to_index[other_proband]]
        proband_kinship = kinship_matrix_c.get_kinship(proband, other_proband)
        # The kinship values may not be equal due to floating precision
        assert abs(proband_kinship - simple_algorithm_kinship) < 0.001

# Use the all pair-wise matrix to get the kinship among all the vertices in the pedigree
for vertex_1, vertex_2 in itertools.combinations(pedigree, 2):
    kinship = direct_kinship_matrix[vertex_to_index[vertex_1], vertex_to_index[vertex_2]]

# Convert the kinship matrix to a numpy matrix and free the previous matrix
# NOTE: The kinship_c_implementation object becomes empty after this
(vertex_to_id, numpy_matrix) = kinship_matrix_c.to_numpy_and_free()
for proband in pedigree.get_sink_vertices():
    for other_proband in pedigree.get_sink_vertices():
        numpy_matrix_kinship = numpy_matrix[vertex_to_id[proband], vertex_to_id[other_proband]]
        simple_algorithm_kinship = direct_kinship_matrix[vertex_to_index[proband], vertex_to_index[other_proband]]
        assert abs(simple_algorithm_kinship - numpy_matrix_kinship) < 0.001
        try:
            # The call must fail because the matrix is empty after the transformation
            kinship_matrix_c.get_kinship(proband, other_proband)
            assert False
        except IndexError:
            continue

print("Kinships have been calculated")
