import csv

from lineagekit.core.Pedigree import Pedigree
from lineagekit.utility.utility import *
import time

filepath = get_file_path("Specify the path to the pedigree:")

pedigree = Pedigree.get_pedigree_graph_from_file(filepath=filepath, separation_symbol=" ",
                                                 missing_parent_notation=["-1"], skip_first_line=True)
founders = set(pedigree.get_top_level_vertices())
founders_and_probands = founders.union(pedigree.get_sink_vertices())
start_time = time.time()
kinship_matrix_c = pedigree.calculate_probands_kinship(probands=founders_and_probands)
end_time = time.time()
elapsed_time = end_time - start_time
print(f"Function took {elapsed_time:.4f} seconds to execute.")

cartegene_filename = 'cartegene_probands_and_top_level_vertices_kinships.csv'
with open(cartegene_filename, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Proband_1_id', 'Proband_2_id', 'Kinship'])
    for proband in founders_and_probands:
        for other_proband in founders_and_probands:
            if proband >= other_proband:
                continue
            proband_kinship = kinship_matrix_c.get_kinship(proband, other_proband)
            writer.writerow([proband, other_proband, proband_kinship])
print("Kinships have been calculated")
