import csv

from basic.Pedigree import Pedigree
from utility.utility import *
import time

filepath = get_file_path("Specify the path to the pedigree:")

pedigree = Pedigree.get_pedigree_graph_from_file(filepath=filepath, separation_symbol=" ",
                                                 missing_parent_notation=["-1"], skip_first_line=True)
start_time = time.time()
kinship_matrix_c = pedigree.calculate_probands_kinship()
end_time = time.time()
elapsed_time = end_time - start_time
print(f"Function took {elapsed_time:.4f} seconds to execute.")

cartegene_filename = 'filtered_individuals.csv'
with open(cartegene_filename, mode='w', newline='') as file:
    writer = csv.writer(file)

    writer.writerow(['Proband_1_id', 'Proband_2_id', 'Kinship'])

    for proband in pedigree.get_sink_vertices():
        for other_proband in pedigree.get_sink_vertices():
            if proband == other_proband:
                continue
            proband_kinship = kinship_matrix_c.get_kinship(proband, other_proband)
            writer.writerow([proband, other_proband, proband_kinship])

print("Kinships have been calculated")
