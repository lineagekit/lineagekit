import os
from pathlib import Path

from lineagekit.core.Pedigree import Pedigree
from lineagekit.utility.utility import get_file_path, get_non_existing_path

results_dir_name = "results"

filepath = get_file_path("Specify the path to the pedigree:")
results_filename = get_non_existing_path("Specify the filename to save the results (without extension):")
pedigree = Pedigree.get_pedigree_graph_from_file(filepath=filepath, separation_symbol=" ",
                                                 missing_parent_notation=["-1"], skip_first_line=True)
founders_contribution_factors = pedigree.calculate_contribution_factors(vertices=pedigree.nodes)
results_filepath = Path(results_dir_name) / f"{results_filename}.csv"
os.makedirs(results_dir_name, exist_ok=True)
with open(results_filepath, "a") as results_file:
    results_file.write("vertex_id,contribution_factor\n")
    for founder, factor in founders_contribution_factors.items():
        results_file.write(f"{founder},{factor}\n")
print("The contribution factors have been calculated")
