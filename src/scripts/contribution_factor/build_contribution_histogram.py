import os

import matplotlib
import matplotlib.pyplot as plt
import csv

from utility.utility import get_file_path, get_non_existing_path

results_dir = "results"
matplotlib.use("TkAgg")


def build_histogram():
    file_path = get_file_path("Specify the path to the data file:")
    result_file = get_non_existing_path("Specify the name of the resulting file (without extension):")
    contribution_factors = []

    with open(file_path, newline='') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        for row in reader:
            contribution_factors.append(float(row[1]))

    # Plot histogram
    plt.figure(figsize=(8, 5))
    plt.hist(contribution_factors, bins=10, edgecolor='black', alpha=0.7,
             log=True, label='Contribution factor')
    # Labels and title
    plt.xlabel("Contribution Factor")
    plt.ylabel("Frequency")
    plt.title("Histogram of Contribution Factors")
    result_filename = f"{result_file}.svg"
    plt.savefig(result_filename, dpi=300)  # High-resolution output
    print(f"Histogram saved as {result_filename}")
    plt.show()


os.chdir(results_dir)
build_histogram()
