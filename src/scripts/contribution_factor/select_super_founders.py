import csv
import os

from utility.utility import get_file_path, get_float_value, get_non_existing_path, \
    get_natural_number_input_in_bounds, get_natural_number_input

results_dir_name = "results"


def filter_by_threshold(reader, threshold):
    """Filters vertices based on a contribution factor threshold"""
    top_founders = []
    for row in reader:
        try:
            vertex_id, contribution_factor = row[0], float(row[1])
            if contribution_factor >= threshold:
                top_founders.append((vertex_id, contribution_factor))
        except (IndexError, ValueError):
            print(f"Skipping invalid row: {row}")
    return top_founders


def filter_top_n(reader, top_n):
    """Selects the top N vertices with the highest contribution factors"""
    try:
        return sorted(
            [(row[0], float(row[1])) for row in reader],
            key=lambda x: x[1],
            reverse=True
        )[:top_n]
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


def write_results(result_path, top_founders):
    """Writes the selected vertices to a CSV file"""
    with open(result_path, "w", newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["vertex_id", "contribution_factor"])
        writer.writerows(top_founders)
    print(f"Filtered results saved to {result_path}")


def process_input_file(reader, mode) -> [int]:
    if mode == 1:
        threshold_input_prompt = "Enter the contribution factor threshold: "
        threshold = get_float_value(threshold_input_prompt)
        return filter_by_threshold(reader, threshold)
    if mode == 2:
        victors_number_input_prompt = "Enter the number of top vertices to select: "
        top_n = get_natural_number_input(victors_number_input_prompt)
        return filter_top_n(reader, top_n)


def select_top_founders():
    file_path = get_file_path("Enter the path to the data file: ")
    result_filename = get_non_existing_path("Specify the name for the resulting file (without extension):")
    result_path = os.path.join(results_dir_name, f"{result_filename}.csv")
    input_prompt = "Choose mode:\n(1) Threshold-based filtering\n(2) Select top N\n"
    mode = get_natural_number_input_in_bounds(input_request=input_prompt, lower_bound=1, upper_bound=2)
    try:
        with open(file_path, newline='') as csvfile:
            reader = csv.reader(csvfile)
            next(reader, None)
            top_founders = process_input_file(reader, mode)
        write_results(result_path, top_founders)
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' does not exist.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


if __name__ == "__main__":
    select_top_founders()
