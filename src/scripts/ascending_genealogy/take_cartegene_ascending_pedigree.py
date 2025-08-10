from lineagekit.core.gen_graph import GenGraph
import csv

individual_column = 'ind'
cartegene_column = 'ego_cag'

# Ask the user to specify the path to the pedigree
genealogy_path = input("Specify the path to the Balsac dataset: ")

# Parse the information about the cartegene participants
probands = []
with open(genealogy_path, mode='r') as file:
    reader = csv.DictReader(file, delimiter=';')
    for row in reader:
        if row[cartegene_column]:
            probands.append(int(row[individual_column]))

# Parse the file
genealogy_graph = GenGraph.get_graph_from_file(filepath=genealogy_path, separation_symbol=";",
                                               skip_first_line=True, missing_parent_notation=[""])

# Remove all the vertices that are not in the ascending graph
genealogy_graph.reduce_to_ascending_graph(probands=probands)
genealogy_graph.save_to_file("cartegene_ascending.pedigree")
