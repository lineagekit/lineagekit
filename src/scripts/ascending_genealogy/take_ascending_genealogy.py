from lineagekit.core.GenGraph import GenGraph

# Ask the user to specify the path to the pedigree
genealogy_path = input("Specify the path to the genealogy file: ")
# Parse the file
genealogy_graph = GenGraph.get_graph_from_file(filepath=genealogy_path)
# Specify the list of probands
probands = [1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144]
# Calculate the vertices in the ascending graph. This operation does not modify the graph
vertices_in_ascending_graph = genealogy_graph.get_ascending_vertices_from_probands(probands=probands)
# Remove all the vertices that are not in the ascending graph
genealogy_graph.reduce_to_ascending_graph(probands=probands)
