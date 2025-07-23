from lineagekit.core.GenGraph import GenGraph
from lineagekit.utility.utility import get_file_path, get_non_existing_path


def run_interactive_session():
    # Ask the user to specify the path to the pedigree
    genealogy_path = get_file_path("Specify the path to the genealogy file: ")
    result_filepath = get_non_existing_path("Specify the path to the result file: ")
    # Parse the file
    genealogy_graph = GenGraph.get_graph_from_file(filepath=genealogy_path, missing_parent_notation=['0'])
    # Specify the list of probands
    probands = [1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144]
    # Calculate the vertices in the ascending graph. This operation does not modify the graph
    vertices_in_ascending_graph = genealogy_graph.get_ascending_vertices_from_probands(probands=probands)
    # Remove all the vertices that are not in the ascending graph
    genealogy_graph.reduce_to_ascending_graph(probands=probands)
    genealogy_graph.save_to_file(filepath=result_filepath, missing_parent_notation='0')


if __name__ == '__main__':
    run_interactive_session()
