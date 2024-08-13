from basic.CoalescentTree import CoalescentTree
from basic.PloidPedigree import PloidPedigree
from utility.utility import *


def save_largest_clade_and_get_probands():
    filepath = get_file_path("Specify the path to the coalescent tree:\n")
    coalescent_tree: CoalescentTree = CoalescentTree.get_coalescent_tree_from_file(filepath=filepath)
    largest_clade = coalescent_tree.get_largest_clade_by_probands()
    coalescent_tree.reduce_to_ascending_graph(largest_clade)
    while True:
        try:
            result_clade_filepath = input("Specify the path to the resulting file:\n")
            coalescent_tree.save_to_file(result_clade_filepath)
            break
        except OSError:
            print("Specify a different file path")
    return set(coalescent_tree.get_sink_vertices())


probands = save_largest_clade_and_get_probands()
pedigree_filepath = get_file_path("Specify the path to the pedigree file:\n")
print("Processing the graph")
probands = PloidPedigree.get_ploids_from_individual_ids(PloidPedigree.get_individual_ids_from_ploids(probands))
genealogical_graph = PloidPedigree.get_ploid_pedigree_from_file(filepath=pedigree_filepath, probands=probands)
result_filepath = input("Specify the path to the resulting file for the ascending genealogy:\n")
print("Saving the ascending genealogy")
genealogical_graph.save_as_diploid(filepath=result_filepath)
