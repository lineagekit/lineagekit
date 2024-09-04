.. _ascending_genealogy:

#################################
Calculating ascending genealogy
#################################

All the classes extending from :class:`GenGraph` have functionality for calculating
the ascending graph for a given set of vertices. In general, it works like this:

.. code-block:: python

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

For instance, we can find this function useful when we work with the BALSAC dataset, but we're only
interested in the Cartegene individuals:

.. code-block:: python

    # Specify the path
    balsac_path = input("Specify the path to BALSAC: ")
    # Parse the BALSAC dataset
    balsac_pedigree = GenGraph.get_graph_from_file(filepath=balsac_path)
    # Get the list of Cartgene participants
    cartegene_participants = # Getting the cartegene individuals ... #
    # Remove all the vertices that are not present in the ascending pedigree of the Cartegene individuals
    balsac_pedigree.reduce_to_ascending_graph(probands=cartegene_participants)

Another example can involve a coalescent tree. Suppose that we want to sub-select
a given list of probands in the tree and work with the result. We can use the :meth:`reduce_to_ascending_graph <GenGraph.GenGraph.reduce_to_ascending_graph>`
function in :class:`GenGraph` and :meth:`remove_unary_nodes <CoalescentTree.CoalescentTree.remove_unary_nodes>` function in :class:`CoalescentTree`
to achieve this.

.. code-block:: python

    # Specify the path to the coalescent tree
    tree_path = input("Specify the path to the genealogy file: ")
    # Parse the coalescent tree
    coalescent_tree = CoalescentTree.get_coalescent_tree_from_file(filepath=tree_path)
    # Get the list of probands
    probands = [1, 2, 4, 11, 34, 156, 1044, 12346, 274668, 12005168]
    # Remove all the vertices that are not in the ascending tree
    coalescent_tree.reduce_to_ascending_graph(probands=probands)
    # After reduction, the tree may contain unary nodes, which are nodes with only a single child.
    # We remove these vertices here
    coalescent_tree.remove_unary_nodes()
