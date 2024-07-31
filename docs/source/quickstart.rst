.. _quickstart:

##################
Quickstart
##################

------------------
Parsing files
------------------

One of the most common use-cases that we have when working with genealogies is parsing pedigree data from a file.
For example, let's assume that we need to parse the Balsac dataset and perform some operations with it.

For example, we can use the **Pedigree** class to do this:

.. code-block:: python

    # Parse the input file given by the specified file path
    balsac = Pedigree.get_pedigree_graph_from_file(filepath=filepath, separation_symbol=";",
                                                   missing_parent_notation=[""], skip_first_line=True)
    # Perform some operations on the obtained pedigree


This code snippet parses the Balsac dataset specified by the given file path and returns an instance of the
**Pedigree** class with which we can work.

The *separation_symbol* parameter specifies which sequence of characters is
used to separate the columns in the input file.

The *missing_parent_notation* variable specifies which sequence of
characters is used to specify that the given parent is unknown.

In most cases, either "." or "-1" is used for this purpose. However, in Balsac, missing parents are specified by
an empty string. For example:

    `1;;; #Other parameters`

means that the individual 1 has no parents.

The *skip_first_line* option tells the function whether the first line of the input file should be ignored. In most
cases, the first line contains some comments about the file that must not be treated as a part of the actual data.

Currently, **GenGraph**, **Pedigree**, **PloidPedigree** and **CoalescentTree** are the available graph classes. You
can refer to the class' documentation to get more information.

------------------------------------
Working with genealogies
------------------------------------
.. _NetworkX: https://networkx.org/documentation/stable/index.html

Once you have created the genealogy, you can start working with it using the various methods available by the library.
The base functionality for all the graphs is located within the **GenGraph** class from which all the other
classes extend. **GenGraph**, in turn, extends from the
`DiGraph <https://networkx.org/documentation/stable/reference/classes/digraph.html#networkx.DiGraph>`_ class defined in
`NetworkX <https://networkx.org/documentation/stable/index.html>`_.

The basic functionality that is available in every class includes adding/removing vertices and edges,
calculating the ascending genealogy, and serializing the obtained graph.

.. code-block:: python

    # Parse the input file given by the specified file path
    balsac = Pedigree.get_pedigree_graph_from_file(filepath=filepath, separation_symbol=";",
                                                   missing_parent_notation=[""], skip_first_line=True)
    # Get the number of vertices in the pedigree
    individuals_number = balsac.get_vertices_number()
    # Check if the given vertex is a founder
    balsac.is_founder(1)
    # Remove the vertex
    balsac.remove_node(1)
    # Remove all the vertices that are not present in the ascending genealogy
    balsac.reduce_to_ascending_genealogy([1, 2, 3, 5, 8, 13, 21, 34])
    # Save the resulting graph to a file
    balsac.save_to_file("balsac_ascending.txt")

You should refer to the corresponding class documentation to read the full function documentation.
