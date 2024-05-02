.. _quickstart:

==================
Quickstart
==================

------------------
Parsing files
------------------

One of the most common use-cases that we have when working with genealogies is parsing pedigree data from a file.
For example, let's assume that we need to parse the Balsac dataset and perform some operations with it.

To do that, we can use the **Pedigree** class as follows:

.. code-block:: python

    # Parses the input file given by the specified file path
    balsac = = Pedigree.get_pedigree_graph_from_file(filepath=filepath, separation_symbol=";",
                                                     missing_parent_notation=[""], skip_first_line=True)
    # Performing some operations on the obtained pedigree


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

For information about other parameters, refer to the documentation of the graph class that you are using.

------------------------------------
Working with genealogies
------------------------------------
.. _NetworkX: https://networkx.org/documentation/stable/index.html

Once you have created the genealogy, you can start working with it using the various methods available by the library.
The base functionality for all the graphs is located within the **GenGraph** class from which all the other graph
classes extend. **GenGraph**, in turn, extends from the
`DiGraph <https://networkx.org/documentation/stable/reference/classes/digraph.html#networkx.DiGraph>`_ class defined in
`NetworkX <https://networkx.org/documentation/stable/index.html>`_.

More information about various functions ...

------------------
Saving the results
------------------

More information regarding saving the results ...