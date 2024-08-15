.. _kinship:

######################
Calculating kinships
######################

In population genetics, the term "kinship" typically refers to the level of genetic relatedness between a pair of
individuals, and is represented as a number between 0 and 1. More on calculating kinships can be found
`here <https://academic.oup.com/bioinformatics/article/35/6/1002/5085372>`_.

There are several methods available in the :class:`Pedigree` class for kinship calculation.

----------------------------------
Calculating all-pairwise kinships
----------------------------------

You can calculate the kinship coefficients between all the pairs of individuals in your pedigree by running the
:meth:`calculate_kinship <Pedigree.Pedigree.calculate_kinship>` function. This function calculates the kinship
coefficients and returns a `numpy <https://numpy.org/>`_ matrix as the result.
Additionally, the function returns a dictionary that allows you to index the matrix:

.. code-block:: python

    # Parse the input file given by the specified file path
    pedigree = Pedigree.get_pedigree_graph_from_file(filepath=small_pedigree_path)
    # Calculate the kinship coefficients
    (vertex_to_index, direct_kinship_matrix) = pedigree.calculate_kinship()
    for vertex_1, vertex_2 in itertools.combinations(pedigree, 2):
        kinship = direct_kinship_matrix[vertex_to_index[vertex_1], vertex_to_index[vertex_2]]


Every kinship value is stored as a 16-bit float value.

Note that the complexity of this function grows quadratically with the size of the pedigree. If your pedigree is large,
you may run out of memory.

----------------------------------
Proband kinship calculation
----------------------------------
In some cases, we only need to calculate the kinship coefficients among a given list of individuals, not the entire
pedigree. If the pedigree is small, we can still use the previous method. However, when the pedigree is large,
this approach may not be sufficient. For this purpose, we can use the
:meth:`calculate_probands_kinship <Pedigree.Pedigree.calculate_probands_kinship>` function, which is a much more
efficient version of the previous algorithm when we are only interested in a few kinship coefficient values.

==================================
API
==================================

The :meth:`calculate_probands_kinship <Pedigree.Pedigree.calculate_probands_kinship>` function takes two additional
parameters.
The first parameter, :attr:`probands`, is the list of individuals for which the kinship coefficients are calculated.
Note that there are no restrictions on the individuals in the list; they simply should be valid vertices in the pedigree.
If not specified, the sink vertices are used as probands.

The second parameter, :attr:`mode`, specifies the running mode of the algorithm. The default option, :attr:`SPEED`,
offers the best running time and should be used when memory is not an issue. In cases where memory restrictions are high,
the :attr:`MEMORY` option can be used to reduce memory usage at the cost of increased running time.
On average, the MEMORY option can save up to 25% of memory and take around 40-45% more running time.

.. code-block:: python

    # Calculating proband kinship coefficients
    kinship_matrix_default = pedigree.calculate_probands_kinship()
    kinship_matrix_custom = pedigree.calculate_probands_kinship(probands={1, 3}, mode=KinshipMode.MEMORY)

==================================
Algorithm description
==================================

The main algorithm's idea is described in section 4.2 of the `mentioned paper <https://academic.oup.com/bioinformatics/article/35/6/1002/5085372>`_. The general idea of the algorithm is as follows:

1. Form a queue of vertices for which the kinship coefficients should be calculated.
2. While the queue is not empty, do the following:

   a. Calculate the kinship coefficient between the chosen vertex and every other vertex in the matrix.

   b. Check whether there are parent vertices of this vertex such that all their children are in the matrix. If there are such non-proband parents, remove them from the kinship matrix.

   c. If there are any children of this vertex such that all their parents are in the matrix, add them to the queue.

This algorithm description is not the most precise, but it explains the general idea of how this function works.
