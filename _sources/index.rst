##################################################
lineagekit: Python library with genealogy methods
##################################################


###################
Introduction
###################
This is the manual for **lineagekit**, a population genetic library that contains various useful tools for working
with genealogies. With lineagekit you can:

* Parse genealogies from files of different formats.
* Build custom genealogies from scratch.
* Modify the obtained genealogies.
* Obtain information by using utility function provided by the library.
* Calculate kinships in large genealogies.

##################
Installation
##################

1. Clone the `repository <https://github.com/lineagekit/lineagekit>`_ and ``cd`` into the root of the project.
2. Run ``pip install .``

##################
Development policy
##################

If you would like to contribute to the project, or if you have found a certain issue with the code,
refer to our
`development policy <https://github.com/lineagekit/lineagekit/tree/dev?tab=readme-ov-file#development-policy>`_
for more information.

##################
Publishing
##################

If you want to publish a new version of `lineagekit` to PyPI (so that it can be installed with `pip`), follow these steps:

1. Update the version number in ``setup.cfg`` located in the root directory of the project. For example, let's say that the new version of the library is 1.1.0.

2. Create a Git tag for the version you specified by running::

     git tag v1.1.0

3. Push the tag to the remote repository with::

     git push origin v1.1.0


Now, you should be able to install the package with pip.

####################
Contents
####################

.. toctree::
   :maxdepth: 2

   quickstart
   gen_graph
   abstract_pedigree
   pedigree
   ploid_pedigree
   coalescent_tree

   Use-cases <use_cases/index>
