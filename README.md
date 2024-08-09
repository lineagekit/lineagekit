

# lineagekit
A python library with genealogy methods. We will seek to integrate code or functionality or both from:

# [Documentation](#documentaion)
You can find the documentation for this project [here](https://lineagekit.github.io/lineagekit/)

## ISGEN (Dominic's code on importance sampling)
https://github.com/DomNelson/ISGen

# Earlier iterations of genealogy alignment project:

Useful data structures built on top of `networkx`:
- `Genealogical`: Generic data structure that includes implementation of useful functions on geneological objects (trees, pedigrees, etc.).
- `Pedigree`
    - Implements Lange's Kinship estimation algorithm with sparse matrices (may need to be optimized)
- `HaplotypeGraph`: Converts from a pedigree to a graph encompassing all potential inheritance paths.
    - Supports `autosomal`, `X`, `Y`, and `mitochondrial` inheritance paths.
    - Has a (non-tested) message passing algorithm for computing kinship.
- `Traversal` (corresponds to trees)

https://github.com/ivan-krukov/aligning-genealogies/tree/master/genealogy_aligner

## genlib (Marie-Hélène Roy-Gagnon's R package). 

## Andrii's aligner

## Luke pedigree simulations

## sgkit

# Installation

1. Clone the repository and `cd` into the root of the project.
2. Run `pip install .`

To use this package, you need to have a C++ compiler installed on your machine with support for at least C++17.
# [Development policy](#development-policy)

If you want to contribute to this project, follow these steps:

1. Create an issue with a description of what changes or improvements you would like to add.
2. Create a separate branch for the issue that you have created. Make sure that the branch source is the **dev** branch.
3. On the branch that you have created, introduce the changes or fixes that you want. 
**Important!** Don't forget to add tests for the code that you've added. If you were working on fixing a bug in the
existing code, try to create a test case for which the old version of the code didn't work properly, but the new version does.
4. Clone the repository, checkout into the newly created branch, install all the requirements mentioned in
**requirements.txt**.
5. Push your changes and create a pull request.

If you would like to request a certain change in the code, or if you have found a bug / unexpected behaviour,
you can just create an issue with the corresponding description. 

If you want to build the documentation locally, `cd` into the `docs` directory and run `make html`.


