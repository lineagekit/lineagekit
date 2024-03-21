# lineagekit
A python library with genealogy methods. We will seek to integrate code or functionality or both from:

# ISGEN (Dominic's code on importance sampling)
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

# genlib (Marie-Hélène Roy-Gagnon's R pacakge. 

# Andrii's aligner

# Luke pedigree simulations

# sgkit
