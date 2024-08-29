from __future__ import annotations

from typing import Iterable

from basic.AbstractPedigree import AbstractPedigree


class Pedigree(AbstractPedigree):
    """
    The general pedigree class where every node represents an individual.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @staticmethod
    def get_pedigree_graph_from_file(filepath: str, probands: Iterable[int] = None,
                                     missing_parent_notation=None, separation_symbol=' ',
                                     skip_first_line: bool = False) -> Pedigree:
        """
        Processes the input graph and builds the corresponding pedigree. There is one vertex
        for an individual in the input file.

        Args:
            filepath (str): The path to the input file.
            probands (Iterable[int]): Optional parameter. The probands for which the ascending genealogy should be
                                      calculated. By default, all the vertices from the input file are stored.
            missing_parent_notation: The list of text sequences representing that the given individual has no parents.
                                     If not specified, the default values "-1" and "." are used (meaning that both are
                                     accepted at the same time).
            separation_symbol (str): The symbol used to separate the values in a line. By default, a space is used.
            skip_first_line (bool): Specifies whether the first line in the file should be skipped. Can be useful if the
                                    header does not start with a '#' symbol.
        """
        gen_graph = AbstractPedigree.get_graph_from_file(filepath=filepath,
                                                         probands=probands,
                                                         parent_number=2,
                                                         missing_parent_notation=missing_parent_notation,
                                                         separation_symbol=separation_symbol,
                                                         skip_first_line=skip_first_line)
        result = Pedigree()
        result.update(gen_graph)
        return result

    def get_vertex_spouses(self, vertex: int):
        """
        Args:
            vertex: The vertex id.

        Returns:
            The vertex's spouses.
        """
        spouses = {parent for child in self.get_children(vertex) for parent in self.get_parents(child)}
        spouses.discard(vertex)
        return spouses

    def get_vertex_parents_number_of_spouses(self, vertex: int) -> int:
        """
        Args:
            vertex: The vertex id.

        Returns:
            The number of vertex's parents' spouses..
        """
        spouses = set()
        for parent in self.get_parents(vertex):
            spouses.update(self.get_vertex_spouses(parent))
        return len(spouses)

    def get_vertices_parents(self, vertices: Iterable[int]):
        """
        Args:
            vertices:  The vertices ids.

        Returns:
            All the parents for the given vertices. That is, a list of vertices where each vertex in the list is
            a child of one of the specified vertices.
        """
        parents = set()
        for vertex in vertices:
            parents.update(self.get_parents(vertex))
        return parents

    def get_vertices_children(self, vertices: Iterable[int]):
        """
        Args:
            vertices: The vertices ids.

        Returns:
            All the children for the given vertices. That is, a list of vertices where each vertex in the list is
            a child of one of the specified vertices.
        """
        children = set()
        for vertex in vertices:
            children.update(self.get_children(vertex))
        return children
