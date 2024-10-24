from __future__ import annotations

import random
import warnings
from typing import Iterable

from basic.AbstractPedigree import AbstractPedigree

from src.utility.utility import random_subselect_poisson


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
            The number of vertex's parents' spouses.
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

    def _introduce_errors(self, error_rate: float):
        """
        Simulates the errors within the pedigree and returns them.
        The errors are simulated using the following approach
        1) A random number of error events is simulated based on the error rate. Specifically, a random Poisson variate
           is taken where the math expectation is n * error_rate where n is the number of individuals in
           the graph with parents
        2) For every selected individual v, a random parent w is selected and a random individual u from
           the same level as w. Then, v is disconnected from w and is connected with u.
           Notice that the actual graph is not changed, this changes will take place once they are applied
        Below you can find the illustration of these changes:
            w_1  #Same level#   w_2
             \
              \
               \
                \
                 \
                  \
                   \
                   v_p
        After applying the error:
            w_1  #Same level#  w_2
                            /
                           /
                          /
                         /
                        /
                       /
                      /
                   v_p
        Args:
            error_rate: The probability of an error per individual
        Returns:
            The function returns a list of tuples where every tuple represents an error.
            The tuple structure is as follows:
            1) The first value - the child id for which the parent has been changed.
            2) The second tuple contains the old parent id. In other words, this is the edge that needs to be removed
            3) The third tuple contains the new parent id. In other words, this is the edge that needs to be added

            The nesting of the ids into separate tuples may seem redundant, but this formate is used to make it
            compatible with the general API.
        """
        vertices_with_parents = [x for x in self if self.has_parents(x)]
        random_individuals = random_subselect_poisson(vertices_with_parents, error_rate)
        errors_result_list = []
        for vertex in random_individuals:
            vertex_parents = self.get_parents(vertex)
            random_parent = random.sample(vertex_parents, 1)[0]
            random_parent_level = self.get_vertex_level(random_parent)
            try:
                new_parent = self._select_new_parent_from_level(random_parent_level, vertex_parents)
            except ValueError:
                warnings.warn("One of the errors has been skipped as there were no other "
                              "candidates for reconnection")
                continue
            errors = (vertex, (random_parent,), (new_parent,))
            errors_result_list.append(errors)
        return errors_result_list
