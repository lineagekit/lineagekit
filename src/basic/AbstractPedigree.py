from enum import Enum
from typing import Iterable

import kinship
import numpy as np

from basic.GenGraph import GenGraph

from abc import ABC, abstractmethod
import random


class KinshipMode(Enum):
    SPEED = 0
    MEMORY = 1


class AbstractPedigree(GenGraph, ABC):
    """
    The abstract base class for all Pedigree abstract classes containing all the common functionality.
    This class is for internal use only.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(parent_number=2, *args, **kwargs)

    def _calculate_self_kinship(self, kinship_matrix: np.ndarray, vertex: int, vertex_to_index: {int: int}):
        """
        Helper function that calculates the self-kinship.

        Args:
            kinship_matrix: The kinship matrix.
            vertex: The vertex id.
            vertex_to_index: The dictionary that maps vertex ids to indices.

        Returns:
            The vertex self-kinship.
        """
        vertex_parents = self.get_parents(vertex)
        vertex_index = vertex_to_index[vertex]
        if len(vertex_parents) != 2:
            kinship_matrix[vertex_index, vertex_index] = 0.5
        else:
            vertex_parents_indices = [vertex_to_index[vertex_parent] for vertex_parent in vertex_parents]
            kinship_matrix[vertex_index, vertex_index] = (
                    (1 + kinship_matrix[vertex_parents_indices[0], vertex_parents_indices[1]]) / 2)

    def _calculate_pair_kinship(self, kinship_matrix: np.ndarray, first_vertex: int,
                                second_vertex: int, vertex_to_index: {int: int}):
        """
        Helper function that calculates the kinship between two vertices.

        Args:
            kinship_matrix: The kinship matrix.
            first_vertex: The first vertex id.
            second_vertex: The second vertex id.
            vertex_to_index: The dictionary that maps vertex ids to indices.

        Returns:
            The kinship between the two vertices.
        """
        first_vertex_index = vertex_to_index[first_vertex]
        second_vertex_index = vertex_to_index[second_vertex]
        vertex_parents_indices = [vertex_to_index[vertex_parent] for vertex_parent in self.get_parents(first_vertex)]
        first_second_kinship_non_normalized = 0
        for vertex_parent_index in vertex_parents_indices:
            first_second_kinship_non_normalized += kinship_matrix[vertex_parent_index, second_vertex_index]
        first_second_kinship = first_second_kinship_non_normalized / 2
        kinship_matrix[first_vertex_index, second_vertex_index] = first_second_kinship
        kinship_matrix[second_vertex_index, first_vertex_index] = first_second_kinship

    def calculate_kinship(self):
        """
        Calculates all-pairwise kinship coefficients.

        Returns:
            tuple:
                A tuple containing:
                    1. A dictionary mapping vertex IDs to their corresponding kinship matrix indices.
                    2. A numpy.ndarray representing the kinship matrix where each element (i, j) is the
                      kinship coefficient between vertex i and vertex j.

        Example:
            >>> # Example usage of the function:
            >>> vertex_to_index, direct_kinship_matrix = pedigree.calculate_kinship()
            >>> vertices_kinship = direct_kinship_matrix[vertex_to_index[vertex_1], vertex_to_index[vertex_2]]

        """
        n = self.get_vertices_number()
        kinship_matrix = np.empty((n, n), dtype=np.float16)
        vertex_to_index = dict()
        current_index = 0
        for level in reversed(self.get_levels()):
            for vertex in level:
                vertex_to_index[vertex] = current_index
                for processed_vertex, index in vertex_to_index.items():
                    if processed_vertex == vertex:
                        self._calculate_self_kinship(kinship_matrix=kinship_matrix, vertex=vertex,
                                                     vertex_to_index=vertex_to_index)
                    else:
                        self._calculate_pair_kinship(kinship_matrix=kinship_matrix, first_vertex=vertex,
                                                     second_vertex=processed_vertex, vertex_to_index=vertex_to_index)
                current_index += 1
        return vertex_to_index, kinship_matrix

    def calculate_probands_kinship(self, probands: set[int] = None, mode: KinshipMode = KinshipMode.SPEED):
        """
        Calculates the all-pairwise kinship coefficients for the given list of vertices.

        Args:
            probands (set[int]): The list of vertices for which the kinship coefficients should be calculated.
                                 If not specified, the sink vertices are used.
            mode (KinshipMode): The running mode of the kinship calculation. By default, the SPEED option is used,
                                which offers the best performance in terms of running time and is suitable when memory
                                is not a concern.
                                Alternatively, the MEMORY option is available for cases where memory usage is a concern.
                                On average, the MEMORY option reduces memory usage by approximately 25% but results in
                                about 40% additional running time.

        Returns:
            The kinship matrix for the probands.

        Example:
            >>> kinship_matrix = pedigree.calculate_probands_kinship()
            >>> probands_kinship = kinship_matrix.get_kinship(proband, other_proband)
        """
        if probands is None:
            probands = frozenset(self.get_sink_vertices())
            children_map = {x: self.get_children(x) for x in self}
            parents_map = {x: self.get_parents(x) for x in self}
        else:
            # Calculate the ascending genealogy if the proband list is custom
            ascending_genealogy: set[int] = self.get_ascending_vertices_from_probands(probands)
            children_map = {x: list(ascending_genealogy.intersection(self.get_children(x))) for x in
                            ascending_genealogy}
            parents_map = {x: self.get_parents(x) for x in ascending_genealogy}
        if mode == KinshipMode.SPEED:
            kinship_sparse_matrix = kinship.calculate_kinship_sparse_speed(
                sink_vertices=probands, children=children_map,
                parents=parents_map
            )
        else:
            kinship_sparse_matrix = kinship.calculate_kinship_sparse_memory(
                sink_vertices=probands, children=children_map,
                parents=parents_map
            )
        return kinship_sparse_matrix

    def _select_new_parent_from_level(self, level_index: int, vertex_parents: Iterable[int]):
        # We need to select a vertex from the level_index level.
        # The obvious solution would be to take all the vertices from that level and remove all the parent vertices.
        # Unfortunately, this approach can be quite slow for large graphs, as we need to filter a big list on
        # every iteration of this loop.
        # We could also try to select len(vertex_parents) + 1 vertices from the level and filter
        # the resulting list, but there are corner cases as well that make the code to be not that elegant.
        # So here, we simply try our luck and select a random vertex from the same level. This usually works
        # on the first try and performs fast. If we fail a large number of times, we reserve to the initial version
        # that was mentioned. It's extremely likely that in that case we don't have any other vertices from that
        # level, so the length of the list will be <= 4 elements.
        max_attempts = len(self.get_levels()[level_index]) * 2
        attempts = 0
        while attempts < max_attempts:
            new_parent = random.choice(self.get_levels()[level_index])
            if new_parent not in vertex_parents:
                break
            attempts += 1
        else:
            new_parent_candidates = [x for x in self.get_levels()[level_index] if x
                                     not in vertex_parents]
            if not new_parent_candidates:
                raise ValueError("No other candidates available")
            new_parent = random.sample(new_parent_candidates, 1)[0]
        return new_parent

    @abstractmethod
    def _introduce_errors(self, error_rate: float):
        pass

    def apply_errors(self, errors):
        """
        Applies the errors simulated by the :meth:`introduce_and_record_errors` function.
        """
        for vertex, remove_parents, add_parents in errors:
            self.remove_edges_from(((remove_parent, vertex) for remove_parent in remove_parents))
            self.add_edges_from(((add_parent, vertex) for add_parent in add_parents))

    def reverse_errors(self, errors):
        """
        Reverses the errors that have been applied by the :meth:`apply_errors` function.
        """
        for vertex, remove_parents, add_parents in errors:
            self.remove_edges_from(((add_parent, vertex) for add_parent in add_parents))
            self.add_edges_from(((remove_parent, vertex) for remove_parent in remove_parents))

    def introduce_and_record_errors(self, error_rate: float, apply_errors: bool = True):
        """
        Simulates the errors within the pedigree and returns them. Additionally, the errors can be automatically
        applied if ``apply_errors`` is set to ``True``. The errors can be also applied manually by calling the
        :meth:`apply_errors` function.

        Args:
            error_rate: The probability of an individual being connected to a different parent
            apply_errors: If set to ``True``, the errors will be applied automatically. Default is ``True``.

        Returns:
            The function returns a list of tuples where every tuple represents an error.
            The tuple structure is as follows:
            1) The first value - the child id for which the parents have been changed.
            2) The second tuple contains the old parent ids. In other words, these are the edges that
            need to be removed
            3) The third tuple contains the new parent ids. In other words, these are the edges that
            need to be added
        """
        errors = self._introduce_errors(error_rate=error_rate)
        if apply_errors:
            self.apply_errors(errors)
        return errors
