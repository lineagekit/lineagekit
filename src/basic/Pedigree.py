from __future__ import annotations

from typing import Iterable

import numpy as np

from src.basic.GenGraph import GenGraph

from enum import Enum
import kinship


class KinshipMode(Enum):
    SPEED = 0
    MEMORY = 1


class Pedigree(GenGraph):

    def __init__(self, *args, **kwargs):
        super().__init__(parent_number=2, *args, **kwargs)

    @staticmethod
    def get_pedigree_graph_from_file(filepath: str, probands: Iterable[int] = None,
                                     missing_parent_notation=None, separation_symbol=' ',
                                     skip_first_line: bool = False) -> Pedigree:
        """!
        @brief This method processes the input graph and builds the corresponding pedigree. There is one vertex
        for an individual in the input file.
        """
        gen_graph = GenGraph.get_graph_from_file(filepath=filepath,
                                                 probands=probands,
                                                 parent_number=2,
                                                 missing_parent_notation=missing_parent_notation,
                                                 separation_symbol=separation_symbol,
                                                 skip_first_line=skip_first_line)
        result = Pedigree()
        result.update(gen_graph)
        return result

    def get_vertex_spouses(self, vertex: int):
        spouses = {parent for child in self.get_children(vertex) for parent in self.get_parents(child)}
        spouses.discard(vertex)
        return spouses

    def get_vertex_parents_number_of_spouses(self, vertex: int) -> int:
        spouses = set()
        for parent in self.get_parents(vertex):
            spouses.update(self.get_vertex_spouses(parent))
        return len(spouses)

    def get_vertices_parents(self, vertices: Iterable[int]):
        parents = set()
        for vertex in vertices:
            parents.update(self.get_parents(vertex))
        return parents

    def get_vertices_children(self, vertices: Iterable[int]):
        children = set()
        for vertex in vertices:
            children.update(self.get_children(vertex))
        return children

    def _calculate_self_kinship(self, kinship_matrix: np.ndarray, vertex: int, vertex_to_index: {int: int}):
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
        if probands is None:
            probands = frozenset(self.get_sink_vertices())
            children_map = {x: self.get_children(x) for x in self}
            parents_map = {x: self.get_parents(x) for x in self}
        else:
            # Calculate the ascending genealogy if the proband list is custom
            ascending_genealogy: set[int] = self.get_ascending_graph_from_vertices(probands)
            children_map = {x: list(ascending_genealogy.intersection(self.get_children(x))) for x in ascending_genealogy}
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
