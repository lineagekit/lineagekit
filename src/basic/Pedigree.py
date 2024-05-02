from __future__ import annotations

from collections import defaultdict
from typing import Iterable
from src.basic.GenGraph import GenGraph
import numpy as np


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

    def _calculate_self_kinship_sparse(self, sparse_kinship_matrix: dict, vertex: int):
        vertex_parents = self.get_parents(vertex)
        if len(vertex_parents) != 2:
            kinship = 0.5
        else:
            [first_parent, second_parent] = vertex_parents
            kinship = (1 + sparse_kinship_matrix[first_parent][second_parent]) / 2
        sparse_kinship_matrix[vertex][vertex] = kinship

    def _calculate_pair_kinship_sparse(self, sparse_kinship_matrix: dict, first_vertex: int, second_vertex: int):
        first_second_kinship_non_normalized = 0
        for vertex_parent in self.get_parents(first_vertex):
            first_second_kinship_non_normalized += sparse_kinship_matrix[vertex_parent][second_vertex]
        first_second_kinship = first_second_kinship_non_normalized / 2
        sparse_kinship_matrix[first_vertex][second_vertex] = first_second_kinship
        sparse_kinship_matrix[second_vertex][first_vertex] = first_second_kinship
        # if self.vertex_to_level_map[first_vertex] == self.vertex_to_level_map[second_vertex]:
        #     sparse_kinship_matrix[first_vertex][second_vertex] = first_second_kinship
        #     sparse_kinship_matrix[second_vertex][first_vertex] = first_second_kinship
        # else:
        #     higher_generation_vertex = (first_vertex if self.vertex_to_level_map[first_vertex] >
        #                                                 self.vertex_to_level_map[second_vertex] else second_vertex)
        #     other_vertex = second_vertex if higher_generation_vertex != second_vertex else first_vertex
        #     sparse_kinship_matrix[higher_generation_vertex][other_vertex] = first_second_kinship

    def get_min_levels(self):
        min_levels = []
        current_level = self.get_levels()[0]
        while current_level:
            next_level = set()
            [next_level.update(self.get_parents(vertex)) for vertex in current_level]
            current_level = sorted(current_level, key=lambda vertex: self.vertex_to_level_map[vertex], reverse=True)
            min_levels.append(current_level)
            for vertex in current_level:
                if vertex in next_level:
                    next_level.remove(vertex)
            current_level = list(next_level)
        return min_levels

    def calculate_kinship_for_probands_sparse(self):
        kinship_sparse_matrix = defaultdict(dict)
        parent_to_remaining_children = {x: set(self.get_children(x)) for x in self}
        reversed_levels = reversed(self.get_levels())
        for level in reversed_levels:
            processed_vertices = list(kinship_sparse_matrix.keys())
            for vertex in level:
                processed_vertices.append(vertex)
                for processed_vertex in processed_vertices:

                    if processed_vertex == vertex:
                        self._calculate_self_kinship_sparse(sparse_kinship_matrix=kinship_sparse_matrix, vertex=vertex)
                    elif processed_vertex not in kinship_sparse_matrix:
                        continue
                    else:
                        self._calculate_pair_kinship_sparse(sparse_kinship_matrix=kinship_sparse_matrix,
                                                            first_vertex=vertex,
                                                            second_vertex=processed_vertex)
                for parent in self.get_parents(vertex):
                    parent_to_remaining_children[parent].remove(vertex)

                    if not parent_to_remaining_children[parent]:
                        kinship_sparse_matrix.pop(parent)
                        for other_vertex in kinship_sparse_matrix:
                            kinship_sparse_matrix[other_vertex].pop(parent)
        return kinship_sparse_matrix

    def calculate_kinship_for_probands_cutting(self):
        pass
