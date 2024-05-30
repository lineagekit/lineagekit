from __future__ import annotations

import heapq
import itertools
import math
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from typing import Iterable

import numpy as np

from src.basic.GenGraph import GenGraph
import kinship

numpy_one_half = np.float16(0.5)
numpy_zero = np.float16(0)


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
        kinship_matrix = np.empty((n, n), dtype=np.float64)
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

    def _calculate_self_kinship_sparse(self, kinship_sparse_matrix: dict, vertex: int):
        # vertex_parents = self.get_parents(vertex)
        vertex_parents = self._pred[vertex]
        if len(vertex_parents) != 2:
            kinship = numpy_one_half
        else:
            first_parent, second_parent = tuple(vertex_parents)
            kinship = np.float16((1 + kinship_sparse_matrix[first_parent][second_parent]) / 2)
        kinship_sparse_matrix[vertex][vertex] = kinship

    def _calculate_half_kinship(self, kinship_matrix: dict, processed_vertex: int, half_vertex: int,
                                half_processed_vertices: {int}):
        processed_parent = [x for x in self.get_parents(half_vertex) if x in kinship_matrix
                            and x not in half_processed_vertices][0]
        assert len([x for x in self.get_parents(half_vertex) if x in kinship_matrix
                    and x not in half_processed_vertices]) == 1
        half_kinship = kinship_matrix[processed_parent][processed_vertex] / 2
        kinship_matrix[half_vertex][processed_vertex] = half_kinship
        kinship_matrix[processed_vertex][half_vertex] = half_kinship

    def _calculate_pair_kinship_sparse(self, kinship_sparse_matrix: dict, first_vertex: int, second_vertex: int):
        first_second_kinship_non_normalized = 0
        for vertex_parent in self._pred[first_vertex]:
            first_second_kinship_non_normalized += kinship_sparse_matrix[second_vertex][vertex_parent]
        first_second_kinship = np.float16(first_second_kinship_non_normalized / 2)
        # first_second_kinship = kinship.calculate_pair_kinship_sparse(kinship_sparse_matrix,
        #                                                              tuple(self._pred[first_vertex]),
        #                                                              second_vertex)
        kinship_sparse_matrix[first_vertex][second_vertex] = first_second_kinship
        kinship_sparse_matrix[second_vertex][first_vertex] = first_second_kinship

    def get_equivalence_classes(self):
        parents_to_children = defaultdict(list)
        for vertex in self:
            vertex_parents = frozenset(self.get_parents(vertex))
            parents_to_children[vertex_parents].append(vertex)
        vertex_to_representative = dict()
        representative_to_class = dict()
        factor_children_map = dict()
        for orphan_vertex in parents_to_children[frozenset([])]:
            vertex_to_representative[orphan_vertex] = orphan_vertex
            representative_to_class[orphan_vertex] = (orphan_vertex,)
        parents_to_children.pop(frozenset([]))
        for vertex_parents, equivalence_class in parents_to_children.items():
            representative = sorted(list(equivalence_class), key=lambda x: len(self.get_children(x)))[0]
            for vertex in equivalence_class:
                vertex_to_representative[vertex] = representative
            representative_to_class[representative] = equivalence_class
        for representative, equivalence_class in representative_to_class.items():
            factor_children = []
            for vertex in equivalence_class:
                factor_children.extend([vertex_to_representative[x] for x in self.get_children(vertex)])
            factor_children_map[representative] = frozenset(factor_children)
        return vertex_to_representative, representative_to_class, factor_children_map

    def _calculate_pair_kinship_sparse_equivalence(self, kinship_sparse_matrix: dict, first_vertex: int,
                                                   second_vertex: int,
                                                   vertex_to_representative: dict,
                                                   representative_to_sibling_kinship: dict):
        assert second_vertex == vertex_to_representative[second_vertex]
        first_second_kinship_non_normalized = 0
        for vertex_parent in self.get_parents(first_vertex):
            vertex_parent_representative = vertex_to_representative[vertex_parent]
            if (vertex_parent != vertex_parent_representative and
                    second_vertex == vertex_parent_representative):
                first_second_kinship_non_normalized += representative_to_sibling_kinship[vertex_parent_representative]
            else:
                first_second_kinship_non_normalized += kinship_sparse_matrix[second_vertex][
                    vertex_parent_representative]
        first_second_kinship = first_second_kinship_non_normalized / 2
        kinship_sparse_matrix[first_vertex][second_vertex] = first_second_kinship
        kinship_sparse_matrix[second_vertex][first_vertex] = first_second_kinship

    def _calculate_internal_equivalence_non_normalized_kinship(self, kinship_sparse_matrix: dict,
                                                               non_representative: int,
                                                               representative: int,
                                                               vertex_to_representative: dict):
        # This function returns non-normalized kinship!
        first_second_kinship_non_normalized = 0
        for vertex_parent in self.get_parents(non_representative):
            vertex_parent_representative = vertex_to_representative[vertex_parent]
            first_second_kinship_non_normalized += kinship_sparse_matrix[representative][vertex_parent_representative]
        return first_second_kinship_non_normalized / 2

    def _calculate_self_kinship_sparse_equivalence(self, kinship_sparse_matrix: dict, vertex: int,
                                                   vertex_to_representative: dict):
        vertex_parents = self.get_parents(vertex)
        if len(vertex_parents) != 2:
            kinship = 0.5
        else:
            [first_parent, second_parent] = [vertex_to_representative[x] for x in vertex_parents]
            if first_parent < second_parent:
                kinship = (1 + kinship_sparse_matrix[first_parent][second_parent]) / 2
            else:
                kinship = (1 + kinship_sparse_matrix[second_parent][first_parent]) / 2
        kinship_sparse_matrix[vertex][vertex] = kinship

    def _process_sublist(self, vertex, sublist, kinship_sparse_matrix):
        for processed_vertex in sublist:
            self._calculate_pair_kinship_sparse(kinship_sparse_matrix=kinship_sparse_matrix,
                                                first_vertex=vertex,
                                                second_vertex=processed_vertex)
    def calculate_sparse_kinship_queue_c(self):
        children_map = {x: self.get_children(x) for x in self}
        parents_map = {x: self.get_parents(x) for x in self}
        probands = frozenset(self.get_sink_vertices())
        kinship_sparse_matrix = defaultdict(dict)
        # Priority queue which only contains those individuals whose parents have been processed
        founders = frozenset({x for x in self if self.is_orphan(x)})
        queue = []
        parent_to_remaining_children = {x: len(self.get_children(x)) for x in self}
        child_to_remaining_parents = {x: len(self.get_parents(x)) for x in self}
        kinship_sparse_matrix = kinship.calculate_kinship_sparse(
            sink_vertices=probands, founders=founders, parent_to_remaining_children=parent_to_remaining_children,
            child_to_remaining_parents=child_to_remaining_parents, children=children_map,
            parents=parents_map, order=self.order(), counter_limit=200
        )
        return kinship_sparse_matrix

    def calculate_sparse_kinship_queue(self):
        probands = frozenset(self.get_sink_vertices())
        kinship_sparse_matrix = defaultdict(dict)
        # Priority queue which only contains those individuals whose parents have been processed
        founders = frozenset({x for x in self if self.is_orphan(x)})
        queue = []
        parent_to_remaining_children = {x: len(self.get_children(x)) for x in self}
        child_to_remaining_parents = {x: len(self.get_parents(x)) for x in self}
        for founder in founders:
            # heapq.heappush(queue, (-len(self.get_vertex_spouses(founder)), (founder,)))
            heapq.heappush(queue, (0, (founder,)))
        vertices_in_the_cut = set(founders)
        processed_vertices = 0
        # TODO: heappushpop is more efficient than heappush followed by heappop
        counter = 0
        counter_limit = 200
        while queue:
            counter += 1
            if counter == counter_limit:
                print(f"The size of the cut: {len(vertices_in_the_cut)}")
                print(f"Queue size: {len(queue)}")
                print(f"Progress {processed_vertices / self.order()}")
                counter = 0
            children_number, vertices = heapq.heappop(queue)
            for vertex in vertices:
                processed_vertices += 1
                vertices_in_the_cut.add(vertex)
                self._calculate_self_kinship_sparse(kinship_sparse_matrix=kinship_sparse_matrix, vertex=vertex)
                for processed_vertex in kinship_sparse_matrix:
                    if processed_vertex == vertex:
                        continue
                    self._calculate_pair_kinship_sparse(kinship_sparse_matrix=kinship_sparse_matrix,
                                                        first_vertex=vertex,
                                                        second_vertex=processed_vertex)
                for parent in self.get_parents(vertex):
                    parent_to_remaining_children[parent] -= 1
                    if parent_to_remaining_children[parent] == 0:
                        if parent not in probands:
                            parent_to_remaining_children.pop(parent)
                            vertices_in_the_cut.remove(parent)
                            kinship_sparse_matrix.pop(parent)
                            for other_vertex in kinship_sparse_matrix:
                                kinship_sparse_matrix[other_vertex].pop(parent)
                children_to_add = set()
                for child in self.get_children(vertex):
                    child_to_remaining_parents[child] -= 1
                    if child_to_remaining_parents[child] == 0:
                        children_to_add.add(child)
                        child_to_remaining_parents.pop(child)
                if children_to_add:
                    additional_space = len(children_to_add)
                    children_parents = self.get_vertices_parents(children_to_add)
                    for child_parent in children_parents:
                        to_be_removed = True
                        for child in self.get_children(child_parent):
                            if (child not in children_to_add and child not in parent_to_remaining_children
                                    and child not in kinship_sparse_matrix):
                                to_be_removed = False
                                break
                        if to_be_removed:
                            additional_space -= 1
                    heapq.heappush(queue, (additional_space, tuple(children_to_add)))
        return kinship_sparse_matrix

    def calculate_sparse_kinship_queue_equivalence(self):
        def get_vertex_parents_equivalence_number_of_spouses(vertex: int):
            spouses = set()
            for parent in self.get_parents(vertex):
                parent_representative = vertex_to_representative[parent]
                for child in factor_children_map[parent_representative]:
                    spouses.update([vertex_to_representative[x] for x in self.get_parents(child)])
            return len(spouses)

        (vertex_to_representative, representative_to_class,
         factor_children_map) = self.get_equivalence_classes()
        probands = frozenset(self.get_sink_vertices())
        kinship_sparse_matrix = defaultdict(dict)
        # Priority queue which only contains those individuals whose parents have been processed
        orphans = [x for x in self if self.is_orphan(x)]
        queue = []
        parent_to_remaining_children = {x: len(factor_children_map[x]) for x in factor_children_map}
        child_to_remaining_parents = {x: len(self.get_parents(x)) for x in factor_children_map}
        representative_to_sibling_kinship = dict()
        for orphan in orphans:
            heapq.heappush(queue, (get_vertex_parents_equivalence_number_of_spouses(orphan), orphan))
            # queue.append(orphan)
        processed_vertices = 0
        while queue:
            print(f"Queue size: {len(queue)}")
            print(f"Progress {processed_vertices / len(factor_children_map)}")
            priority, vertex = heapq.heappop(queue)
            # vertex = queue.pop(0)
            processed_vertices += 1
            self._calculate_self_kinship_sparse_equivalence(kinship_sparse_matrix=kinship_sparse_matrix, vertex=vertex,
                                                            vertex_to_representative=vertex_to_representative)
            for processed_vertex in kinship_sparse_matrix:
                if processed_vertex == vertex:
                    continue
                self._calculate_pair_kinship_sparse_equivalence(kinship_sparse_matrix=kinship_sparse_matrix,
                                                                first_vertex=vertex,
                                                                second_vertex=processed_vertex,
                                                                vertex_to_representative=vertex_to_representative,
                                                                representative_to_sibling_kinship=
                                                                representative_to_sibling_kinship)
            if len(representative_to_class[vertex]) > 1:
                # Calculate the kinship within the equivalence class
                non_representative = [x for x in representative_to_class[vertex] if x != vertex][0]
                kinship = self._calculate_internal_equivalence_non_normalized_kinship(
                    kinship_sparse_matrix=kinship_sparse_matrix,
                    non_representative=non_representative,
                    representative=vertex,
                    vertex_to_representative=vertex_to_representative)
                representative_to_sibling_kinship[vertex] = kinship
            factor_parents = [vertex_to_representative[x] for x in self.get_parents(vertex)]
            for parent in factor_parents:
                parent_to_remaining_children[parent] -= 1
                if parent_to_remaining_children[parent] == 0:
                    if parent not in probands:
                        kinship_sparse_matrix.pop(parent)
                        for other_vertex in kinship_sparse_matrix:
                            kinship_sparse_matrix[other_vertex].pop(parent)
            for child in factor_children_map[vertex]:
                child_to_remaining_parents[child] -= 1
                if child_to_remaining_parents[child] == 0:
                    # queue.append(child)
                    heapq.heappush(queue, (get_vertex_parents_equivalence_number_of_spouses(child), child))
        return kinship_sparse_matrix, vertex_to_representative

    # def get_min_levels(self):
    #     min_levels = []
    #     current_level = self.get_levels()[0]
    #     while current_level:
    #         next_level = set()
    #         [next_level.update(self.get_parents(vertex)) for vertex in current_level]
    #         current_level = sorted(current_level, key=lambda vertex: self.vertex_to_level_map[vertex], reverse=True)
    #         min_levels.append(current_level)
    #         for vertex in current_level:
    #             if vertex in next_level:
    #                 next_level.remove(vertex)
    #         current_level = list(next_level)
    #     return min_levels
    #
    # def calculate_kinship_for_probands_sparse(self):
    #     kinship_sparse_matrix = defaultdict(dict)
    #     parent_to_remaining_children = {x: len(self.get_children(x)) for x in self}
    #     reversed_levels = reversed(self.get_levels())
    #     for i, level in enumerate(reversed_levels):
    #         print(f"Level {i}: {len(level)}")
    #         # processed_vertices = list(kinship_sparse_matrix.keys())
    #         processed_vertices = []
    #         for vertex in level:
    #             print(f"Processing vertex: {vertex}")
    #             print(f"The size of the kinship matrix is {len(processed_vertices)}")
    #             processed_vertices.append(vertex)
    #             for processed_vertex in processed_vertices:
    #                 if processed_vertex == vertex:
    #                     self._calculate_self_kinship_sparse(kinship_sparse_matrix=kinship_sparse_matrix, vertex=vertex)
    #                 elif processed_vertex not in kinship_sparse_matrix:
    #                     continue
    #                 else:
    #                     self._calculate_pair_kinship_sparse(kinship_sparse_matrix=kinship_sparse_matrix,
    #                                                         first_vertex=vertex,
    #                                                         second_vertex=processed_vertex)
    #             for parent in self.get_parents(vertex):
    #                 parent_to_remaining_children[parent] -= 1
    #                 if not parent_to_remaining_children[parent]:
    #                     print("Removing a parent")
    #                     parent_to_remaining_children.pop(parent)
    #                     kinship_sparse_matrix.pop(parent)
    #                     for other_vertex in kinship_sparse_matrix:
    #                         kinship_sparse_matrix[other_vertex].pop(parent)
    #     return kinship_sparse_matrix
    #
    # def _calculate_pair_kinship_from_parents(self, sparse_kinship_matrix: dict, first_vertex: int, second_vertex: int):
    #     first_vertex_parents = self.get_parents(first_vertex)
    #     second_vertex_parents = self.get_parents(second_vertex)
    #     non_normalized_kinship = 0
    #     for first_vertex_parent in first_vertex_parents:
    #         for second_vertex_parent in second_vertex_parents:
    #             non_normalized_kinship += sparse_kinship_matrix[first_vertex_parent][second_vertex_parent]
    #     kinship = non_normalized_kinship / 4
    #     sparse_kinship_matrix[first_vertex][second_vertex] = kinship
    #     sparse_kinship_matrix[second_vertex][first_vertex] = kinship
    #
    # def calculate_kinship_by_levels(self):
    #     kinship_sparse_matrix = defaultdict(dict)
    #     probands = self.get_sink_vertices()
    #     previous_level = []
    #     current_level = self.get_levels()[-1]
    #     while current_level:
    #         print("Processing next level")
    #         processed_vertices = []
    #         for vertex in current_level:
    #             if vertex in previous_level:
    #                 continue
    #                 # vertex_parents = self.get_parents(vertex)
    #                 # if len(vertex_parents) == 2 and vertex_parents[0] not in kinship_sparse_matrix[vertex_parents[1]]:
    #                 #     [first_parent, second_parent] = vertex_parents
    #                 #     previous_level_vertex = first_parent if first_parent in previous_level else second_parent
    #                 #     other_vertex = second_parent if second_parent != previous_level_vertex else first_parent
    #                 #     self._calculate_pair_kinship_from_parents(sparse_kinship_matrix=kinship_sparse_matrix,
    #                 #                                               first_vertex=other_vertex,
    #                 #                                               second_vertex=previous_level_vertex)
    #             self._calculate_self_kinship_sparse(kinship_sparse_matrix, vertex)
    #             for processed_vertex in kinship_sparse_matrix:
    #                 if processed_vertex == vertex:
    #                     continue
    #                 self._calculate_pair_kinship_sparse(kinship_sparse_matrix=kinship_sparse_matrix,
    #                                                     first_vertex=vertex,
    #                                                     second_vertex=processed_vertex)
    #                 # for vertex_parent in self.get_parents(vertex):
    #                 #     if processed_vertex not in kinship_sparse_matrix[vertex_parent]:
    #                 #         self._calculate_pair_kinship_sparse(sparse_kinship_matrix=kinship_sparse_matrix,
    #                 #                                             first_vertex=processed_vertex,
    #                 #                                             second_vertex=vertex_parent)
    #                 # self._calculate_pair_kinship_sparse(sparse_kinship_matrix=kinship_sparse_matrix,
    #                 #                                     first_vertex=vertex,
    #                 #                                     second_vertex=processed_vertex)
    #             processed_vertices.append(vertex)
    #         for vertex in previous_level:
    #             if vertex not in current_level:
    #                 kinship_sparse_matrix.pop(vertex)
    #                 for other_vertex in kinship_sparse_matrix:
    #                     kinship_sparse_matrix[other_vertex].pop(vertex)
    #         next_level = set()
    #         [next_level.update(self.get_children(vertex)) for vertex in current_level]
    #         if not next_level:
    #             break
    #         [next_level.add(x) for x in current_level if x in probands]
    #         next_level = sorted(next_level, key=lambda v: self.vertex_to_level_map[v], reverse=True)
    #         previous_level = set(current_level)
    #         current_level = next_level
    #     return kinship_sparse_matrix
    #
    # def _calculate_pair_kinship_half_version(self, kinship_sparse_matrix: dict, first_vertex: int, second_vertex: int):
    #     first_second_kinship_non_normalized = 0
    #     for vertex_parent in self.get_parents(first_vertex):
    #         if vertex_parent in kinship_sparse_matrix[second_vertex]:
    #             first_second_kinship_non_normalized += kinship_sparse_matrix[second_vertex][vertex_parent]
    #     first_second_kinship = first_second_kinship_non_normalized / 2
    #     kinship_sparse_matrix[first_vertex][second_vertex] = first_second_kinship
    #     kinship_sparse_matrix[second_vertex][first_vertex] = first_second_kinship
    #
    # def calculate_kinship_for_probands_iteratively(self):
    #     def process_vertex(vertex: int):
    #         self._calculate_self_kinship_sparse(kinship_sparse_matrix, vertex)
    #         for processed_vertex in kinship_sparse_matrix:
    #             if processed_vertex == vertex:
    #                 continue
    #             self._calculate_pair_kinship_sparse(kinship_sparse_matrix, vertex, processed_vertex)
    #
    #     probands = self.get_sink_vertices()
    #     current_level = {x for x in self if self.is_orphan(x)}
    #     kinship_sparse_matrix = defaultdict(dict)
    #     kinship_children_map = {x: set(self.get_children(x)) for x in self}
    #     for vertex in current_level:
    #         process_vertex(vertex)
    #     while current_level:
    #         next_level = set()
    #         for vertex in current_level:
    #             children_to_remove = []
    #             for child in kinship_children_map[vertex]:
    #                 if len([x for x in self.get_parents(child) if x in kinship_sparse_matrix]) != len(self.get_parents(child)):
    #                     next_level.add(vertex)
    #                 else:
    #                     next_level.add(child)
    #                     children_to_remove.append(child)
    #             for child in children_to_remove:
    #                 kinship_children_map[vertex].remove(child)
    #             if not kinship_children_map[vertex] and vertex not in probands:
    #                 kinship_children_map.pop(vertex)
    #         print(f"Current level length {len(current_level)}")
    #         print(f"Next level length {len(next_level)}")
    #         level_difference = current_level - next_level
    #         print(f"Level difference {len(level_difference)}")
    #         processed_vertices = []
    #         for vertex in next_level:
    #             if vertex in current_level:
    #                 for processed_vertex in processed_vertices:
    #                     self._calculate_pair_kinship_sparse(kinship_sparse_matrix, processed_vertex, vertex)
    #             else:
    #                 process_vertex(vertex)
    #             processed_vertices.append(vertex)
    #         print(f"Matrix size before cleanup: {len(kinship_sparse_matrix)}")
    #         for vertex in level_difference:
    #             if vertex in probands:
    #                 continue
    #             kinship_sparse_matrix.pop(vertex)
    #             for other_vertex in kinship_sparse_matrix:
    #                 kinship_sparse_matrix[other_vertex].pop(vertex)
    #         print(f"Matrix size after cleanup: {len(kinship_sparse_matrix)}")
    #         # matrix_size = objsize.get_deep_size(kinship_sparse_matrix)
    #         # print(f"Matrix memory allocated {matrix_size} bytes")
    #         current_level = next_level
    #     return kinship_sparse_matrix
    #
    # def calculate_kinship_for_probands_half(self):
    #     def find_unprocessed_parents(current_vertex_parents):
    #         vertex_unprocessed_parents = []
    #         for vertex_parent in current_vertex_parents:
    #             if vertex_parent not in kinship_sparse_matrix or vertex_parent in half_processed_vertices:
    #                 vertex_unprocessed_parents.append(vertex_parent)
    #         return vertex_unprocessed_parents
    #
    #     def notify_parent_child_processed(parent):
    #         memory_vertices = []
    #         parent_to_remaining_children[parent] -= 1
    #         if not parent_to_remaining_children[parent]:
    #             parent_to_remaining_children.pop(parent)
    #             vertex_to_forced_skipped_vertices.pop(parent, [])
    #             parent_children = self.get_children(parent)
    #             vertices_to_keep = []
    #             for parent_child in parent_children:
    #                 # Optimization: Can only perform one look-up
    #                 if parent_child in half_processed_vertices:
    #                     other_parent = [x for x in self.get_parents(parent_child) if x != parent][0]
    #                     other_parent_ancestors = self.get_ascending_graph_from_vertices((other_parent,))
    #                     vertices_to_keep.extend(other_parent_ancestors)
    #                     for other_parent_ancestor in other_parent_ancestors:
    #                         vertex_to_forced_skipped_vertices[other_parent_ancestor].add(parent)
    #             for other_vertex in list(kinship_sparse_matrix.keys()):
    #                 if parent == other_vertex:
    #                     kinship_sparse_matrix[parent].pop(parent)
    #                     continue
    #                 if other_vertex not in vertices_to_keep:
    #                     kinship_sparse_matrix[parent].pop(other_vertex, 0)
    #                     kinship_sparse_matrix[other_vertex].pop(parent, 0)
    #                     if not kinship_sparse_matrix[other_vertex]:
    #                         kinship_sparse_matrix.pop(other_vertex)
    #             if len(kinship_sparse_matrix[parent]) == 0:
    #                 kinship_sparse_matrix.pop(parent)
    #             memory_vertices.append(parent)
    #         return memory_vertices
    #
    #     def update_vertex_processed(processed_vertex):
    #         memory_vertices = []
    #         for parent in self.get_parents(processed_vertex):
    #             memory_vertices.extend(notify_parent_child_processed(parent))
    #         return memory_vertices
    #
    #     parent_to_remaining_children = {x: len(self.get_children(x)) for x in self if x not in self.get_sink_vertices()}
    #     kinship_sparse_matrix = defaultdict(dict)
    #     vertices_to_skip = []
    #     vertex_to_forced_skipped_vertices = defaultdict(set)
    #     half_processed_vertices = dict()
    #     kinship_queue = [x for x in self if not self.get_parents(x)]
    #     while kinship_queue:
    #         vertex = kinship_queue.pop(0)
    #         kinship_entries = {x for x in kinship_sparse_matrix.keys() if x not in vertices_to_skip}
    #         kinship_entries.update(vertex_to_forced_skipped_vertices.get(vertex, []))
    #         # If the vertex has been processed, continue
    #         if vertex in kinship_sparse_matrix:
    #             if vertex in half_processed_vertices:
    #                 # Update the kinships
    #                 last_processed_parent = half_processed_vertices.pop(vertex)
    #                 self._calculate_self_kinship_sparse(kinship_sparse_matrix=kinship_sparse_matrix, vertex=vertex)
    #                 for processed_vertex in kinship_entries:
    #                     if processed_vertex == vertex:
    #                         continue
    #                     kinship_part = kinship_sparse_matrix[processed_vertex][last_processed_parent] / 2
    #                     kinship_sparse_matrix[vertex][processed_vertex] += kinship_part
    #                     kinship_sparse_matrix[processed_vertex][vertex] += kinship_part
    #                 # kinship_queue.extend(self.get_children(vertex))
    #                 kinship_queue[0:0] = self.get_children(vertex)
    #                 vertices_to_skip.extend(notify_parent_child_processed(last_processed_parent))
    #             continue
    #         vertex_parents = self.get_parents(vertex)
    #         unprocessed_parents = find_unprocessed_parents(vertex_parents)
    #         if unprocessed_parents:
    #             # Process as half-processed vertex
    #             assert len(unprocessed_parents) == 1
    #             unprocessed_parent = unprocessed_parents[0]
    #             processed_parent = [x for x in self.get_parents(vertex) if x != unprocessed_parent][0]
    #             half_processed_vertices[vertex] = unprocessed_parent
    #             for processed_vertex in kinship_entries:
    #                 self._calculate_half_kinship(kinship_matrix=kinship_sparse_matrix,
    #                                              processed_vertex=processed_vertex,
    #                                              half_vertex=vertex,
    #                                              half_processed_vertices=half_processed_vertices)
    #             vertices_to_skip.extend(notify_parent_child_processed(processed_parent))
    #             continue
    #         # Process as full vertex
    #         self._calculate_self_kinship_sparse(kinship_sparse_matrix=kinship_sparse_matrix, vertex=vertex)
    #         for processed_vertex in kinship_entries:
    #             if processed_vertex == vertex:
    #                 continue
    #             self._calculate_pair_kinship_half_version(kinship_sparse_matrix=kinship_sparse_matrix,
    #                                                       first_vertex=vertex,
    #                                                       second_vertex=processed_vertex)
    #         kinship_queue[0:0] = self.get_children(vertex)
    #         # kinship_queue.extend(self.get_children(vertex))
    #         vertices_to_skip.extend(update_vertex_processed(vertex))
    #     return kinship_sparse_matrix


def calculate_kinship_for_probands_cutting(self):
    pass
