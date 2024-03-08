"""!
@file genealogical_graph.py
@brief This file contains the realization of the genealogical graph that assigns the vertices to their genealogical
levels and provides a level-by-level framework for preprocessing the graph and the utility class for processing
coalescent trees.
"""

from __future__ import annotations
import itertools

import warnings
import networkx as nx
from matplotlib import pyplot as plt

from src.basic.descendant_memory_cache import DescendantMemoryCache
from src.basic.simple_graph import *


class GenealogicalGraph(SimpleGraph):
    """!
    This class represents a genealogical graph that inherits from the /ref graph::Graph class and additionally
    assigns all the vertices within the graph to their levels. The level of a vertex is defined recursively as follows:
    1) All the proband individuals are assigned to level 0.
    2) If the maximal level among a vertex's children is n, then the vertex's level is n + 1.
    In other words, the level of a vertex is the length of the longest path from a proband to this vertex in a graph.
    """

    def __init__(self, pedigree: SimpleGraph = None, probands: [int] = None, initialize_levels: bool = True):
        if pedigree is None:
            super().__init__()
        else:
            self.parents_map = pedigree.parents_map
            self.children_map = pedigree.children_map
            self.vertices = pedigree.vertices
        if probands is None:
            probands = {x for x in self.parents_map if x not in self.children_map}
        self.probands = probands
        self.levels = []
        self.vertex_to_level_map = dict()
        self.levels_valid = True
        if initialize_levels:
            self.initialize_vertex_to_level_map()
        self.descendant_writer = DescendantMemoryCache()
        self.descendants_valid = False

    def initialize_vertex_to_level_map(self):
        """!
        @brief Assigns the vertices to their levels. Refer to the docstring for this class to understand how a level of
        a vertex is defined.
        """
        self.levels = []
        current_level = 1
        self.vertex_to_level_map = {x: 0 for x in self.probands}
        current_level_vertices = self.probands
        while current_level_vertices:
            current_level_vertices = set(itertools.chain.from_iterable(
                [self.parents_map[x] for x in current_level_vertices if x in self.parents_map])
            )
            for vertex in current_level_vertices:
                self.vertex_to_level_map[vertex] = current_level
            current_level += 1
            self.levels.append([])
        for vertex, level in self.vertex_to_level_map.items():
            self.levels[level].append(vertex)
        self.levels_valid = True

    def get_levels(self):
        """!
        @brief Checks whether the levels data is outdated (which can be caused by removing vertices/edges without
        recalculating the levels). If so, recalculates the levels. Then, returns the levels of the graph.
        """
        if not self.levels_valid:
            self.initialize_vertex_to_level_map()
        return self.levels

    def get_vertices_for_given_level(self, level):
        """!
        @brief Returns the vertices belonging to the specified level.
        @param level The level to be used.
        """
        levels_number = len(self.get_levels())
        if level <= -levels_number or level >= levels_number:
            raise ValueError("Invalid level index")
        if level < 0:
            warnings.warn("The level index is negative.", Warning)
        return self.get_levels()[level]

    def get_proband_descendants(self, vertex_id: int):
        """!
        @param vertex_id The vertex for which the proband descendants should be found.
        @return The descendants of the given vertex which are the probands.
        """
        return [x for x in self.get_vertex_descendants(vertex_id) if x in self.get_probands()]

    @staticmethod
    def get_graph_from_tree(tree: Tree, probands=None) -> GenealogicalGraph:
        """!
        @brief Utility function for creating a SimpleGraph from a tskit tree and creating a GenealogicalGraph.
        @param tree The tree to be processed
        @param probands The probands set to be used while calculating the levels.
        """
        pedigree = SimpleGraph.get_graph_from_tree(tree)
        # TODO: Verify if this is still necessary
        reformatted_parent_dict = dict()
        for key, value in pedigree.parents_map.items():
            reformatted_parent_dict[key] = [value]
        pedigree.parents_map = reformatted_parent_dict
        graph = GenealogicalGraph(pedigree, probands)
        return graph

    def get_vertex_excluded_descendants(self, vertex_id: int, excluded_vertex: int) -> [int]:
        """!
        @brief Returns the descendants of the vertex that can climb to the given vertex by-passing the excluded vertex.
        @param vertex_id The vertex for which the excluded descendants should be calculated
        @param excluded_vertex The vertex that is banned for all the resulting descendants to pass.
        """
        result = set()
        for child in self.children_map[vertex_id]:
            if child != excluded_vertex:
                result.update(self.get_vertex_descendants(child))
        return result

    def get_vertex_descendants(self, vertex_id: int) -> [int]:
        """!
        @brief Returns the descendants of the given vertex. If the descendant map has not been calculated yet,
        calculates it before returning the result.
        """
        if not self.descendants_valid:
            self.process_descendant_map()
            self.descendants_valid = True
        return self.descendant_writer.get_vertex_descendants(vertex_id)

    def is_founder(self, vertex: int):
        """!
        @brief Returns whether the vertex belongs to the top level of the graph
        """
        return self.get_vertex_level(vertex) == len(self.levels) - 1

    def get_vertex_level(self, vertex: int):
        if not self.levels_valid:
            self.initialize_vertex_to_level_map()
        return self.vertex_to_level_map[vertex]

    def get_top_level_vertices(self) -> [int]:
        """!
        @brief Returns the vertices at the top level of the graph.
        """
        return self.get_levels()[-1]

    def get_probands(self):
        """!
        @brief Returns the graph's probands.
        """
        return self.probands

    def remove_vertex(self, vertex: int) -> bool:
        """!
        @brief Removes the vertex from the graph.
        @param vertex The vertex to be removed.
        @returns Whether the vertex was present in the graph.
        """
        removed = super().remove_vertex(vertex)
        # Vertex usually must belong to this map, this can only happen if we try to remove the same vertex twice
        if removed:
            vertex_level = self.vertex_to_level_map[vertex]
            self.vertex_to_level_map.pop(vertex)
            level: [int] = self.levels[vertex_level]
            level.remove(vertex)
            if vertex in self.probands:
                self.probands.remove(vertex)
                self.levels_valid = False
                self.descendants_valid = False
        return removed

    def remove_vertices(self, vertices: [int]) -> bool:
        """!
        @brief Removes the given vertices from the graph.
        @param vertices The vertices to be removed.
        @returns Whether the all the vertices were present in the graph.
        """
        result = True
        for vertex in vertices:
            result = result and self.remove_vertex(vertex)
            self.initialize_vertex_to_level_map()
        return result

    def remove_edge(self, parent: int, child: int) -> bool:
        """!
        @brief Removes the given edge from the graph.
        @param parent The parent vertex.
        @param child The child vertex.
        @return Returns whether the edge was present in the graph.
        """
        result = super().remove_edge(parent=parent, child=child)
        self.descendants_valid = False
        self.levels_valid = False
        return result

    def get_ascending_genealogy_from_vertices_by_levels(self, vertices: [int]):
        """!
        @brief Returns the ascending genealogy for the given list of vertices ordered by levels.
        @param vertices The vertices for which the ascending genealogy should be calculated.
        """
        ascending_genealogy = self.get_ascending_graph_from_vertices(vertices)
        return [[x for x in level if x in ascending_genealogy] for level in self.get_levels()]

    def write_levels_as_diploid(self, file, levels: [[int]]):
        """!
        @brief Writes the given levels of the graph to a file. Assumes that the graph vertices represent ploids, and
        saves the corresponding individuals (diploid organisms) to the file.
        @param file The file to which the content should be written.
        @param levels The levels that should be written to the file.
        """
        processed_ids = set()
        for level in levels:
            for vertex in level:
                vertex_id = vertex // 2
                if vertex_id in processed_ids:
                    continue
                else:
                    processed_ids.add(vertex_id)
                [first_parent_id, second_parent_id] = [-1, -1]
                if vertex in self.parents_map:
                    ploid_id = 2 * vertex_id
                    if ploid_id in self.parents_map:
                        [first_parent, _] = self.parents_map[ploid_id]
                        first_parent_id = first_parent // 2
                    ploid_id += 1
                    if ploid_id in self.parents_map:
                        [second_parent, _] = self.parents_map[ploid_id]
                        second_parent_id = second_parent // 2
                file.write(f"{vertex_id} {first_parent_id} {second_parent_id}\n")

    def draw_graph(self, filename: str):
        """!
        @brief Draws the graph using the networkx drawing library and multipartite layout algorithm. Notice
        that the resulting image will be impossible to read for relatively large graphs.
        @param filename The resulting filename of the image.
        """
        g = nx.DiGraph()
        for level_index, level in enumerate(self.get_levels()):
            for vertex in level:
                g.add_node(vertex, time=level_index)
                vertex_parents = self.parents_map.get(vertex, [])
                for parent in vertex_parents:
                    g.add_edge(vertex, parent)
        fig, ax = plt.subplots(figsize=(655, 655))
        node_size = 1000
        pos = nx.multipartite_layout(g, subset_key="time", align="horizontal")
        nx.draw_networkx(g, pos, ax=ax, node_size=node_size, with_labels=True)
        plt.savefig(filename)

    def save_ascending_genealogy_as_diploid(self, filepath: str, vertices: [int]):
        """!
        @brief Saves the ascending genealogy for the given list of vertices
        treating every vertex as a ploid of a diploid organism.
        @param filepath The path to the file.
        @param vertices The vertices for which the ascending genealogy should be saved.
        """
        levels = self.get_ascending_genealogy_from_vertices_by_levels(vertices)
        file = open(filepath, 'w')
        self.write_levels_as_diploid(file, levels)
        file.close()

    def save_as_diploid(self, filepath: str):
        """!
        @brief Saves the graph treating every vertex as a ploid of a diploid organism.
        @param filepath The path to the file to be written to.
        """
        file = open(filepath, 'w')
        self.write_levels_as_diploid(file, self.get_levels())
        file.close()

    def process_descendant_map(self):
        """!
        @brief This method preprocesses the graph in a level-by-level order calling the descendants cache functions
        for storing the results. Additionally, calls the methods self.process_proband_vertex and
        self.process_level_vertex which can implement additional custom behaviour.
        """
        for x in self.probands:
            self.descendant_writer.record_proband(x)
            self.process_proband_vertex(x)
        counter = 1
        for current_level_vertices in self.get_levels()[1:]:
            for parent in current_level_vertices:
                if parent in self.children_map:
                    children = self.children_map[parent]
                    for child in children:
                        self.process_level_vertex(parent, child, counter)
                        self.descendant_writer.record_child_descendants(parent_id=parent, child_id=child)
            counter += 1
        return self

    def process_level_vertex(self, parent: int, child: int, level: int) -> int:
        """!
        @brief Processes the given parent-child relationship at the given level and stores the results. This method is
        an additional callback that is called within the process_descendant_map method which can be
        overridden by a child
        class for implementing additional level-by-level processing behaviour.
        @param parent The parent vertex.
        @param child The child vertex.
        @param level The level of parent.
        @return The parent vertex
        """
        # Default behaviour
        return parent

    def process_proband_vertex(self, proband: int):
        """! @brief This method is an additional callback that is called within the process_descendant_map method which
             can be overridden by a child class for implementing additional level-by-level processing behaviour.
             @param proband The proband vertex
        """
        # Default behaviour
        self.descendant_writer.record_proband(proband)

    @staticmethod
    def get_diploid_graph_from_file(filename: str, max_parent_number: int = 2,
                                    missing_parent_notation=None, separation_symbol=' ',
                                    skip_first_line: bool = False) -> GenealogicalGraph:
        """!
        @brief Utility function that can be used for getting a diploid genealogical graph from a file.
        @param filename The filename from which the graph will be read.
        @param max_parent_number The maximum number of parents an individual can posses.
        The value must be either 1 or 2.
        @param missing_parent_notation The list of text sequences representing that the given individual has no parents.
        Refer to the documentation for SimpleGraph.get_graph_from_file
        @param separation_symbol The separation sequence used in the file. Refer to the documentation for
        SimpleGraph.get_graph_from_file
        @param skip_first_line Specifies whether the first line in the file should be skipped. Can be useful if the
        header does not start with a '#' symbol.
        @returns The parsed graph.
        """
        pedigree: SimpleGraph = SimpleGraph.get_diploid_graph_from_file(filename=filename,
                                                                        max_parent_number=max_parent_number,
                                                                        missing_parent_notation=missing_parent_notation,
                                                                        separation_symbol=separation_symbol,
                                                                        skip_first_line=skip_first_line)
        return GenealogicalGraph(pedigree=pedigree)

    @staticmethod
    def get_haploid_graph_from_file(filename: str, max_parent_number: int = 2,
                                    missing_parent_notation=None, separation_symbol=' ',
                                    skip_first_line: bool = False) -> GenealogicalGraph:
        """!
        @brief Utility function that can be used for getting a haploid genealogical graph from a file.
        @param filename The filename from which the graph will be read.
        @param max_parent_number The maximum number of parents an individual can posses.
        The value must be either 1 or 2.
        @param missing_parent_notation The list of text sequences representing that the given individual has no parents.
        Refer to the documentation for SimpleGraph.get_graph_from_file
        @param separation_symbol The separation sequence used in the file. Refer to the documentation for
        SimpleGraph.get_graph_from_file
        @param skip_first_line Specifies whether the first line in the file should be skipped. Can be useful if the
        header does not start with a '#' symbol.
        @returns The parsed graph.
        """
        pedigree: SimpleGraph = SimpleGraph.get_haploid_graph_from_file(filename=filename,
                                                                        max_parent_number=max_parent_number,
                                                                        missing_parent_notation=missing_parent_notation,
                                                                        separation_symbol=separation_symbol,
                                                                        skip_first_line=skip_first_line)
        return GenealogicalGraph(pedigree=pedigree)


class CoalescentTree(GenealogicalGraph):
    """!
    This is a helper class that is responsible for working with coalescent trees. Apart from the functionality
    of the GenealogicalGraph, it calculates the connected components (clades) of the graph.
    """

    def __init__(self, pedigree: SimpleGraph):
        super().__init__(pedigree=pedigree)
        self.descendant_writer = DescendantMemoryCache()
        self.process_descendant_map()

    def get_largest_clade_by_size(self) -> [int]:
        """!
        @brief Returns the largest clade in the tree by its total size.
        """
        clades = self.get_connected_components()
        largest_clade = max(clades, key=len)
        return largest_clade

    def get_largest_clade_by_probands(self) -> [int]:
        """!
        @brief Returns the largest clade in the tree by the number of probands.
        """

        def intersection_size(clade):
            return len(self.probands.intersection(clade))

        clades = self.get_connected_components()
        largest_clade = max(clades, key=lambda clade: intersection_size(clade))
        return largest_clade

    def get_root_for_clade(self, clade: [int]) -> int:
        """!
        @brief Returns the root of the given clade.
        """
        max_level_vertex = max(clade, key=lambda x: self.get_vertex_level(x))
        max_level = self.get_vertex_level(max_level_vertex)
        root_vertices = [x for x in clade if x in self.get_levels()[max_level]]
        if len(root_vertices) != 1:
            raise Exception("Invalid clade value")
        assert root_vertices[0] == max_level_vertex
        return root_vertices[0]

    def get_subgraph_from_vertices(self, vertices: [int]) -> CoalescentTree:
        """!
        @brief Returns a new CoalescentTree object which is a subtree of the original tree.
        @param vertices
        @return The subtree
        """
        return CoalescentTree(pedigree=super().get_subgraph_from_vertices(vertices=vertices))

    def remove_unary_nodes(self):
        """!
        @brief Removes all the unary nodes in the coalescent tree and recalculates the levels of the coalescent tree.
        """
        for level in self.get_levels()[1:].__reversed__():
            intermediate_nodes = []
            for vertex in level:
                # Since the first level is omitted, all the vertices processed here must have children
                children = self.children_map[vertex]
                if len(children) == 1:
                    child = vertex
                    while len(children) == 1:
                        intermediate_nodes.append(child)
                        [child] = children
                        if child not in self.children_map:
                            break
                        children = self.children_map[child]
                    if vertex in self.parents_map:
                        [parent] = self.parents_map[vertex]
                        self.children_map[parent].append(child)
                        self.parents_map[child] = [parent]
            self.remove_vertices(intermediate_nodes)
        self.initialize_vertex_to_level_map()
        assert not [x for x in self.children_map if len(self.children_map[x]) == 1]

    def write_levels(self, file, levels):
        """!
        @brief Writes the given levels to a file.
        @param file The file to which the levels will be written to.
        @param levels The levels to be saved.
        """
        for level in levels:
            for vertex in level:
                if vertex in self.parents_map:
                    [parent] = self.parents_map[vertex]
                    file.write(f"{vertex} {parent}\n")

    @staticmethod
    def get_coalescent_tree_from_file(filename) -> CoalescentTree:
        """!
        @brief Utility function to get a coalescent tree from a file.
        """
        pedigree: SimpleGraph = SimpleGraph.get_haploid_graph_from_file(filename)
        result = CoalescentTree(pedigree=pedigree)
        return result

    @staticmethod
    def get_coalescent_tree(tree: Tree, probands=None) -> CoalescentTree:
        """
        @brief Utility function to get a coalescent tree from a tskit tree.
        @param tree The tskit tree.
        @param probands The probands for which the ascending genealogy should be calculated. Do not specify the value
        if all the vertices without children should be considered probands.
        @return The resulting tree.
        """
        genealogical_graph = GenealogicalGraph.get_graph_from_tree(tree, probands)
        return CoalescentTree(genealogical_graph)
