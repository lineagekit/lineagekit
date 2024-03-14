from __future__ import annotations

import itertools

import networkx as nx
from tskit import Tree
import warnings
from typing import Iterable


class GenGraph(nx.DiGraph):
    def __init__(self, *args, **kwargs):
        super(GenGraph, self).__init__(*args, **kwargs)
        self.vertex_to_level_map = dict()
        self.levels = []
        self.levels_valid = False

    def copy(self, as_view=False):
        # TODO: Implement the non-copying version
        copied_graph = super().copy(as_view=False)
        copied_graph.vertex_to_level_map = dict(self.vertex_to_level_map)
        copied_graph.levels = list(self.levels)
        copied_graph.levels_valid = self.levels_valid
        return copied_graph

    def initialize_vertex_to_level_map(self):
        """!
        @brief Assigns the vertices to their levels. Refer to the docstring for this class to understand how a level of
        a vertex is defined.
        The level of a vertex is defined recursively as follows:
        1) All the proband individuals are assigned to level 0.
        2) If the maximal level among a vertex's children is n, then the vertex's level is n + 1.
        In other words, the level of a vertex is the length of the longest path from a proband to this vertex
        in a graph.
        """
        probands = self.get_sink_vertices()
        self.levels = []
        current_level = 1
        self.vertex_to_level_map = {x: 0 for x in probands}
        current_level_vertices = probands
        while current_level_vertices:
            current_level_vertices = set(itertools.chain.from_iterable(
                [self.get_parents(x) for x in current_level_vertices])
            )
            for vertex in current_level_vertices:
                self.vertex_to_level_map[vertex] = current_level
            current_level += 1
            self.levels.append([])
        for vertex, level in self.vertex_to_level_map.items():
            self.levels[level].append(vertex)
        self.levels_valid = True

    def _invalidate_levels(self):
        """!
        @brief Invalidates the levels, causing them to be recalculated on the next request
        """
        if self.levels_valid:
            self.levels = None
            self.vertex_to_level_map = None
            self.levels_valid = False

    def get_levels(self):
        """!
        @brief Checks whether the levels data is outdated (which can be caused by removing vertices/edges without
        recalculating the levels). If so, recalculates the levels. Then, returns the levels of the graph.
        """
        if not self.levels_valid:
            self.initialize_vertex_to_level_map()
        return self.levels

    def get_ascending_genealogy_from_vertices_by_levels(self, vertices: [int]):
        """!
        @brief Returns the ascending genealogy for the given list of vertices ordered by levels.
        @param vertices The vertices for which the ascending genealogy should be calculated.
        """
        ascending_genealogy = self.get_ascending_graph_from_vertices(vertices)
        return [[x for x in level if x in ascending_genealogy] for level in self.get_levels()]

    def get_top_level_vertices(self) -> [int]:
        """!
        @brief Returns the vertices at the top level of the graph.
        """
        return self.get_levels()[-1]

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

    def get_vertex_level(self, vertex: int):
        if not self.levels_valid:
            self.initialize_vertex_to_level_map()
        return self.vertex_to_level_map[vertex]

    def is_founder(self, vertex: int):
        """!
        @brief Returns whether the vertex belongs to the top level of the graph
        """
        vertex_level = self.get_vertex_level(vertex)
        return vertex_level == len(self.levels) - 1

    def add_edge(self, child: int, parent: int, **attr) -> None:
        """!
        @brief Add an edge to the graph
        """
        super().add_edge(parent, child, **attr)
        self._invalidate_levels()

    def remove_edge(self, parent, child):
        """!
        @brief Removes the edge and invalidates the levels
        """
        super().remove_edge(parent, child)
        self._invalidate_levels()

    def remove_edges_from(self, ebunch):
        """
        @brief Removes the edges and invalidates the levels
        """
        super().remove_edges_from(ebunch)
        self._invalidate_levels()

    def remove_node(self, n):
        """
        @brief Removes the node from the graph and invalidates the levels
        """
        super().remove_node(n)
        self._invalidate_levels()

    def remove_nodes_from(self, nodes):
        """!
        @brief Removes all the given nodes and invalidates the levels
        """
        super().remove_nodes_from(nodes)
        self._invalidate_levels()

    def remove_edges_to_children(self, node: int):
        """!
        @brief Removes all the edges going to the children
        """
        self.remove_edges_from(list(self.out_edges(node)))

    def remove_edges_to_parents(self, node: int):
        """!
        @brief Removes all the edges going to this vertex from the parents
        """
        self.remove_edges_from(list(self.in_edges(node)))

    def add_children(self, parent: int, children: Iterable[int], **attr) -> None:
        """!
        @brief Adds the children for the specified vertex
        """
        super().add_edges_from([(parent, child) for child in children], **attr)
        self._invalidate_levels()

    def add_parents(self, child: int, parents: Iterable[int], **attr) -> None:
        """!
        @brief Adds the given parents for the specified vertex
        """
        super().add_edges_from([(parent, child) for parent in parents], **attr)
        self._invalidate_levels()

    def is_parent(self, parent: int, child: int) -> bool:
        """!
        @brief Verifies whether the given relationship is present in the graph
        """
        return super().has_edge(parent, child)

    def get_parents(self, vertex: int):
        """!
        @brief Returns an iterable view over the vertex's parents
        """
        return list(self.predecessors(vertex))

    def get_children(self, vertex: int):
        """!
        @brief Returns an iterable view over the vertex's children
        """
        return list(self.successors(vertex))

    def has_children(self, vertex):
        """!
        @brief Returns whether the given vertex has children
        """
        return any(self.successors(vertex))

    def has_parents(self, vertex):
        """!
        @brief Returns whether the given vertex has parents
        """
        return any(self.predecessors(vertex))

    def is_orphan(self, vertex: int):
        """!
        @brief Returns whether the specified vertex has no parents
        """
        return self.in_degree(vertex) == 0

    def get_connected_components(self):
        """!
        @brief Returns the connected components of the graph.
        """
        return list(nx.weakly_connected_components(self))

    def get_sink_vertices(self):
        """!
        @brief Returns the sink vertices in the graph (that is, the individuals don't have children).
        """
        return [node for node, out_degree in self.out_degree if out_degree == 0]

    def get_orphans(self):
        """!
        @brief Returns the vertices that don't have parents specified.
        """
        return [node for node, in_degree in self.in_degree if in_degree == 0]

    def get_vertices_number(self):
        """!
        @brief Returns the number of vertices in the graph
        """
        return len(self.nodes)

    def get_ascending_graph_from_vertices(self, vertices: Iterable[int]) -> {int}:
        """!
        @brief This method returns all the vertices in the ascending genealogy for the given list of vertices.
        @param vertices The vertices for which the ascending genealogy should be calculated.
        """
        result = set(vertices)
        for vertex in vertices:
            result.update(nx.ancestors(self, vertex))
        return result

    def get_connected_component_for_vertex(self, vertex: int) -> [int]:
        """!
        @brief The function finds all the vertices that are in the same connected component as the passed vertex.
        @param vertex The vertex for which the connected component should be calculated.
        @return The set of vertices that are in the same connected component as the passed vertex.
        """
        visited = set()
        component = []

        def dfs(source_vertex):
            if source_vertex not in visited:
                visited.add(source_vertex)
                component.append(source_vertex)
                for neighbour in self.get_parents(source_vertex):
                    dfs(neighbour)
                for neighbour in self.get_children(source_vertex):
                    dfs(neighbour)

        dfs(vertex)
        return component

    def verify_max_parents_number(self, max_parents_number: int) -> bool:
        """!
        @brief Returns whether every vertex has no more than the given max_parents_number parents.
        """
        for node in self.nodes():
            in_degree = self.in_degree(node)
            if in_degree > max_parents_number:
                return False
        return True

    def reduce_to_ascending_genealogy(self, probands: Iterable[int]):
        """!
        @brief Reduces the given graph, so that only the information for the ascending genealogy for the given
        probands.
        """
        ascending_genealogy = self.get_ascending_graph_from_vertices(probands)
        self.remove_nodes_from(set(self.nodes()).difference(ascending_genealogy))

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
                if self.has_parents(vertex):
                    ploid_id = 2 * vertex_id
                    parents = self.get_parents(ploid_id)
                    if parents:
                        [first_parent, _] = parents
                        first_parent_id = first_parent // 2
                    ploid_id += 1
                    parents = self.get_parents(ploid_id)
                    if parents:
                        [second_parent, _] = parents
                        second_parent_id = second_parent // 2
                file.write(f"{vertex_id} {first_parent_id} {second_parent_id}\n")

    def save_ascending_genealogy_as_diploid(self, filepath: str, vertices: Iterable[int]):
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

    @staticmethod
    def get_graph_from_tree(tree: Tree, probands: Iterable[int] = None) -> GenGraph:
        """!
        @brief This function builds a simple graph from the given tskit tree. Every node in the tree is treated as
        haploid.
        """
        result_graph = GenGraph(tree.parent_dict.items())
        if probands is not None:
            result_graph.reduce_to_ascending_genealogy(probands=probands)
        return result_graph

    @staticmethod
    def get_graph_from_file(filepath: str, ploidy: int, max_parent_number: int = 2, probands: Iterable[int] = None,
                            missing_parent_notation=None, separation_symbol=' ', skip_first_line: bool = False) \
            -> GenGraph:
        """!
        @brief Parses the genealogical graph from the file specified by the path. The ploidy parameter specifies
        whether the organisms in the file should be treated as haploid or diploid.
        Notice that the every line of the file must contain at least ploidy + 1 ids, but it can optionally have
        some metadata which is ignored by this class.
        @param filepath The path to the file to be used. The file can optionally start with 1 comment line starting with
        the '#' symbol.
        @param max_parent_number The maximum number of parents an individual can posses.
        The value must be either 1 or 2.
        @param probands The probands for which the ascending genealogy should be calculated. By default, all the vertices
        from the input file are stored
        @param ploidy: The number of ploids that an organism possesses. Must be either 1 or 2.
        @param separation_symbol The symbol used to separate the values in a line. By default, a space is used.
        @param missing_parent_notation The list of text sequences representing that the given individual has no parents.
        @param skip_first_line Specifies whether the first line in the file should be skipped. Can be useful if the
        header does not start with a '#' symbol.
        If not specified, the default values "-1" and "." are used (meaning that both are accepted at the same time).
        @return The processed pedigree.

        """
        if ploidy != 1 and ploidy != 2:
            raise Exception(f"The ploidy must be either 1 or 2, found {ploidy} specified")
        pedigree: GenGraph = GenGraph()

        def process_line(file_line: str):
            if ploidy == 1:
                pedigree.add_haploid_line(line=file_line, max_parent_number=max_parent_number,
                                          missing_parent_notation=missing_parent_notation,
                                          separation_symbol=separation_symbol)
            else:
                pedigree.add_line_from_pedigree(line=file_line,
                                                max_parent_number=max_parent_number,
                                                missing_parent_notation=missing_parent_notation,
                                                separation_symbol=separation_symbol)

        with open(filepath, 'r') as file:
            first_line = file.readline()
            if not skip_first_line and not first_line.__contains__('#'):
                process_line(first_line)
            for line in file:
                process_line(line)
        file.close()
        if probands is not None:
            pedigree.reduce_to_ascending_genealogy(probands=probands)
        return pedigree

    @staticmethod
    def get_haploid_graph_from_file(filepath: str, max_parent_number: int = 2, probands: Iterable[int] = None,
                                    missing_parent_notation=None, separation_symbol=' ',
                                    skip_first_line: bool = False) -> GenGraph:
        """!
        @brief This method processes the input graph considering that every individual is diploid.
        """
        return GenGraph.get_graph_from_file(filepath=filepath, ploidy=1,
                                            probands=probands,
                                            max_parent_number=max_parent_number,
                                            missing_parent_notation=missing_parent_notation,
                                            separation_symbol=separation_symbol,
                                            skip_first_line=skip_first_line)

    @staticmethod
    def get_diploid_graph_from_file(filepath: str, max_parent_number: int = 2, probands: Iterable[int] = None,
                                    missing_parent_notation=None, separation_symbol=' ',
                                    skip_first_line: bool = False) -> GenGraph:
        """!
        @brief Parses the pedigree from the file specified by the filename. Every individual is treated as a diploid
        organism.
        @param max_parent_number The maximum number of parents an individual can posses.
        The value must be either 1 or 2.
        @param probands The probands for which the ascending genealogy should be calculated. By default, all the vertices
        from the input file are stored
        @param filepath The path to the file to be used. The file can optionally start with 1 comment line starting with
        the '#' symbol.
        @param separation_symbol The symbol used to separate the values in a line. By default, a space is used.
        @param missing_parent_notation The list of text sequences representing that the given individual has no parents.
        If not specified, the default values "-1" and "." are used (meaning that both are accepted at the same time).
        @param skip_first_line Specifies whether the first line in the file should be skipped. Can be useful if the
        header does not start with a '#' symbol.
        @return The processed pedigree.
        """
        return GenGraph.get_graph_from_file(filepath=filepath, ploidy=2, probands=probands,
                                            max_parent_number=max_parent_number,
                                            missing_parent_notation=missing_parent_notation,
                                            separation_symbol=separation_symbol,
                                            skip_first_line=skip_first_line)

    @staticmethod
    def parse_line(line: str, max_parent_number: int, missing_parent_notation: [str], separation_symbol=' '):
        return list(map(lambda name: int(name) if name not in missing_parent_notation else name,
                        line.strip('\n').split(separation_symbol)[:max_parent_number + 1]))

    def _on_multiple_vertex_definition(self, vertex: int, new_parents):
        self.remove_edges_to_parents(vertex)
        warnings.warn(f"Individual {vertex} is specified multiple times in the graph."
                      f"The previous parents are {self.get_parents(vertex)}, new values: {new_parents}", UserWarning)

    def add_haploid_line(self, line: str, max_parent_number: int, separation_symbol=' ', missing_parent_notation=None):
        """!
        @brief Processes the given line and updated the graph, treating every individual as haploid.
        """
        if missing_parent_notation is None:
            missing_parent_notation = ("-1", '.')
        try:
            child, *parents = GenGraph.parse_line(line=line,
                                                  missing_parent_notation=missing_parent_notation,
                                                  max_parent_number=max_parent_number,
                                                  separation_symbol=separation_symbol)
        except ValueError:
            raise Exception("Invalid line")
        if child in self and self.has_parents(child):
            self._on_multiple_vertex_definition(child, parents)
        for parent in parents:
            if parent not in missing_parent_notation:
                self.add_edge(child=child, parent=parent)

    def add_line_from_pedigree(self, line: str, max_parent_number: int,
                               missing_parent_notation=None, separation_symbol=' '):
        """!
        @brief This function processes a single line from a pedigree and updates the graph accordingly.
        It treats every id as a diploid individual, so if the individual's id is x, then the resulting graph
        will have two vertices 2 * x and 2 * x + 1 for their ploids. If you want to find the individual id by its ploid,
        you can use integer division and divide the id by 2.
        For example, if the given line is "1 2 3" (representing that the individual with id 1 has two parents with
        ids 2 and 3), then the resulting graph will have the following information:
        parents_map[2] = 4, 5
        parents_map[3] = 6, 7
        @param line The line to be parsed. The line must consists of at least three integer values separated by
        the separation symbol.
        @param missing_parent_notation The list of text sequences representing that the given individual has no parents.
        If not specified, the default values "-1" and "." are used (meaning that both are accepted at the same time).
        @param separation_symbol: The symbol used to separate the integers in the line. By default, a space is used.
        @param max_parent_number The maximum number of parents a vertex can have. Must be either 1 or 2
        """
        if missing_parent_notation is None:
            missing_parent_notation = ("-1", '.')
        child, *parents = GenGraph.parse_line(line=line,
                                              missing_parent_notation=missing_parent_notation,
                                              max_parent_number=max_parent_number,
                                              separation_symbol=separation_symbol)
        child_ploid = 2 * int(child)
        if child_ploid in self and self.has_parents(child_ploid):
            self._on_multiple_vertex_definition(vertex=child_ploid, new_parents=parents)
        for parent in parents:
            if parent not in missing_parent_notation:
                parent = int(parent)
                self.add_edge(child=child_ploid, parent=2 * parent)
                self.add_edge(child=child_ploid, parent=2 * parent + 1)
                child_ploid += 1
