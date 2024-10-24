from __future__ import annotations

import itertools

import networkx as nx
import pandas as pd
from tskit import Tree
import warnings
from typing import Iterable, Callable, TextIO


class GenGraph(nx.DiGraph):
    """
    The general class for a genealogical directed graph.
    """

    def __init__(self, parent_number: int = 2, *args, **kwargs):
        super(GenGraph, self).__init__(*args, **kwargs)
        self.parent_number = parent_number
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

    def _initialize_vertex_to_level_map(self):
        """
        Assigns the vertices to their levels. Refer to the docstring for this class to understand how a level of
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
        """
        Invalidates the levels, causing them to be recalculated on the next request.
        """
        if self.levels_valid:
            self.levels = None
            self.vertex_to_level_map = None
            self.levels_valid = False

    def get_levels(self):
        """
        Checks whether the levels data is outdated (which can be caused by removing vertices/edges without
        recalculating the levels). If so, recalculates the levels. Then, returns the levels of the graph.

        Returns:
            The graph's levels.
        """
        if not self.levels_valid:
            self._initialize_vertex_to_level_map()
        return self.levels

    def get_ascending_graph_from_vertices_by_levels(self, vertices: Iterable[int]):
        """
        Returns the ascending graph for the given list of vertices ordered by levels.

        Args:
            vertices (Iterable[int]): The vertices for which the ascending graph should be calculated.

        Returns:
            The ascending graph for the given list of vertices ordered by levels.
        """
        ascending_graph = self.get_ascending_vertices_from_probands(vertices)
        return [[x for x in level if x in ascending_graph] for level in self.get_levels()]

    def get_top_level_vertices(self) -> [int]:
        """
        Returns:
            The vertices at the top level of the graph.
        """
        return self.get_levels()[-1]

    def get_vertices_for_given_level(self, level):
        """
        Args:
            level (int): The level to be used.

        Returns:
            The vertices belonging to the specified level.
        """
        levels_number = len(self.get_levels())
        if level <= -levels_number or level >= levels_number:
            raise ValueError("Invalid level index")
        if level < 0:
            warnings.warn("The level index is negative.", Warning)
        return self.get_levels()[level]

    def get_vertex_level(self, vertex: int):
        """
        Args:
            vertex (int): The vertex for which the level should be calculated.

        Returns:
            The level of the specified vertex.
        """
        if not self.levels_valid:
            self._initialize_vertex_to_level_map()
        return self.vertex_to_level_map[vertex]

    def is_founder(self, vertex: int):
        """
        Args:
            vertex (int): The vertex id.

        Returns:
             True if the vertex belongs to the top level of the graph, False otherwise.
        """
        vertex_level = self.get_vertex_level(vertex)
        return vertex_level == len(self.levels) - 1

    def add_edge(self, child: int, parent: int, **attr) -> None:
        """
        Adds an edge to the graph. If either child or parent vertex is not present in the graph, adds the vertex to
        the graph.

        Args:
            child (int): The child id.
            parent (int): The parent id.
            attr (dict): The attributes of the edge.
        """
        super().add_edge(parent, child, **attr)
        self._invalidate_levels()

    def remove_edge(self, parent, child):
        """
        Removes the edge and invalidates the levels.

        Args:
            parent (int): The parent id.
            child (int): The child id.
        """
        super().remove_edge(parent, child)
        self._invalidate_levels()

    def remove_edges_from(self, ebunch):
        """
        Removes the edges and invalidates the levels.

        Args:
            ebunch (Iterable[int]): The edges to be removed.
        """
        super().remove_edges_from(ebunch)
        self._invalidate_levels()

    def remove_node(self, vertex):
        """
        Removes the node from the graph and invalidates the levels.

        Args:
            vertex (int): The node to be removed.
        """
        super().remove_node(vertex)
        self._invalidate_levels()

    def remove_nodes_from(self, nodes):
        """
        Removes all the given nodes and invalidates the levels.

        Args:
            nodes (Iterable[int]): The nodes to be removed.
        """
        super().remove_nodes_from(nodes)
        self._invalidate_levels()

    def remove_edges_to_children(self, node: int):
        """
        Removes all the edges going to the children.

        Args:
            node (int): The node to be removed.
        """
        self.remove_edges_from(list(self.out_edges(node)))

    def remove_edges_to_parents(self, node: int):
        """
        Removes all the edges going to this vertex from the parents.

        Args:
            node (int): The node to be removed.
        """
        self.remove_edges_from(list(self.in_edges(node)))

    def add_children(self, parent: int, children: Iterable[int], **attr) -> None:
        """
        Adds the children for the specified vertex.

        Args:
            parent (int): The parent id.
            children (Iterable[int]): The children to be added.
            attr (dict): The attributes of the edges.
        """
        super().add_edges_from([(parent, child) for child in children], **attr)
        self._invalidate_levels()

    def add_parents(self, child: int, parents: Iterable[int], **attr) -> None:
        """
        Adds the given parents for the specified vertex.

        Args:
            child (int): The parent id.
            parents (Iterable[int]): The parents to be added.
            attr (dict): The attributes of the edges.
        """
        super().add_edges_from([(parent, child) for parent in parents], **attr)
        self._invalidate_levels()

    def is_parent(self, parent: int, child: int) -> bool:
        """
        Returns:
            True if the given relationship is present in the graph, False otherwise.
        """
        return super().has_edge(parent, child)

    def get_parents(self, vertex: int):
        """
        Returns:
            List of the vertex's parents.
        """
        return list(self.predecessors(vertex))

    def get_children(self, vertex: int):
        """
        Returns:
            List of the vertex's children.
        """
        return list(self.successors(vertex))

    def has_children(self, vertex):
        """
        Returns:
            True if the given vertex has children, False otherwise.
        """
        return any(self.successors(vertex))

    def has_parents(self, vertex):
        """
        Returns:
            True if the given vertex has parents, False otherwise.
        """
        return any(self.predecessors(vertex))

    def is_orphan(self, vertex: int):
        """
        Returns:
            True if the specified vertex has no parents, False otherwise.
        """
        return self.in_degree(vertex) == 0

    def get_connected_components(self):
        """
        Returns:
            The connected components of the graph.
        """
        return list(nx.weakly_connected_components(self))

    def get_sink_vertices(self):
        """
        Returns:
            The sink vertices in the graph (that is, the individuals don't have children).
        """
        return [node for node, out_degree in self.out_degree if out_degree == 0]

    def get_orphans(self):
        """
        Returns:
            The vertices that don't have parents.
        """
        return [node for node, in_degree in self.in_degree if in_degree == 0]

    def get_vertices_number(self) -> int:
        """
        Returns:
            The number of vertices in the graph.
        """
        return len(self.nodes)

    def get_ascending_vertices_from_probands(self, probands: Iterable[int]) -> {int}:
        """
        This method returns all the vertices in the ascending graph for the given list of vertices.

        Args:
            probands (Iterable[int]): The vertices for which the ascending graph should be calculated.

        Returns:
            The vertices in the ascending graph.
        """
        result = set()
        for vertex in probands:
            if vertex in self.nodes:
                result.add(vertex)
                result.update(nx.ancestors(self, vertex))
        return result

    def get_connected_component_for_vertex(self, vertex: int) -> [int]:
        """
        The function finds all the vertices that are in the same connected component as the passed vertex.

        Args:
            vertex (int): The vertex for which the connected component should be calculated.

        Returns:
            The set of vertices that are in the same connected component as the passed vertex.
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
        """
        Returns whether every vertex has no more than the given max_parents_number parents.

        Returns:
            bool: True if the vertex has no more than the given max_parents_number parents, False otherwise.
        """
        for node in self.nodes():
            in_degree = self.in_degree(node)
            if in_degree > max_parents_number:
                return False
        return True

    def reduce_to_ascending_graph(self, probands: Iterable[int]):
        """
        Reduces the given graph, so that only the information for the ascending graph for the given
        probands remains in the end.

        Args:
            probands (Iterable[int]): The probands.
        """
        ascending_graph = self.get_ascending_vertices_from_probands(probands)
        self.remove_nodes_from(set(self.nodes()).difference(ascending_graph))

    @staticmethod
    def get_graph_from_tree(tree: Tree, probands: Iterable[int] = None) -> GenGraph:
        """
        This function builds a simple graph from the given tskit tree. Every node in the tree is treated as
        haploid.

        Args:
            tree (Tree): The tskit tree to be used.
            probands (Iterable[int]): Optional parameter. If specified, the resulting graph is reduced, so that only
                                      the ascending graph for the given list of probands is returned.

        Returns:
            The obtained GenGraph.
        """
        result_graph = GenGraph(tree.parent_dict.items())
        if probands is not None:
            result_graph.reduce_to_ascending_graph(probands=probands)
        return result_graph

    @staticmethod
    def get_graph_from_file(filepath: str, parent_number: int = 2, probands: Iterable[int] = None,
                            missing_parent_notation=None, separation_symbol=' ', skip_first_line: bool = False) \
            -> GenGraph:
        """
        Parses the genealogical graph from the file specified by the path.
        Notice that the every line of the file must contain at most max_parent_number + 1 ids, but it can
        optionally have some metadata which is ignored by this class.

        Args:
            filepath (str): The path to the file to be used. The file can optionally start with 1 comment line starting
                            with the '#' symbol.
            parent_number (int): The maximum number of parents an individual can posses. The default is 2.
            probands (Iterable[int]): Optional parameter. The probands for which the ascending graph should be
                                      calculated. By default, all the vertices from the input file are stored.
            missing_parent_notation: The list of text sequences representing that the given individual has no parents.
                                     If not specified, the default values "-1" and "." are used (meaning that both are
                                    accepted at the same time).
            separation_symbol (str): The symbol used to separate the values in a line. By default, a space is used.
            skip_first_line (bool): Specifies whether the first line in the file should be skipped. Can be useful if the
                                    header does not start with a '#' symbol.

        Returns:
            The processed pedigree.
        """
        pedigree: GenGraph = GenGraph(parent_number=parent_number)

        def process_line(file_line: str):
            pedigree._add_haploid_line(line=file_line, max_parent_number=parent_number,
                                       missing_parent_notation=missing_parent_notation,
                                       separation_symbol=separation_symbol)

        pedigree._read_file_and_parse_lines(filepath=filepath, skip_first_line=skip_first_line,
                                            parse_operation=process_line)
        if probands is not None:
            pedigree.reduce_to_ascending_graph(probands=probands)
        return pedigree

    @staticmethod
    def parse_line(line: str, max_parent_number: int, missing_parent_notation: [str], separation_symbol=' '):
        return list(map(lambda name: int(name) if name not in missing_parent_notation else name,
                        line.strip('\n').split(separation_symbol)[:max_parent_number + 1]))

    def _on_multiple_vertex_definition(self, vertex: int, new_parents):
        if new_parents == self.get_parents(vertex):
            warnings.warn(f"Individual {vertex} is specified multiple times with the same parents", UserWarning)
            return
        warnings.warn(f"Individual {vertex} is specified multiple times in the graph."
                      f"The previous parents are {self.get_parents(vertex)}, new values: {new_parents}", UserWarning)
        self.remove_edges_to_parents(vertex)

    @staticmethod
    def _read_file_and_parse_lines(filepath: str, skip_first_line: bool,
                                   parse_operation: Callable[[str], None]):
        with open(filepath, 'r') as file:
            first_line = file.readline()
            if not skip_first_line and not first_line.__contains__('#'):
                parse_operation(first_line)
            for line in file:
                parse_operation(line)
        file.close()

    def _add_haploid_line(self, line: str, max_parent_number: int, separation_symbol=' ', missing_parent_notation=None):
        """
        Processes the given line and updated the graph, treating every individual as haploid.
        """
        if missing_parent_notation is None:
            missing_parent_notation = ("-1", '.')
        try:
            child, *parents = GenGraph.parse_line(line=line,
                                                  missing_parent_notation=missing_parent_notation,
                                                  max_parent_number=max_parent_number,
                                                  separation_symbol=separation_symbol)
        except ValueError:
            raise Exception(f"Invalid line {line}")
        if child in self and self.has_parents(child):
            self._on_multiple_vertex_definition(child, parents)
        for parent in parents:
            if parent not in missing_parent_notation:
                self.add_edge(child=child, parent=parent)

    def save_to_file(self, filename: str, separator: str = ' ', missing_parent_notation: str = "-1"):
        """
        Saves the graph to a file.

        Args:
            filename (str): The path of the resulting file.
            separator (str): The string used to separate the columns in the file. The default is ' '.
            missing_parent_notation (str): The string used to indicate that the parent is not known.
                                           The default is "-1".
        """
        file = open(filename, 'w')
        self._save_graph_to_file(file=file, separator=separator, missing_parent_notation=missing_parent_notation)
        file.close()

    def _save_graph_to_file(self, file: TextIO, separator: str = ' ', missing_parent_notation: str = "-1"):
        """
        Saves the graph to a file.
        """
        for vertex in self.nodes:
            parents = self.get_parents(vertex)
            vertices_to_write = [vertex] + parents
            vertices_to_write = [str(vertex) for vertex in vertices_to_write]
            if missing_parent_notation is not None and len(parents) < self.parent_number:
                vertices_to_write += [missing_parent_notation] * (self.parent_number - len(parents))
            file.write(f"{separator.join(vertices_to_write)}\n")

    def has_edge(self, parent: int, child: int):
        """
        Returns: Whether the edge is present in the graph
        """
        return super().has_edge(parent, child)

    def to_data_frame(self):
        """
        Convert a GenGraph object to a DataFrame. Formatted for workflows with libraries such as sgkit.
        """
        data = []
        for vertex in self.nodes:
            # Use predecessors to get parents (assuming the first parent is the sire and the second is the dam)
            parents = list(self.predecessors(vertex))[:2]
            # Ensure two entries for parents, filling missing values with "."
            while len(parents) < 2:
                parents.append(".")
            # Add the vertex and its parents to the data list
            data.append([vertex] + parents)

        df = pd.DataFrame(data, columns=["ID", "SIRE", "DAM"])
        # Replace missing parent indicators as needed to align with sgkit
        df.replace({"": ".", None: "."}, inplace=True)
        return df
