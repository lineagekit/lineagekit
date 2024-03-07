"""!
@file graph.py
@brief The file contains the realization of the basic graph class that is responsible for parsing a file with a
genealogical graph (i.e. a pedigree or a coalescent tree) and building the parents and children maps for this graph.
"""
from __future__ import annotations

from tskit import Tree
from collections import defaultdict
from typing import Sequence


class SimpleGraph:
    """!
    This class represents a simple genealogical graph and stores the information about the vertices and the edges
    (children-parent relationships) in a graph
    """

    def __init__(self, children_map: defaultdict = None, parents_map: defaultdict = None):
        if parents_map is None:
            parents_map = defaultdict(list)
        if children_map is None:
            children_map = defaultdict(list)
        self.children_map = children_map
        self.parents_map = parents_map
        self.vertices = set(self.children_map) | set(self.parents_map)

    def get_connected_components(self) -> [[int]]:
        """!
        @brief This method returns the connected components of the graph.
        """
        visited = set()
        connected_components = []
        for vertex in self.parents_map:
            if vertex in visited:
                continue
            connected_component = self.get_connected_component_for_vertex(vertex)
            visited.update(connected_component)
            connected_components.append(connected_component)
        return connected_components

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
                for neighbour in self.parents_map.get(source_vertex, []):
                    dfs(neighbour)
                for neighbour in self.children_map.get(source_vertex, []):
                    dfs(neighbour)

        dfs(vertex)
        return component

    def get_ascending_genealogy_from_vertices(self, vertices: [int]) -> {int}:
        """!
        @brief This method returns all the vertices in the ascending genealogy for the given list of vertices.
        @param vertices The vertices for which the ascending genealogy should be calculated.
        """
        result = set(vertices)
        current_level_vertices = set(result)
        while current_level_vertices:
            next_level_vertices = set()
            for vertex in current_level_vertices:
                if vertex not in self.parents_map:
                    continue
                parents = self.parents_map[vertex]
                result.update(parents)
                next_level_vertices.update(parents)
            current_level_vertices = next_level_vertices
        return result

    def remove_vertices(self, vertices: [int]) -> bool:
        """!
        @brief This function removes the vertices from the graph.
        @param vertices The vertices to be removed.
        @return Returns whether the all the vertices were present in the graph.
        """
        result = True
        for vertex in vertices:
            result = result and self.remove_vertex(vertex)
        return result

    def remove_vertex(self, vertex: int) -> bool:
        """!
        @brief This function removes the vertex from the graph.
        @param vertex The vertex to be removed.
        @return Returns whether the vertex was present in the graph
        """
        complimentary_dictionaries = [(self.parents_map, self.children_map), (self.children_map, self.parents_map)]
        for (first, second) in complimentary_dictionaries:
            if vertex in first:
                keys = first[vertex]
                for key in keys:
                    if key in second:
                        second_values_list = second[key]
                        try:
                            second_values_list.remove(vertex)
                        except ValueError:
                            pass
                        if not second_values_list:
                            second.pop(key)
                first.pop(vertex)
        try:
            self.vertices.remove(vertex)
        except KeyError:
            return False
        return True

    def remove_edge(self, parent: int, child: int) -> bool:
        """!
        @brief Removes the given edge from the graph.
        @param parent The parent vertex.
        @param child The child vertex.
        @return Returns whether the edge was present in the graph.
        """
        try:
            self.parents_map[child].remove(parent)
            self.children_map[parent].remove(child)
        except ValueError:
            return False
        return True

    def add_edge(self, parent: int, child: int, update_vertices: bool = True):
        """!
        @brief This function updated the children map for the given parent-child relationship.
        @param parent The parent id.
        @param child The child id.
        @param update_vertices:
        """
        self.children_map[parent].append(child)
        self.parents_map[child].append(parent)
        if update_vertices:
            self.vertices.update((parent, child))

    def add_vertex(self, vertex: int, exists_ok: bool = True):
        """!
        @brief Adds the given vertices to the graph.
        @param vertex The vertex to be added to the graph.
        @param exists_ok Specifies whether this method should throw an exception if the given vertex is already
        in the graph.
        """
        vertices_before = len(self.vertices)
        self.vertices.add(vertex)
        vertices_after = len(self.vertices)
        if not exists_ok and vertices_before + 1 != vertices_after:
            raise ValueError("The specified vertices is already in the graph")

    def add_vertices(self, vertices: Sequence[int], exists_ok: bool = True):
        """!
        @brief Adds the given vertices to the graph.
        @param vertices The vertices to be added to the graph.
        @param exists_ok Specifies whether this method should throw an exception if one the given vertices is already
        in the graph.
        """
        vertices_before = len(self.vertices)
        self.vertices.update(vertices)
        vertices_after = len(self.vertices)
        if not exists_ok and vertices_before + len(vertices) != vertices_after:
            raise ValueError("One of the specified vertices is already in the graph")

    def get_vertices_number(self) -> int:
        """!
        @brief Returns the number of vertices in the graph
        """
        return len(self.vertices)

    @staticmethod
    def get_graph_from_tree(tree: Tree) -> SimpleGraph:
        """!
        @brief This function builds a simple graph from the given tskit tree. Every node in the tree is treated as
        haploid.
        """
        graph = SimpleGraph()
        graph.parents_map = tree.parent_dict
        for (child, parent) in tree.parent_dict.items():
            graph.add_edge(parent=parent, child=child)
            graph.add_vertices((child, parent))
        return graph

    @staticmethod
    def get_graph_from_file(filename: str, ploidy: int, missing_parent_notation=None, separation_symbol=' ') \
            -> SimpleGraph:
        """!
        @brief Parses the genealogical graph from the file specified by the filename. The ploidy parameter specifies
        whether the organisms in the file should be treated as haploid or diploid.
        Notice that the every line of the file must contain at least ploidy + 1 ids, but it can optionally have
        some metadata which is ignored by this class.
        @param filename The path to the file to be used. The file can optionally start with 1 comment line starting with
        the '#' symbol.
        @param ploidy: The number of ploids that an organism possesses. Must be either 1 or 2.
        @param separation_symbol The symbol used to separate the values in a line. By default, a space is used.
        @param missing_parent_notation The list of text sequences representing that the given individual has no parents.
        If not specified, the default values "-1" and "." are used (meaning that both are accepted at the same time).
        @return The processed pedigree.
        """
        if ploidy == 1:
            if missing_parent_notation is not None:
                raise ValueError("If the specified ploidy is 1, missing_parent_notation must be None, "
                                 "as every vertex mentioned in the file must have a parent")
        elif ploidy != 2:
            raise Exception(f"The ploidy must be either 1 or 2, found {ploidy} specified")
        pedigree: SimpleGraph = SimpleGraph()

        def process_line(file_line: str):
            if ploidy == 1:
                pedigree.add_haploid_line(file_line, separation_symbol)
            else:
                pedigree.add_line_from_pedigree(file_line, missing_parent_notation, separation_symbol)

        file = open(filename, 'r')
        lines = file.readlines()
        if not lines[0].__contains__('#'):
            process_line(lines[0])
        lines.pop(0)
        for line in lines:
            process_line(line)
        file.close()
        assert set(pedigree.parents_map.keys()).issubset(pedigree.vertices)
        assert set(pedigree.children_map.keys()).issubset(pedigree.vertices)
        return pedigree

    @staticmethod
    def get_haploid_graph_from_file(filename: str, separation_symbol=' ') -> SimpleGraph:
        """!
        @brief This method processes the input graph considering that every individual is diploid.
        """
        return SimpleGraph.get_graph_from_file(filename=filename, ploidy=1, separation_symbol=separation_symbol)

    @staticmethod
    def get_diploid_graph_from_file(filename, missing_parent_notation=None, separation_symbol=' ') -> SimpleGraph:
        """!
        @brief Parses the pedigree from the file specified by the filename. Every individual is treated as a diploid
        organism.
        @param filename The path to the file to be used. The file can optionally start with 1 comment line starting with
        the '#' symbol.
        @param separation_symbol The symbol used to separate the values in a line. By default, a space is used.
        @param missing_parent_notation The list of text sequences representing that the given individual has no parents.
        If not specified, the default values "-1" and "." are used (meaning that both are accepted at the same time).
        @return The processed pedigree.
        """
        return SimpleGraph.get_graph_from_file(filename=filename, ploidy=2,
                                               missing_parent_notation=missing_parent_notation,
                                               separation_symbol=separation_symbol)

    def add_haploid_line(self, line: str, separation_symbol=' '):
        """!
        @brief Processes the given line and updated the graph, treating every individual as haploid.
        """
        try:
            (child, parent) = list(map(lambda name: int(name), line.strip('\n').split(separation_symbol)))[:2]
        except ValueError:
            raise Exception("Invalid line")
        self.add_edge(parent=parent, child=child)

    def add_line_from_pedigree(self, line: str, missing_parent_notation=None, separation_symbol=' '):
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
        """
        if missing_parent_notation is None:
            missing_parent_notation = ("-1", '.')
        (child, mother, father) = list(map(lambda name: str(name), line.strip('\n').split(separation_symbol)))[:3]
        child_ploid = 2 * int(child)
        for parent in (mother, father):
            if parent not in missing_parent_notation:
                if child_ploid in self.parents_map:
                    raise ValueError("The same individual is specified multiple times in the input file")
                parent = int(parent)
                self.add_edge(2 * parent, child_ploid, update_vertices=False)
                self.add_edge(2 * parent + 1, child_ploid, update_vertices=False)
                self.add_vertices((2 * parent, 2 * parent + 1, child_ploid))
                child_ploid += 1

    def get_subgraph_from_vertices(self, vertices: [int]) -> SimpleGraph:
        """!
        @brief Returns a new graph containing only the specified vertices
        @param vertices The vertices that the resulting graph should contain
        @return The resulting subgraph
        """

        def narrow_function(map_to_narrow: dict):
            return {key: map_to_narrow[key] for key in vertices if key in map_to_narrow}

        children_map_vertices = narrow_function(self.children_map)
        parents_map_vertices = narrow_function(self.parents_map)
        return SimpleGraph(children_map=defaultdict(list, children_map_vertices),
                           parents_map=defaultdict(list, parents_map_vertices))

    def get_sink_vertices(self):
        """!
        @brief Returns the sink vertices in the graph (that is, the individuals that have parents,
        but don't have children).
        """
        return [x for x in self.parents_map.keys() if x not in self.children_map.keys()]

    def get_orphans(self):
        """!
        @brief Returns the vertices that don't have parents specified.
        """
        return [x for x in self.children_map if x not in self.parents_map]
