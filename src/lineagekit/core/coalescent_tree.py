from __future__ import annotations

from pathlib import Path

import networkx as nx
from typing import Iterable

from tskit import Tree, TreeSequence

from lineagekit.core.gen_graph import GenGraph


class CoalescentTree(GenGraph):
    """
    Special class designed for working with coalescent trees.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(parent_number=1, *args, **kwargs)

    def copy(self, as_view=False):
        copied_graph = super().copy(as_view=False)
        copied_graph._parent_number = self._parent_number
        return copied_graph

    def get_largest_clade_by_size(self) -> [int]:
        """
        Returns:
             The largest clade in the tree by its total size.
        """
        clades = self.get_connected_components()
        largest_clade = max(clades, key=len)
        return largest_clade

    def get_largest_clade_by_probands(self) -> [int]:
        """
        Returns:
             The largest clade in the tree by the number of probands.
        """
        probands = frozenset(self.get_sink_vertices())

        def intersection_size(clade):
            return len(probands.intersection(clade))

        clades = self.get_connected_components()
        largest_clade = max(clades, key=lambda clade: intersection_size(clade))
        return largest_clade

    def get_root_vertex(self) -> int:
        """
        Returns:
            The root vertex of the tree consisting of a single clade. If the root cannot be determined, the
            method throws an exception.

        Raises:
            ValueError: If the tree has multiple source vertices or the tree is empty.
        """
        founders = self.get_founders()
        if len(founders) > 1:
            raise ValueError("There are multiple source vertices in the tree")
        if not founders:
            assert self.get_vertices_number() == 0
            raise ValueError("The tree is empty")
        return founders[0]

    def get_root_for_clade(self, clade: [int]) -> int:
        """
        Returns:
             The root of the given clade.
        """
        max_level_vertex = max(clade, key=lambda x: self.get_vertex_level(x))
        max_level = self.get_vertex_level(max_level_vertex)
        root_vertices = [x for x in clade if x in self.get_levels()[max_level]]
        if len(root_vertices) != 1:
            raise Exception("Invalid clade value")
        assert root_vertices[0] == max_level_vertex
        return root_vertices[0]

    def remove_unary_nodes(self):
        """
        Removes all the unary nodes in the coalescent tree and recalculates the levels of the coalescent tree.
        """
        for level in self.get_levels()[1:].__reversed__():
            intermediate_nodes = []
            for vertex in level:
                # Since the first level is omitted, all the vertices processed here must have children
                if vertex not in self:
                    continue
                children = self.get_children(vertex)
                if len(children) == 1:
                    child = vertex
                    while len(children) == 1:
                        intermediate_nodes.append(child)
                        [child] = children
                        children = self.get_children(child)
                        if not children:
                            break
                    parents = self.get_parents(vertex)
                    if parents:
                        [parent] = parents
                        self.add_edge(child=child, parent=parent)
                        self.remove_edges_to_parents(child)
                        self.add_edge(child=child, parent=parent)
            self.remove_nodes_from(intermediate_nodes)
        assert not [x for x in self.nodes if len(self.get_children(x)) == 1]

    def get_path_between_vertices(self, source: int, target: int) -> list[int]:
        """
        This function finds the path between the given two vertices. If there is no path between the given
        two vertices, None is returned.
        Important! For performance purposes, this function does not verify that there is no more than one path between
        the vertices (which must be true for a tree). Instead, a random path is returned

        Args:
            source (int): The source vertex.
            target (int): The target vertex.

        Returns:
            The path between the given vertices.

        Raises:
            NodeNotFound if either the source or target vertex does not exist
        """
        return nx.shortest_path(self, source=source, target=target)

    def merge_edge(self, parent: int, child: int):
        """
        This method merges the specified edge in the tree. The purpose of this function is to model a
        common error in an ARG, where what was intended to be an edge is, in reality, just a vertex.
        Args:
            parent: The parent vertex of the edge.
            child: The child vertex of the edge
        """
        child_children = self.get_children(child)
        if not child_children:
            raise ValueError(f"The specified {child}-{parent} edge is a proband edge")
        self.add_edges_from((parent, child_child) for child_child in child_children)
        self.remove_edges_from((child, child_child) for child_child in child_children)
        self.remove_edge(parent=parent, child=child)

    def unmerge_edge(self, child: int) -> int:
        def get_new_vertex_id():
            return max(self.nodes()) + 1

        child_parent = self.get_parents(child)
        if not child_parent:
            raise Exception(f"The specified vertex {child} does not have a parent")
        child_parent = child_parent[0]
        child_parent_children = self.get_children(child_parent)
        if len(child_parent_children) < 3:
            raise Exception(f"The specified vertex {child} is not a part of a polytomy")
        self.remove_edge(parent=child_parent, child=child)
        new_vertex_id = get_new_vertex_id()
        self.add_edge(parent=new_vertex_id, child=child)
        child_parent_parent = self.get_parents(child_parent)
        if child_parent_parent:
            child_parent_parent = child_parent_parent[0]
            self.remove_edge(parent=child_parent_parent, child=child_parent)
            self.add_edge(parent=child_parent_parent, child=new_vertex_id)
        else:
            self.remove_edges_to_parents(child_parent)
        self.add_edge(parent=new_vertex_id, child=child_parent)
        return new_vertex_id

    def get_vertex_parent(self, vertex_id: int) -> int | None:
        """
        This function returns the unique parent vertex of the given vertex.
        Args:
            vertex_id: The child vertex id.

        Returns:
            The parent vertex id.
        """
        vertex_parents = self.get_parents(vertex_id)
        if len(vertex_parents) > 1:
            raise Exception(f"The tree is invalid, the vertex {vertex_id} has multiple parents")
        if not vertex_parents:
            return None
        return vertex_parents[0]

    def subdivide_tree(self, edge_child_vertex: int) -> (CoalescentTree, CoalescentTree):
        """
        This method creates the two trees by removing the specified edge.

        Args:
            edge_child_vertex: The child vertex of the edge. Notice that in a coalescent tree, the edge
            can be identified by the id of its child vertex.

        Returns:
            The two trees that are obtained by removing the edge. The first tree corresponds to the
            upper subtree, the second tree corresponds to the lower subtree.
        """
        bottom_tree_vertices = list(self.get_descendants_for_vertex(edge_child_vertex))
        bottom_tree_vertices.append(edge_child_vertex)
        other_vertices = set(self.nodes).difference(bottom_tree_vertices)
        bottom_tree: CoalescentTree = self.copy()
        upper_tree: CoalescentTree = self.copy()
        bottom_tree.remove_nodes_from(other_vertices)
        upper_tree.remove_nodes_from(bottom_tree_vertices)
        upper_tree.remove_unary_nodes()
        bottom_tree.remove_unary_nodes()
        return upper_tree, bottom_tree

    @staticmethod
    def get_coalescent_tree(tree: Tree, probands: Iterable[int] = None) -> CoalescentTree:
        """
        Utility function to get a coalescent tree from a tskit tree.

        Args:
            tree: The tskit tree.
            probands: The probands for which the ascending genealogy should be calculated. Do not specify the value
                      if all the vertices without children should be considered probands.
        Returns:
            The resulting tree.
        """
        result = CoalescentTree()
        for child, parent in tree.parent_dict.items():
            result.add_edge(parent=parent, child=child)
        if probands:
            result.reduce_to_ascending_graph(probands)
        return result

    @staticmethod
    def get_arg(tree_sequence: TreeSequence):
        # TODO: Consider creating a separate ARG class
        first_tree = tree_sequence.first()
        coalescent_tree = CoalescentTree.get_coalescent_tree(first_tree)
        present_edges = {(key, value) for key, value in first_tree.parent_dict.items()}
        for tree in tree_sequence.trees():
            for (child, parent) in tree.parent_dict.items():
                if not (child, parent) in present_edges:
                    coalescent_tree.add_edge(child=child, parent=parent)
        return coalescent_tree

    @staticmethod
    def get_coalescent_tree_from_file(filepath: str | Path, probands: Iterable[int] = None,
                                      missing_parent_notation=None, separation_symbol=' ',
                                      skip_first_line: bool = False) -> CoalescentTree:
        """
        Utility function to get a coalescent tree from a file.

        The format of the file is as follows:
        child parent
        """
        pedigree: GenGraph = GenGraph.get_graph_from_file(filepath=filepath, probands=probands, parent_number=1,
                                                          missing_parent_notation=missing_parent_notation,
                                                          separation_symbol=separation_symbol,
                                                          skip_first_line=skip_first_line)
        result = CoalescentTree()
        result.update(pedigree)
        return result
