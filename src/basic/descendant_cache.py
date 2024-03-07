"""!
@file descendant_cache.py
@brief This file contains the interface for a descendant cache.
"""

from abc import abstractmethod


class DescendantCache:
    """Represents the interface a graphe cache storing the descendants for every vertex in the graph"""

    def record_descendants(self, left_parent_id: int, right_parent_id: int, child_id: int):
        self.record_child_descendants(left_parent_id, child_id)
        self.record_child_descendants(right_parent_id, child_id)

    @abstractmethod
    def record_proband(self, proband: int):
        """!
        @brief Stores all the important information about a proband vertex.
        @param proband The proband vertex.
        """
        pass

    @abstractmethod
    def get_vertex_descendants(self, vertex_id: int):
        """!
        @brief Returns the descendants of a vertex.
        """
        pass

    @abstractmethod
    def record_child_descendants(self, parent_id, child_id):
        """!
        @brief Records the parent-child relationship.
        """
        pass

    @abstractmethod
    def remove_vertex(self, vertex: int):
        """!
        @brief Removes the vertex from the descendants map.
        """
        pass
