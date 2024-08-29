from __future__ import annotations
from typing import Iterable

from basic.AbstractPedigree import AbstractPedigree


class PloidPedigree(AbstractPedigree):
    """
    Pedigree class where every node represents a ploid (half an individual). The maternal ploid is connected with
    both the mother's ploid (if present), and the paternal ploid is connected with the father's ploids (if present).
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @staticmethod
    def get_ploid_pedigree_from_file(filepath: str, probands: Iterable[int] = None,
                                     missing_parent_notation=None, separation_symbol=' ',
                                     skip_first_line: bool = False) -> PloidPedigree:
        """
        Parses the genealogical graph from the file specified by the path, creating two vertices per
        individual in the input file.
        If the individual's id is x, then the resulting graph will have two vertices 2 * x and 2 * x + 1 for their
        ploids.
        For example, if the given line is "1 2 3" (representing that the individual with id 1 has two parents with
        ids 2 and 3), then the resulting graph will have the following information:
        parents_map[2] = 4, 5
        parents_map[3] = 6, 7

        Args:
            filepath (str): The path to the file to be used. The file can optionally start with 1 comment line
            starting with the '#' symbol.
            probands (Iterable[int]): The probands for which the ascending genealogy should be calculated.
            By default, all the vertices from the input file are stored
            separation_symbol (str): The symbol used to separate the values in a line. By default, a space is used.
            missing_parent_notation (Iterable[str]): The list of text sequences representing that the given individual
                                    has no parents. If not specified, the default values "-1" and "." are used
                                    (meaning that both are accepted at the same time).
            skip_first_line (bool): Specifies whether the first line in the file should be skipped. Can be useful if the
                            header does not start with a '#' symbol.

        Returns:
            The processed pedigree.
        """
        pedigree: PloidPedigree = PloidPedigree()

        def process_line(file_line: str):
            pedigree.add_line_from_pedigree(line=file_line,
                                            max_parent_number=2,
                                            missing_parent_notation=missing_parent_notation,
                                            separation_symbol=separation_symbol)

        AbstractPedigree._read_file_and_parse_lines(filepath=filepath, skip_first_line=skip_first_line,
                                                    parse_operation=process_line)
        if probands is not None:
            pedigree.reduce_to_ascending_graph(probands=probands)
        return pedigree

    def add_line_from_pedigree(self, line: str, max_parent_number: int,
                               missing_parent_notation=None, separation_symbol=' '):
        """
        This function processes a single line from a pedigree and updates the graph accordingly.

        Args:
            line: The line to be parsed. The line must consists of at least three integer values separated by
                 the separation symbol.
            missing_parent_notation: The list of text sequences representing that the given individual has no parents.
                                    If not specified, the default values "-1" and "." are used
                                    (meaning that both are accepted at the same time).
            separation_symbol: The symbol used to separate the integers in the line. By default, a space is used.
            max_parent_number: The maximum number of parents a vertex can have. Must be either 1 or 2.
        """
        if missing_parent_notation is None:
            missing_parent_notation = ("-1", '.')
        child, *parents = AbstractPedigree.parse_line(line=line,
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

    def _write_levels_as_diploid(self, file, levels: [[int]]):
        """
        Writes the given levels of the graph to a file. Assumes that the graph vertices represent ploids, and
        saves the corresponding individuals (diploid organisms) to the file.

        Args:
            file: The file to which the content should be written.
            levels: The levels that should be written to the file.
        """
        # TODO: Refactor the code to avoid code duplication
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
                    if ploid_id in self:
                        parents = self.get_parents(ploid_id)
                        if len(parents) == 2:
                            [first_parent, _] = parents
                            first_parent_id = first_parent // 2
                        else:
                            print(f"{ploid_id} is in the graph, but has no parents\n")
                    ploid_id += 1
                    if ploid_id in self:
                        parents = self.get_parents(ploid_id)
                        if len(parents) == 2:
                            [second_parent, _] = parents
                            second_parent_id = second_parent // 2
                        else:
                            print(f"{ploid_id} is in the graph, but has no parents\n")
                file.write(f"{vertex_id} {first_parent_id} {second_parent_id}\n")

    def save_ascending_genealogy_as_diploid(self, filepath: str, vertices: Iterable[int]):
        """
        Saves the ascending genealogy for the given list of vertices
        treating every vertex as a ploid of a diploid organism.

        Args:
            filepath: The path to the file.
            vertices: The vertices for which the ascending genealogy should be saved.
        """
        levels = self.get_ascending_graph_from_vertices_by_levels(vertices)
        file = open(filepath, 'w')
        self._write_levels_as_diploid(file, levels)
        file.close()

    def save_as_diploid(self, filepath: str):
        """
        Saves the graph treating every vertex as a ploid of a diploid organism.

        Args:
            filepath: The path to the file to be written to.
        """
        file = open(filepath, 'w')
        self._write_levels_as_diploid(file, self.get_levels())
        file.close()

    @staticmethod
    def get_individual_ids_from_ploids(ploid_ids: Iterable[int]):
        """
        Transforms the given list of ploid ids into the individual ids.

        Returns:
            The individual ids.
        """
        return {x // 2 for x in ploid_ids}

    @staticmethod
    def get_ploids_from_individual_ids(individual_ids: Iterable[int]):
        """
        Transforms the given list of individual ids into the corresponding ploids' ids.

        Returns:
            The ploids' ids.
        """
        return [2 * x for x in individual_ids] + [2 * x + 1 for x in individual_ids]
