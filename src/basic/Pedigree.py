from __future__ import annotations
from typing import Iterable
from src.basic.GenGraph import GenGraph


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
