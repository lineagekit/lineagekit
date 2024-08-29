from typing import Tuple, Dict, Set
import numpy as np

class TimeSparseMatrix:
    def get_kinship(self, key1: int, key2: int) -> float: ...
    def to_numpy_and_free(self) -> Tuple[Dict[int, int], np.ndarray]: ...

class MemorySparseMatrix:
    def get_kinship(self, key1: int, key2: int) -> float: ...
    def to_numpy_and_free(self) -> Tuple[Dict[int, int], np.ndarray]: ...

def calculate_kinship_sparse_speed(
    children: Dict[int, int],
    parents: Dict[int, int],
    sink_vertices: Set[int]
) -> TimeSparseMatrix: ...

def calculate_kinship_sparse_memory(
    children: Dict[int, int],
    parents: Dict[int, int],
    sink_vertices: Set[int]
) -> MemorySparseMatrix: ...
