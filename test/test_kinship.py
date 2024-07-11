import itertools

import pytest
import os
import glob
from src.basic.Pedigree import *

accuracy_precision = 0.001


@pytest.fixture
def test_pedigrees():
    return os.path.join(os.path.dirname(__file__), "kinship_test")


@pytest.fixture
def parsed_pedigrees(test_pedigrees):
    pattern = os.path.join(test_pedigrees, "*.pedigree")
    files = glob.glob(pattern)
    pedigrees = {filepath: Pedigree.get_pedigree_graph_from_file(filepath=filepath, separation_symbol=" ",
                                                                 missing_parent_notation=["-1"], skip_first_line=True)
                 for filepath in files}
    return pedigrees


def test_kinship(parsed_pedigrees):
    for filepath, pedigree in parsed_pedigrees.items():
        print(f"Testing kinship calculations on {filepath}")
        pedigree: Pedigree
        (vertex_to_index, direct_kinship_matrix) = pedigree.calculate_kinship()
        kinship_matrix_c_speed = pedigree.calculate_probands_kinship()
        kinship_matrix_c_memory = pedigree.calculate_probands_kinship(mode=KinshipMode.MEMORY)
        for proband in pedigree.get_sink_vertices():
            for other_proband in pedigree.get_sink_vertices():
                simple_algorithm_kinship = direct_kinship_matrix[
                    vertex_to_index[proband], vertex_to_index[other_proband]]
                kinship_c_speed = kinship_matrix_c_speed.get_kinship(proband, other_proband)
                kinship_c_memory = kinship_matrix_c_memory.get_kinship(proband, other_proband)
                for (first, second) in itertools.combinations([simple_algorithm_kinship,
                                                               kinship_c_speed, kinship_c_memory], 2):
                    if abs(first - second) > accuracy_precision:
                        print(f"{proband} {other_proband}: {first} {second}")
                        assert False
