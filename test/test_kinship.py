import itertools

import pytest
import os
import glob
from basic.Pedigree import *

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


@pytest.fixture
def kinship_test_2_path(test_pedigrees):
    return os.path.join(test_pedigrees, "kinship_test_2.pedigree")


@pytest.fixture
def kinship_test_2(kinship_test_2_path):
    return Pedigree.get_pedigree_graph_from_file(filepath=kinship_test_2_path, separation_symbol=" ",
                                                 missing_parent_notation=["-1"], skip_first_line=True)


def test_kinship_custom_probands(kinship_test_2):
    kinship_test_2: Pedigree
    (vertex_to_index, direct_kinship_matrix) = kinship_test_2.calculate_kinship()
    kinship_direct = direct_kinship_matrix[vertex_to_index[1], vertex_to_index[3]]
    kinship_matrix_c_speed = kinship_test_2.calculate_probands_kinship(probands={1, 3})
    kinship_matrix_c_memory = kinship_test_2.calculate_probands_kinship(probands={1, 3}, mode=KinshipMode.MEMORY)
    kinship_speed = kinship_matrix_c_speed.get_kinship(1, 3)
    kinship_matrix = kinship_matrix_c_memory.get_kinship(1, 3)
    assert abs(kinship_direct - kinship_speed) <= accuracy_precision
    assert abs(kinship_direct - kinship_matrix) <= accuracy_precision
    assert abs(kinship_speed - kinship_matrix) <= accuracy_precision
    # Testing that only the passed probands are present in the resulting map
    for first, second in itertools.combinations(set(kinship_test_2.get_sink_vertices()).difference([1, 3]), 2):
        with pytest.raises(IndexError):
            kinship_matrix_c_speed.get_kinship(first, second)
        with pytest.raises(IndexError):
            kinship_matrix_c_speed.get_kinship(first, second)


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
