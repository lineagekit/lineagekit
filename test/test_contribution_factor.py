import os

import numpy
import pytest

from lineagekit.core.pedigree import Pedigree


@pytest.fixture
def test_pedigree():
    return os.path.join(os.path.dirname(__file__), "test_data", "1000_8.pedigree")


def test_contribution_factor_calculation(test_pedigree):
    pedigree = Pedigree.get_pedigree_graph_from_file(filepath=test_pedigree, separation_symbol=" ",
                                                     missing_parent_notation=["-1"], skip_first_line=True)
    founders = pedigree.get_founders()
    founders_contribution_factors = pedigree.calculate_contribution_factors()
    for founder in founders:
        founder_factor = pedigree.calculate_vertex_contribution_factor(founder)
        assert numpy.isclose(founder_factor, founders_contribution_factors[founder])
