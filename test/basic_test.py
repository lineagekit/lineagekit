import pytest
import os
from src.basic.simple_graph import SimpleGraph


@pytest.fixture
def test_pedigrees():
    return os.path.join(os.path.dirname(__file__), "test_pedigrees")


def test_simple_graph(test_pedigrees):
    filepath = f"{test_pedigrees}/simple_1.txt"
    graph = SimpleGraph.get_diploid_graph_from_file(filename=filepath)
    assert graph.get_vertices_number() == 22, f"Expected 22 vertices, got {graph.get_vertices_number()}"
    ascending_graph = graph.get_ascending_graph_from_vertices([4, 5])
    assert len(ascending_graph) == 10, f"Expected 10 vertices, got {len(ascending_graph)}"
