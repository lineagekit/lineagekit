import itertools
from typing import Iterable

import networkx
import pytest
import os
from basic.GenGraph import GenGraph
from basic.Pedigree import Pedigree
from basic.PloidPedigree import PloidPedigree
from basic.CoalescentTree import CoalescentTree

from src.basic.AbstractPedigree import AbstractPedigree


@pytest.fixture
def test_pedigrees():
    return os.path.join(os.path.dirname(__file__), "test_pedigrees")


@pytest.fixture()
def simple_1_missing_parent_notation():
    return ["-1", "."]


@pytest.fixture
def parse_simple_1(test_pedigrees, simple_1_missing_parent_notation) -> PloidPedigree:
    filepath = f"{test_pedigrees}/simple_1.txt"

    return PloidPedigree.get_ploid_pedigree_from_file(filepath=filepath,
                                                      missing_parent_notation=simple_1_missing_parent_notation)


@pytest.fixture
def parse_simple_1_haploid(test_pedigrees, simple_1_missing_parent_notation) -> Pedigree:
    filepath = f"{test_pedigrees}/simple_1.txt"

    return Pedigree.get_pedigree_graph_from_file(filepath=filepath,
                                                 missing_parent_notation=simple_1_missing_parent_notation)


@pytest.fixture
def parse_simple_1_haploid_ascending_proband() -> Iterable[int]:
    return 1, 3, -123123


@pytest.fixture
def parse_simple_1_haploid_ascending(test_pedigrees, parse_simple_1_haploid_ascending_proband,
                                     simple_1_missing_parent_notation) -> Pedigree:
    filepath = f"{test_pedigrees}/simple_1.txt"

    return Pedigree.get_pedigree_graph_from_file(filepath=filepath,
                                                 probands=parse_simple_1_haploid_ascending_proband,
                                                 missing_parent_notation=simple_1_missing_parent_notation)


@pytest.fixture
def simple_1(parse_simple_1) -> PloidPedigree:
    return parse_simple_1.copy()


@pytest.fixture
def simple_1_haploid(parse_simple_1_haploid) -> Pedigree:
    return parse_simple_1_haploid.copy()


@pytest.fixture
def simple_1_haploid_ascending(parse_simple_1_haploid_ascending) -> Pedigree:
    return parse_simple_1_haploid_ascending.copy()


@pytest.fixture
def parse_coalescent_tree_1(test_pedigrees) -> CoalescentTree:
    filepath = f"{test_pedigrees}/coalescent_tree_1.txt"
    return CoalescentTree.get_coalescent_tree_from_file(filepath=filepath)


@pytest.fixture
def coalescent_tree_1(parse_coalescent_tree_1) -> CoalescentTree:
    return parse_coalescent_tree_1.copy()


@pytest.fixture
def parse_coalescent_tree_2(test_pedigrees) -> CoalescentTree:
    filepath = f"{test_pedigrees}/coalescent_tree_2.txt"
    return CoalescentTree.get_coalescent_tree_from_file(filepath=filepath)


@pytest.fixture
def coalescent_tree_2(parse_coalescent_tree_2) -> CoalescentTree:
    return parse_coalescent_tree_2.copy()


def two_parents_defined_correctly(graph: GenGraph, child: int, first_parent: int, second_parent: int) -> bool:
    child_ploid_set = {2 * child, 2 * child + 1}
    if not set(graph.get_parents(2 * child)).union(graph.get_parents(2 * child + 1)) == {2 * first_parent,
                                                                                         2 * first_parent + 1,
                                                                                         2 * second_parent,
                                                                                         2 * second_parent + 1}:
        return False

    def ploid_children_map_is_consistent(individual_id):
        return set(graph.get_children(2 * individual_id)) == set(graph.get_children(2 * individual_id + 1))

    if not ploid_children_map_is_consistent(first_parent) or not ploid_children_map_is_consistent(second_parent):
        return False
    for (first_ploid, second_ploid) in itertools.product([2 * first_parent, 2 * first_parent + 1],
                                                         [2 * second_parent, 2 * second_parent + 1]):
        if not child_ploid_set.issubset(
                set(graph.get_children(first_ploid)).union(graph.get_children(second_ploid))):
            return False
    return True


def diploid_child_parent_relationship_consistent(graph: GenGraph, child: int, parent: int):
    if set(graph.get_children(2 * parent)) != set(graph.get_children(2 * parent + 1)):
        return False
    parent_ploid_set = {2 * parent, 2 * parent + 1}
    ploids_related_with_parent = 0
    for child_ploid in [2 * child, 2 * child + 1]:
        if parent_ploid_set == set(graph.get_parents(child_ploid)):
            ploids_related_with_parent += 1
            if child_ploid not in graph.get_children(2 * parent):
                return False
    return ploids_related_with_parent == 1


def transform_individual_ids_into_ploids(values: Iterable[int]) -> set[int]:
    return {val for x in values for val in [2 * x, 2 * x + 1]}


def test_subgraph_consistency(simple_1):
    simple_1.remove_nodes_from(list(simple_1.nodes()))
    assert simple_1.get_vertices_number() == 0, (f"Expected the graph to be empty after removing all the vertices,"
                                                 f" got {simple_1.get_vertices_number()} vertices instead")


def test_ascending_genealogy_does_not_throw_exceptions_for_non_existing_vertices(simple_1_haploid):
    ascending_genealogy = simple_1_haploid.get_ascending_vertices_from_probands([-1, -12123, 1, 2])
    assert set(ascending_genealogy) == {1, 2, 6, 7, 10, 11}


def test_diploid_graph_parsing_and_basic_functions(test_pedigrees, simple_1, simple_1_missing_parent_notation):
    filepath = f"{test_pedigrees}/simple_1.txt"
    missing_parent_notation = simple_1_missing_parent_notation
    graph = simple_1
    assert graph.verify_max_parents_number(2)
    assert not graph.verify_max_parents_number(1)
    lines = open(filepath, 'r').readlines()
    lines.pop(0)
    for line in lines:
        [child, first_parent, second_parent] = line.strip().split(' ')
        child = int(child)
        missing_parents = 0
        for parent in [first_parent, second_parent]:
            if parent in missing_parent_notation:
                missing_parents += 1
            else:
                parent = int(parent)
                assert diploid_child_parent_relationship_consistent(graph, child, parent)
        assert missing_parents == len([x for x in [2 * child, 2 * child + 1] if not graph.get_parents(x)])
    assert graph.get_vertices_number() == 22, f"Expected 22 vertices, got {graph.get_vertices_number()}"
    assert set(graph.nodes()) == set(range(2, 24)), (f"Expected nodes' ids to be from 1 to 22, "
                                                     f"got {set(graph.nodes())}")
    ascending_graph = graph.get_ascending_vertices_from_probands([4, 5])
    ploids_ascending_graph = transform_individual_ids_into_ploids([2, 6, 7, 10, 11])
    assert ascending_graph == ploids_ascending_graph, (f"Expected {ploids_ascending_graph},"
                                                       f" got {ascending_graph}")

    graph.remove_nodes_from(transform_individual_ids_into_ploids([8]))
    ascending_graph = graph.get_ascending_vertices_from_probands(transform_individual_ids_into_ploids([3, 4, 5]))
    ascending_graph_correct = transform_individual_ids_into_ploids([3, 4, 5, 9])
    assert ascending_graph == ascending_graph_correct, (f"Expected {ascending_graph_correct},"
                                                        f" got {ascending_graph}")
    graph.add_edges_from([(12, 6), (13, 6), (14, 9), (15, 9)])
    ascending_graph = graph.get_ascending_vertices_from_probands(transform_individual_ids_into_ploids([3, 4, 5]))
    ascending_graph_correct = transform_individual_ids_into_ploids([3, 4, 5, 6, 7, 9, 10, 11])
    assert ascending_graph == ascending_graph_correct, (f"Expected {ascending_graph_correct},"
                                                        f" got {ascending_graph}")
    graph.add_edges_from([(6, 2), (7, 3), (8, 4), (9, 4)])
    ascending_graph = graph.get_ascending_vertices_from_probands(transform_individual_ids_into_ploids([3, 4, 5]))
    assert ascending_graph == ascending_graph_correct, (f"Expected {ascending_graph_correct},"
                                                        f" got {ascending_graph}")
    ascending_graph = graph.get_ascending_vertices_from_probands(transform_individual_ids_into_ploids([1, 2]))
    ascending_graph_correct = transform_individual_ids_into_ploids([1, 2, 3, 4, 6, 7, 9, 10, 11])
    assert ascending_graph == ascending_graph_correct, (f"Expected {ascending_graph_correct},"
                                                        f" got {ascending_graph}")
    sink_vertices = set(graph.get_sink_vertices())
    sink_vertices_valid = transform_individual_ids_into_ploids([1, 2, 5])
    assert sink_vertices_valid == sink_vertices, (f"Expected sink vertices to be {sink_vertices_valid},"
                                                  f" got {sink_vertices}")
    for sink_vertex in graph.get_sink_vertices():
        assert not graph.get_children(sink_vertex)
        assert not graph.has_children(sink_vertex)
    graph_nodes = set(graph.nodes())
    for non_sink_vertex in graph_nodes.difference(sink_vertices):
        assert graph.get_children(non_sink_vertex)
        assert graph.has_children(non_sink_vertex)
    orphans = set(graph.get_orphans())
    # One of the ploid for 4 and 5 have no parents.Since we don't want to rely on the implementation,
    # we will verify it this way
    orphans_superset = transform_individual_ids_into_ploids([4, 5, 7, 9, 10, 11])
    assert orphans_superset.issuperset(orphans), f"Expected orphans to be {orphans_superset} to contain {orphans}"
    assert len(transform_individual_ids_into_ploids([4, 5]).intersection(orphans)) == 2, (f"Expected only two ploids "
                                                                                          f"from 4 and 5 not "
                                                                                          f"to have parents")
    for orphan in orphans:
        assert not graph.get_parents(orphan)
        assert not graph.has_parents(orphan)
        assert graph.is_orphan(orphan)
    for non_orphan in graph_nodes.difference(orphans):
        assert graph.get_parents(non_orphan)
        assert graph.has_parents(non_orphan)
        assert not graph.is_orphan(non_orphan)


def test_ascending_genealogy_reduction(simple_1_haploid):
    graph = simple_1_haploid
    graph.reduce_to_ascending_graph([1, 3])
    assert set(graph.nodes()) == {1, 6, 7, 10, 11, 3, 8, 9}
    assert graph.get_children(8) == [3]
    graph.reduce_to_ascending_graph([1])
    assert set(graph.nodes()) == {1, 6, 7, 10, 11}


def test_ascending_genealogy_parsing(simple_1_haploid, simple_1_haploid_ascending,
                                     parse_simple_1_haploid_ascending_proband):
    graph_reduced = simple_1_haploid_ascending
    graph_full = simple_1_haploid
    graph_full.reduce_to_ascending_graph(parse_simple_1_haploid_ascending_proband)
    assert networkx.is_isomorphic(graph_reduced, graph_full)


def test_connected_components(simple_1_haploid):
    graph = simple_1_haploid
    graph.remove_node(6)
    connected_component_first = set(graph.get_connected_component_for_vertex(8))
    connected_component_valid = {3, 4, 5, 8, 9, 10, 11}
    assert connected_component_first == connected_component_valid, (f"Expected the connected component to be"
                                                                    f" {connected_component_valid},"
                                                                    f"got {connected_component_first}")
    connected_component_second = set(graph.get_connected_component_for_vertex(1))
    connected_component_valid = {1, 2, 7}
    assert connected_component_second == connected_component_valid, (f"Expected the connected component to be"
                                                                     f" {connected_component_valid},"
                                                                     f"got {connected_component_valid}")
    connected_components = [connected_component_first, connected_component_second]
    connected_components_valid = graph.get_connected_components()
    # The connected components are of different sizes, we can use the length as the key
    assert sorted(connected_components_valid, key=len) == sorted(connected_components, key=len), \
        (f"Expected the connected components to be "
         f"{connected_components_valid}, "
         f"got {connected_components} instead")


def haploid_child_parent_relationship_consistent(graph: GenGraph, child: int, parent: int):
    return parent in graph.get_parents(child) and child in graph.get_children(parent) \
        and graph.is_parent(parent=parent, child=child)


def test_ploid_parsing(test_pedigrees, simple_1_missing_parent_notation, simple_1_haploid):
    assert simple_1_haploid.get_vertices_number() == 11, (f"Expected 11 vertices in the graph, "
                                                          f"got {simple_1_haploid.get_vertices_number()}")
    assert set(simple_1_haploid.nodes()) == set(range(1, 12)), (f"Expected vertices from 1 to 11, "
                                                                f"got {simple_1_haploid.nodes()}")
    filepath = f"{test_pedigrees}/simple_1.txt"
    missing_parent_notation = simple_1_missing_parent_notation
    graph = simple_1_haploid
    lines = open(filepath, 'r').readlines()
    lines.pop(0)
    for line in lines:
        [child, first_parent, second_parent] = line.strip().split(' ')
        child = int(child)
        missing_parents = 0
        for parent in [first_parent, second_parent]:
            if parent in missing_parent_notation:
                missing_parents += 1
            else:
                parent = int(parent)
                assert haploid_child_parent_relationship_consistent(graph, child, parent)
        assert 2 - missing_parents == len(graph.get_parents(child)), (f"Expected {2 - missing_parents} parents "
                                                                      f"for {child},"
                                                                      f"got {len(graph.get_parents(child))} instead")


def test_levels_assignment(simple_1_haploid):
    graph = simple_1_haploid
    valid_levels = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 1, 7: 1, 8: 1, 9: 1, 10: 2, 11: 2}
    for vertex, level in valid_levels.items():
        assert vertex in graph
        assert graph.get_vertex_level(vertex) == level
    assert graph.is_founder(10)
    assert graph.is_founder(11)
    graph.remove_edges_to_parents(6)
    graph.remove_edges_to_children(6)
    graph.remove_edges_to_children(8)
    graph.remove_edges_to_parents(8)
    valid_levels = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 8: 0, 10: 0, 11: 0, 7: 1, 9: 1}
    for vertex, level in valid_levels.items():
        assert vertex in graph
        assert graph.get_vertex_level(vertex) == level
    assert graph.is_founder(7)
    assert graph.is_founder(9)


def test_saving_to_file(simple_1_haploid, simple_1_missing_parent_notation):
    graph = simple_1_haploid
    filepath = "test.txt"
    graph.save_to_file("test.txt", separator=" ;", missing_parent_notation="!")
    assert os.path.exists(filepath)
    try:
        parsed_graph = Pedigree.get_graph_from_file(filepath=filepath,
                                                    separation_symbol=" ;",
                                                    missing_parent_notation=["!"])
    except Exception as ex:
        os.remove(filepath)
        pytest.fail(f"Could not parse the saved file, the exception: {ex}")

    assert networkx.is_isomorphic(graph, parsed_graph)
    os.remove(filepath)


def test_saving_as_diploid(simple_1):
    filepath = "test.txt"
    simple_1.save_as_diploid(filepath=filepath)
    try:
        parsed_simple_1 = PloidPedigree.get_ploid_pedigree_from_file(filepath=filepath)
    except Exception as ex:
        os.remove(filepath)
        pytest.fail(f"Could not parse the saved file, the exception: {ex}")
    assert networkx.is_isomorphic(simple_1, parsed_simple_1)
    os.remove(filepath)
    individuals_ids = frozenset([1, 3, 4])
    ploid_ids = PloidPedigree.get_ploids_from_individual_ids(individuals_ids)
    assert {2, 3, 6, 7, 8, 9} == frozenset(ploid_ids)
    assert frozenset(PloidPedigree.get_individual_ids_from_ploids(ploid_ids)) == individuals_ids
    simple_1.save_ascending_genealogy_as_diploid(filepath=filepath, vertices=ploid_ids)
    simple_1.reduce_to_ascending_graph(ploid_ids)
    try:
        parsed_graph = PloidPedigree.get_ploid_pedigree_from_file(filepath=filepath)
    except Exception as ex:
        os.remove(filepath)
        pytest.fail(f"Could not parse the saved file, the exception: {ex}")
    assert networkx.is_isomorphic(simple_1, parsed_graph)
    os.remove(filepath)


@pytest.mark.parametrize("data", ["parse_simple_1", "parse_simple_1_haploid",
                                  "simple_1", "simple_1_haploid"])
@pytest.mark.parametrize("run", range(30))
def test_error_simulation(data, run, request):
    pedigree: AbstractPedigree = request.getfixturevalue(data)
    simulated_errors = pedigree.introduce_and_record_errors(error_rate=0.5, apply_errors=False)

    def verify_edges(child_vertex, disconnected_vertices, connected_vertices):
        for disconnected_vertex in disconnected_vertices:
            assert not pedigree.has_edge(child=child_vertex, parent=disconnected_vertex), \
                (f"The edge {child_vertex}—{disconnected_vertex} is supposed to be removed from the graph,"
                 f" but it's present")
        for connected_vertex in connected_vertices:
            assert pedigree.has_edge(child=child_vertex, parent=connected_vertex), \
                (f"The edge {child_vertex}—{connected_vertex} "
                 f"is supposed to be present in the graph, but it's not")

    for vertex, removed_vertices, add_vertices in simulated_errors:
        verify_edges(vertex, add_vertices, removed_vertices)
        # Verify that the new parent is not an old parent
        assert not [x for x in add_vertices if x in pedigree.get_parents(vertex)]
        # Verify that the vertices are on the same level
        for add_vertex, remove_vertex in itertools.product(add_vertices, removed_vertices):
            assert pedigree.get_vertex_level(add_vertex) == pedigree.get_vertex_level(remove_vertex)
    pedigree.apply_errors(simulated_errors)
    for vertex, removed_vertices, add_vertices in simulated_errors:
        verify_edges(vertex, removed_vertices, add_vertices)
    pedigree.reverse_errors(simulated_errors)
    for vertex, removed_vertices, add_vertices in simulated_errors:
        verify_edges(vertex, add_vertices, removed_vertices)


def test_coalescent_tree_parsing(coalescent_tree_1):
    assert len(coalescent_tree_1.nodes) == 8
    assert coalescent_tree_1.get_root_for_clade(coalescent_tree_1.nodes) == 8
    copy_tree = coalescent_tree_1.copy()
    coalescent_tree_1.remove_unary_nodes()
    assert networkx.is_isomorphic(copy_tree, coalescent_tree_1)
    assert (frozenset(coalescent_tree_1.get_largest_clade_by_size()) ==
            frozenset(coalescent_tree_1.get_largest_clade_by_probands()))


def test_coalescent_tree_2(coalescent_tree_2):
    largest_clade_by_size = coalescent_tree_2.get_largest_clade_by_size()
    assert frozenset(largest_clade_by_size) == frozenset(range(1, 9))
    largest_clade_by_probands = coalescent_tree_2.get_largest_clade_by_probands()
    assert frozenset(largest_clade_by_probands) == frozenset(range(9, 13))
    assert coalescent_tree_2.get_root_for_clade(largest_clade_by_size) == 8
    assert coalescent_tree_2.get_root_for_clade(largest_clade_by_probands) == 12
    coalescent_tree_2.remove_unary_nodes()
    assert coalescent_tree_2.get_parents(1) == [8]
    assert coalescent_tree_2.get_parents(2) == [8]
    for removed_vertex in range(3, 8):
        assert removed_vertex not in coalescent_tree_2
    largest_clade_by_size_new = coalescent_tree_2.get_largest_clade_by_size()
    assert frozenset(largest_clade_by_size_new) == frozenset(largest_clade_by_probands)
