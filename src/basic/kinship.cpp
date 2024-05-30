#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <map>
#include <unordered_map>
#include <queue>
#include <set>
#include <vector>
#include <iostream>

namespace py = pybind11;

typedef std::map<int, std::map<int, float>> KinshipSparseMatrix;

using QueueElement = std::pair<int, std::vector<int>>;

void analyze_cut(KinshipSparseMatrix& matrix)
{
    int count = 0;
    for (const auto& row : matrix)
    {
        count += row.second.size();
    }
    int expected_count = (matrix.size() * (matrix.size() + 1)) / 2;
    std::cout << "Expected count " << expected_count << ", actual count " << count << std::endl;
}

void swap(int& first, int& second)
{
    const int temp = first;
    first = second;
    second = temp;
}

void calculate_self_kinship_sparse(KinshipSparseMatrix& kinship_sparse_matrix, const int vertex,
                                   const std::unordered_map<int, std::vector<int>>& parents)
{
    float kinship = 0;
    auto it = parents.find(vertex);
    if (it == parents.end())
    {
        kinship = 0.5;
    }
    else
    {
        const std::vector<int>& vertex_parents = it->second;
        if (vertex_parents.size() != 2)
        {
            kinship = 0.5;
        } else
        {
            int first_parent = vertex_parents[0];
            int second_parent = vertex_parents[1];
            if (first_parent > second_parent)
            {
                swap(first_parent, second_parent);
            }
            kinship = (1 + kinship_sparse_matrix.at(first_parent).at(second_parent)) / 2.0f;
        }
    }
    kinship_sparse_matrix[vertex][vertex] = kinship;
}

void calculate_pair_kinship_sparse(KinshipSparseMatrix& kinship_sparse_matrix, int first_vertex,
                                    int second_vertex, const std::unordered_map<int, std::vector<int>>& parents)
{
    auto it = parents.find(first_vertex);
    float first_second_kinship_non_normalized = 0.0;
    int second_vertex_copy = second_vertex;
    if (it != parents.end())
    {
        const std::vector<int>& first_parents = it->second;
        if (!first_parents.empty())
        {
            int parent = first_parents[0];
            if (second_vertex_copy > parent)
            {
                swap(parent, second_vertex_copy);
            }
            first_second_kinship_non_normalized += kinship_sparse_matrix.at(second_vertex_copy).at(parent);
        }
        if (first_parents.size() > 1)
        {
            int parent = first_parents[1];
            if (second_vertex_copy > parent)
            {
                swap(parent, second_vertex_copy);
            }
            first_second_kinship_non_normalized += kinship_sparse_matrix.at(second_vertex_copy).at(parent);
        }
        first_second_kinship_non_normalized /= 2.0f;
    }
    if (first_vertex > second_vertex)
    {
        swap(first_vertex, second_vertex);
    }
    kinship_sparse_matrix[first_vertex][second_vertex] = first_second_kinship_non_normalized;
    // kinship_sparse_matrix[second_vertex][first_vertex] = first_second_kinship_non_normalized;
}

KinshipSparseMatrix calculate_kinship_sparse(const std::unordered_set<int>& sink_vertices,
                                             const std::unordered_set<int>& founders,
                                             std::unordered_map<int, int>& parent_to_remaining_children,
                                             std::unordered_map<int, int>& child_to_remaining_parents,
                                             std::unordered_map<int, std::vector<int>>& children,
                                             std::unordered_map<int, std::vector<int>>& parents,
                                             const int order, const int counter_limit)
{
    KinshipSparseMatrix kinship_sparse_matrix;
    std::priority_queue<QueueElement, std::vector<QueueElement>, std::greater<QueueElement>> queue;

    for (const auto& founder : founders)
    {
        queue.push(std::make_pair(1, std::vector<int>{founder}));
    }

    int processed_vertices = 0;
    int counter = 0;

    while (!queue.empty())
    {
        auto [_, vertices] = queue.top();
        queue.pop();
        for (const auto& vertex : vertices)
        {
            counter++;
            if (counter == counter_limit)
            {
                std::cout << "The size of the cut: " << kinship_sparse_matrix.size() << std::endl;
                analyze_cut(kinship_sparse_matrix);
                // std::cout << "The number of buckets: " << kinship_sparse_matrix.bucket_count() << std::endl;
                std::cout << "Queue size: " << queue.size() << std::endl;
                std::cout << "Progress " << static_cast<float>(processed_vertices) / order << std::endl;
                counter = 0;
            }
            processed_vertices++;

            // Calculate self kinship
            calculate_self_kinship_sparse(kinship_sparse_matrix, vertex, parents);

            for (const auto &[processed_vertex, _] : kinship_sparse_matrix)
            {
                if (processed_vertex == vertex)
                {
                    continue;
                }
                calculate_pair_kinship_sparse(kinship_sparse_matrix, vertex, processed_vertex, parents);
            }

            for (const auto& parent : parents.at(vertex))
            {
                parent_to_remaining_children[parent]--;
                if (parent_to_remaining_children[parent] == 0 && sink_vertices.find(parent) == sink_vertices.end())
                {
                    parent_to_remaining_children.erase(parent);
                    kinship_sparse_matrix.erase(parent);
                    for (auto& other_vertex : kinship_sparse_matrix)
                    {
                        other_vertex.second.erase(parent);
                    }
                }
            }

            std::unordered_set<int> children_to_add;
            for (const auto& child : children.at(vertex))
            {
                child_to_remaining_parents[child]--;
                if (child_to_remaining_parents[child] == 0)
                {
                    children_to_add.insert(child);
                    child_to_remaining_parents.erase(child);
                }
            }

            if (!children_to_add.empty())
            {
                float additional_space = children_to_add.size();
                std::unordered_set<int> children_parents;
                for (const auto &child : children_to_add)
                {
                    for (const auto& parent : parents.at(child))
                    {
                        children_parents.insert(parent);
                    }
                }
                for (const auto &child_parent : children_parents)
                {
                    bool to_be_removed = true;
                    for (const auto& child : children.at(child_parent))
                    {
                        if (children_to_add.find(child) == children_to_add.end() &&
                            parent_to_remaining_children.find(child) == parent_to_remaining_children.end() &&
                            kinship_sparse_matrix.find(child) == kinship_sparse_matrix.end())
                        {
                            to_be_removed = false;
                            break;
                        }
                    }
                    if (to_be_removed)
                    {
                        additional_space -= 0.5f;
                    }
                }
                queue.push(std::make_pair(additional_space, std::vector<int>(children_to_add.begin(), children_to_add.end())));
            }
        }
    }
    return kinship_sparse_matrix;
}

PYBIND11_MODULE(kinship, m)
{
    m.def("calculate_kinship_sparse", &calculate_kinship_sparse, "Calculate kinship sparse matrix",
          py::arg("sink_vertices"), py::arg("founders"), py::arg("parent_to_remaining_children"),
          py::arg("child_to_remaining_parents"), py::arg("children"), py::arg("parents"),
          py::arg("order"), py::arg("counter_limit"));
}
