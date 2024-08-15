#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <queue>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include "parallel_hashmap/phmap.h"
#include "sparsepp/spp.h"

namespace py = pybind11;

using phmap::flat_hash_map;
using phmap::parallel_flat_hash_map;
using spp::sparse_hash_map;
using phmap::flat_hash_set;

typedef sparse_hash_map<int, sparse_hash_map<int, float>> MemorySparseMatrix;
typedef flat_hash_map<int, flat_hash_map<int, float>> TimeSparseMatrix;

template<typename T>
struct inner_map_type;

template<typename K, typename V>
struct inner_map_type<sparse_hash_map<K, sparse_hash_map<K, V>>> {
    using type = sparse_hash_map<K, V>;
};

template<typename K, typename V>
struct inner_map_type<flat_hash_map<K, flat_hash_map<K, V>>> {
    using type = flat_hash_map<K, V>;
};

template<typename KinshipMatrix>
void calculate_pair_kinships_sparse(KinshipMatrix& kinship_sparse_matrix, const int vertex,
    const std::unordered_map<int, std::vector<int>>& parents_map)
{
    using inner_map = typename inner_map_type<KinshipMatrix>::type;
    auto it = parents_map.find(vertex);
    inner_map& vertex_map = kinship_sparse_matrix[vertex];
    inner_map* first_parent_map = nullptr;
    inner_map* second_parent_map = nullptr;

    int first_parent = -1;
    int second_parent = -1;
    float self_kinship = 0.5f;
    if (it != parents_map.end())
    {
        const std::vector<int>& parents = it->second;
        if (!parents.empty())
        {
            first_parent = parents[0];
            first_parent_map = &kinship_sparse_matrix.at(first_parent);
            if (parents.size() > 1)
            {
                second_parent = parents[1];
                second_parent_map = &kinship_sparse_matrix.at(second_parent);
                if (first_parent > second_parent)
                {
                    self_kinship = (1 + kinship_sparse_matrix.at(second_parent).at(first_parent)) / 2.0f;
                }
                else
                {
                    self_kinship = (1 + kinship_sparse_matrix.at(first_parent).at(second_parent)) / 2.0f;
                }
            }
        }
    }
    vertex_map[vertex] = self_kinship;
    for (auto& [second_vertex, second_vertex_map] : kinship_sparse_matrix)
    {
        if (second_vertex == vertex)
        {
            continue;
        }
        float first_second_kinship_non_normalized = 0.0;
        if (first_parent_map)
        {
            if (second_vertex > first_parent)
            {
                first_second_kinship_non_normalized += first_parent_map->at(second_vertex);
            }
            else
            {
                first_second_kinship_non_normalized += second_vertex_map.at(first_parent);
            }
            if (second_parent_map)
            {
                if (second_vertex > second_parent)
                {
                    first_second_kinship_non_normalized += second_parent_map->at(second_vertex);
                }
                else
                {
                    first_second_kinship_non_normalized += second_vertex_map.at(second_parent);
                }
            }
            first_second_kinship_non_normalized /= 2.0f;
        }
        if (vertex > second_vertex)
        {
            second_vertex_map[vertex] = first_second_kinship_non_normalized;
        }
        else
        {
            vertex_map[second_vertex] = first_second_kinship_non_normalized;
        }
    }
}

using QueueElement = std::pair<float, std::vector<int>>;

struct CompareOnlyFirst {
    bool operator()(const QueueElement& a, const QueueElement& b) {
        return a.first > b.first; // Min-heap: smallest priority at the top
    }
};

template<typename KinshipMatrix>
KinshipMatrix calculate_kinship_sparse(
    std::unordered_map<int, std::vector<int>>& children,
    std::unordered_map<int, std::vector<int>>& parents,
    std::unordered_set<int>& sink_vertices)
{
    // Find the founders
    flat_hash_set<int> founders;
    for (const auto& entry : parents)
    {
        if (entry.second.empty())
        {
            founders.insert(entry.first);
        }
    }
    // Initialize the maps that will keep track of how many children/parents need to be processed
    flat_hash_map<int, int> parent_to_remaining_children;
    flat_hash_map<int, int> child_to_remaining_parents;
    parent_to_remaining_children.reserve(children.size());
    child_to_remaining_parents.reserve(parents.size());
    for (const auto& entry : children)
    {
        parent_to_remaining_children[entry.first] = entry.second.size();
    }
    for (const auto& entry : parents)
    {
        child_to_remaining_parents[entry.first] = entry.second.size();
    }
    // Initialize the queue with the founders
    KinshipMatrix kinship_sparse_matrix;
    std::priority_queue<QueueElement, std::vector<QueueElement>, CompareOnlyFirst> queue;
    for (const auto& founder : founders)
    {
        queue.emplace(1, std::vector<int>{founder});
    }
    // Calculate the kinships
    while (!queue.empty())
    {
        auto [priority, vertices] = queue.top();
        queue.pop();
        for (const int vertex : vertices)
        {
            calculate_pair_kinships_sparse(kinship_sparse_matrix, vertex, parents);
            for (const auto& parent : parents.at(vertex))
            {
                // Updating the counter of unprocessed children for the parent
                if (--parent_to_remaining_children[parent] == 0 && sink_vertices.find(parent) == sink_vertices.end())
                {
                    // Erasing the parent's kinship information because it isn't a proband
                    // and all of its children have been processed
                    parent_to_remaining_children.erase(parent);
                    kinship_sparse_matrix.erase(parent);
                    for (auto& other_vertex : kinship_sparse_matrix)
                    {
                        if (parent > other_vertex.first)
                        {
                            other_vertex.second.erase(parent);
                        }
                    }
                }
            }

            flat_hash_set<int> children_to_add;
            for (const auto& child : children.at(vertex))
            {
                if (--child_to_remaining_parents[child] == 0)
                {
                    children_to_add.insert(child);
                    child_to_remaining_parents.erase(child);
                }
            }

            if (!children_to_add.empty())
            {
                float additional_space = (float) children_to_add.size();
                flat_hash_set<int> children_parents;
                for (const int child : children_to_add)
                {
                    for (const int parent : parents.at(child))
                    {
                        children_parents.insert(parent);
                    }
                }
                for (const int child_parent : children_parents)
                {
                    const int children_unprocessed = parent_to_remaining_children.at(child_parent);
                    int counter = 0;
                    for (const int child : children.at(child_parent))
                    {
                        if (children_to_add.find(child) == children_to_add.end())
                        {
                            counter++;
                        }
                    }
                    if (children_unprocessed != counter)
                    {
                        additional_space -= 1.0f;
                    }
                }
                queue.emplace(additional_space, std::vector<int>(children_to_add.begin(), children_to_add.end()));
            }
        }
    }
    return kinship_sparse_matrix;
}

TimeSparseMatrix calculate_kinship_sparse_speed(
    std::unordered_map<int, std::vector<int>>& children,
    std::unordered_map<int, std::vector<int>>& parents,
    std::unordered_set<int>& sink_vertices)
{
    return calculate_kinship_sparse<TimeSparseMatrix>(children, parents, sink_vertices);
}

MemorySparseMatrix calculate_kinship_sparse_memory(
    std::unordered_map<int, std::vector<int>>& children,
    std::unordered_map<int, std::vector<int>>& parents,
    std::unordered_set<int>& sink_vertices)
{
    return calculate_kinship_sparse<MemorySparseMatrix>(children, parents, sink_vertices);
}

template<typename KinshipMatrix>
std::pair<std::unordered_map<size_t, size_t>, py::array_t<float>> convert_to_numpy_and_free(KinshipMatrix& self)
{
    // Create a mapping of keys to contiguous indices
    std::unordered_map<size_t, size_t> key_to_index;
    int index = 0;
    for (const auto& row : self)
    {
        key_to_index[row.first] = index++;
    }

    size_t size = key_to_index.size();

    // Create a NumPy array of shape (size, size)
    py::array_t<float> numpy_matrix({size, size});
    auto buffer = numpy_matrix.mutable_unchecked<2>();

    for (const auto& row : self)
    {
        int row_index = key_to_index[row.first];
        for (const auto& col : row.second)
        {
            int col_index = key_to_index[col.first];
            buffer(row_index, col_index) = col.second;
            buffer(col_index, row_index) = col.second;
        }
    }
    self.clear();
    return {key_to_index, numpy_matrix};
}


PYBIND11_MODULE(kinship, m)
{
    py::class_<TimeSparseMatrix>(m, "TimeSparseMatrix")
        .def("get_kinship", [](const TimeSparseMatrix& self, int key1, int key2)
        {
            if (key1 < key2)
            {
                return self.at(key1).at(key2);
            }
            return self.at(key2).at(key1);
        }, "Get the float value for two integers", py::arg("key1"), py::arg("key2"))
        .def("to_numpy_and_free", [](TimeSparseMatrix& self)
        {
            auto [key_to_index, numpy_matrix] = convert_to_numpy_and_free(self);
            return std::make_tuple(key_to_index, numpy_matrix);
        }, "Convert the sparse matrix to a NumPy array and free the memory");
    py::class_<MemorySparseMatrix>(m, "MemorySparseMatrix")
        .def("get_kinship", [](const MemorySparseMatrix& self, int key1, int key2)
        {
            if (key1 < key2)
            {
                return self.at(key1).at(key2);
            }
            return self.at(key2).at(key1);
        }, "Get the float value for two integers", py::arg("key1"), py::arg("key2"))
         .def("to_numpy_and_free", [](MemorySparseMatrix& self)
        {
            auto [key_to_index, numpy_matrix] = convert_to_numpy_and_free(self);
            return std::make_tuple(key_to_index, numpy_matrix);
        }, "Convert the sparse matrix to a NumPy array and free the memory");
    m.def("calculate_kinship_sparse_speed", &calculate_kinship_sparse_speed,
          "Calculate kinship sparse matrix (running time preference)",
          py::arg("children"), py::arg("parents"), py::arg("sink_vertices"));
    m.def("calculate_kinship_sparse_memory", &calculate_kinship_sparse_memory,
          "Calculate kinship sparse matrix (memory preference)",
          py::arg("children"), py::arg("parents"), py::arg("sink_vertices"));
}
