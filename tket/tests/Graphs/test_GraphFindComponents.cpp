// Copyright 2019-2022 Cambridge Quantum Computing
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <catch2/catch_test_macros.hpp>
#include <set>

#include "Graphs/AdjacencyData.hpp"
#include "Graphs/GraphRoutines.hpp"
#include "Utils/RNG.hpp"

using std::map;
using std::set;
using std::vector;

namespace tket {
namespace graphs {
namespace tests {

// For testing the connected component function
struct ComponentsTestData {
  vector<vector<size_t>> raw_adjacency_data;
  vector<set<size_t>> components;
};

// Used to test the connected component functions:
// create a graph with known components
struct ComponentsParameters {
  size_t max_tree_size = 10;
  size_t max_number_of_extra_edges_to_add_to_tree = 2;
  size_t number_of_singletons = 10;
  size_t number_of_multivertex_components = 10;

  ComponentsTestData get_test_data(RNG& rng) const {
    ComponentsTestData data;

    // KEY: the vertex VALUE: the neighbours
    // (simpler to use a map rather than a vector).
    map<size_t, vector<size_t>> neighbours;

    size_t next_vertex = 0;

    // Singletons - no edges
    for (; next_vertex < number_of_singletons; ++next_vertex) {
      data.components.emplace_back(set<size_t>{next_vertex});
    }

    // Trees
    set<size_t> tree_vertices;

    // Each vertex is a tree node which may grow some children
    vector<size_t> tree_vertex_indices_to_grow;

    for (size_t i = 0; i < number_of_multivertex_components; ++i) {
      // The initial seed node.
      tree_vertices.insert(next_vertex);
      tree_vertex_indices_to_grow.push_back(next_vertex);
      ++next_vertex;

      // Keep growing a random tree node.
      while (tree_vertices.size() < max_tree_size) {
        const size_t node = rng.get_element(tree_vertex_indices_to_grow);
        tree_vertex_indices_to_grow.push_back(next_vertex);
        tree_vertices.insert(next_vertex);
        neighbours[node].emplace_back(next_vertex);
        ++next_vertex;
      }
      // Now, we've made a tree.
      data.components.emplace_back(tree_vertices);
      tree_vertex_indices_to_grow.clear();

      // Should be contiguous vertex indices.
      REQUIRE(!tree_vertices.empty());
      const size_t min_index = *tree_vertices.cbegin();
      const size_t max_index = *tree_vertices.crbegin();
      REQUIRE(tree_vertices.size() == max_index - min_index + 1);
      REQUIRE(max_index + 1 == next_vertex);
      tree_vertices.clear();

      // Add a few extra edges to the tree
      for (size_t i = 0; i < max_number_of_extra_edges_to_add_to_tree; ++i) {
        const size_t vertex1 = rng.get_size_t(min_index, max_index);
        const size_t vertex2 = rng.get_size_t(min_index, max_index);
        if (vertex1 != vertex2) {
          neighbours[vertex1].emplace_back(vertex2);
        }
      }
    }

    // FINALLY: we've got components stuck together and singletons;
    // let's permute the vertex indices.
    // Element[i] tells the the NEW vertex index for the OLD index i.
    const auto new_indices = rng.get_permutation(next_vertex);

    data.raw_adjacency_data.resize(next_vertex);
    for (const auto& entry : neighbours) {
      auto& new_neighbours = data.raw_adjacency_data[new_indices[entry.first]];
      for (size_t old_index : entry.second) {
        new_neighbours.emplace_back(new_indices[old_index]);
      }
    }

    for (auto& old_component : data.components) {
      const auto old_index_set = old_component;
      old_component.clear();
      for (size_t old_index : old_index_set) {
        old_component.insert(new_indices[old_index]);
      }
    }
    return data;
  }
};

SCENARIO("Correctly calculates graph components") {
  RNG rng;
  ComponentsParameters components_parameters;

  for (components_parameters.max_tree_size = 3;
       components_parameters.max_tree_size < 20;
       components_parameters.max_tree_size += 5) {
    for (components_parameters.max_number_of_extra_edges_to_add_to_tree = 0;
         components_parameters.max_number_of_extra_edges_to_add_to_tree < 5;
         ++components_parameters.max_number_of_extra_edges_to_add_to_tree) {
      for (int counter = 0; counter < 5; ++counter) {
        const auto test_data = components_parameters.get_test_data(rng);
        const AdjacencyData cleaned_adjacency_data(
            test_data.raw_adjacency_data);

        const auto calculated_components =
            GraphRoutines::get_connected_components(cleaned_adjacency_data);

        // The components may come in a different order.
        // Since they're a partition, we can do a (lowest index) -> components
        // map.
        INFO(
            "counter="
            << counter
            << ", max_tree_size= " << components_parameters.max_tree_size
            << ", extra_edges added to tree = "
            << components_parameters.max_number_of_extra_edges_to_add_to_tree);
        REQUIRE(test_data.components.size() == calculated_components.size());

        // KEY: the smallest vertex index in a component.
        // VALUE: the index in the test_data vector of components.
        map<size_t, size_t> expected_components;

        for (size_t i = 0; i < test_data.components.size(); ++i) {
          const auto& single_component_set = test_data.components[i];
          REQUIRE(!single_component_set.empty());
          const size_t lowest_index = *single_component_set.crbegin();
          REQUIRE(expected_components.count(lowest_index) == 0);
          expected_components[lowest_index] = i;
        }

        // Now check the detailed components.
        for (const auto& calculated_component : calculated_components) {
          REQUIRE(!calculated_component.empty());
          const size_t lowest_index = *calculated_component.crbegin();
          REQUIRE(expected_components.count(lowest_index) != 0);
          const auto& expected_component_set =
              test_data.components[expected_components[lowest_index]];

          // The two sets must be identical.
          for (size_t exp_index : expected_component_set) {
            REQUIRE(calculated_component.count(exp_index) != 0);
          }
          for (size_t calc_index : calculated_component) {
            REQUIRE(expected_component_set.count(calc_index) != 0);
          }
        }
      }
    }
  }
}

}  // namespace tests
}  // namespace graphs
}  // namespace tket
