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
#include "Graphs/LargeCliquesResult.hpp"
#include "Utils/RNG.hpp"

using std::set;
using std::vector;

namespace tket {
namespace graphs {
namespace tests {

struct MaxCliqueTestData {
  vector<vector<size_t>> raw_adjacency_data;
  // It's not 100% guaranteed that there are no other cliques
  // of equal or larger size,
  // since by pure chance adding random edges might produce a larger clique.
  vector<set<size_t>> cliques;
};

// Used to test the maximum size clique function.
struct MaxCliqueParameters {
  size_t number_of_vertices = 20;
  size_t max_clique_size = 4;
  size_t approx_number_of_cliques = 3;
  size_t approx_number_of_extra_edges = 5;

  MaxCliqueTestData get_test_data(RNG& rng) const;
};

static void attempt_to_add_clique(
    const MaxCliqueParameters& parameters, RNG& rng, MaxCliqueTestData& data,
    set<size_t>& clique_vertices) {
  clique_vertices.clear();
  bool should_add_clique = false;
  for (size_t i = 0; i < 2 * parameters.max_clique_size; ++i) {
    clique_vertices.insert(rng.get_size_t(parameters.number_of_vertices - 1));
    if (clique_vertices.size() == parameters.max_clique_size) {
      should_add_clique = true;
      break;
    }
  }
  if (should_add_clique) {
    data.cliques.emplace_back(clique_vertices);
    // Add the actual edges to make this a clique.
    for (size_t i : clique_vertices) {
      for (size_t j : clique_vertices) {
        if (i != j) {
          data.raw_adjacency_data[i].push_back(j);
        }
      }
    }
  }
}

MaxCliqueTestData MaxCliqueParameters::get_test_data(RNG& rng) const {
  MaxCliqueTestData data;
  data.raw_adjacency_data.resize(number_of_vertices);
  set<size_t> clique_vertices;

  for (size_t counter = 0; counter < 2 * approx_number_of_cliques; ++counter) {
    attempt_to_add_clique(*this, rng, data, clique_vertices);
    if (data.cliques.size() >= approx_number_of_cliques) {
      break;
    }
  }
  // Finally, add extra edges.
  for (size_t counter = 0; counter < approx_number_of_extra_edges; ++counter) {
    const auto i = rng.get_size_t(number_of_vertices - 1);
    const auto j = rng.get_size_t(number_of_vertices - 1);
    if (i != j) {
      data.raw_adjacency_data[i].push_back(j);
    }
  }
  return data;
}

static bool is_clique(
    const set<size_t>& vertices, const AdjacencyData& cleaned_adjacency_data) {
  if (vertices.size() < 2) {
    return true;
  }
  // simplest: just count the edges
  size_t edge_count = 0;
  for (size_t i : vertices) {
    for (size_t j : vertices) {
      if (i < j && cleaned_adjacency_data.edge_exists(i, j)) {
        ++edge_count;
      }
    }
  }
  return edge_count * 2 == vertices.size() * (vertices.size() - 1);
}

static bool set_is_present(
    const set<size_t>& vertices, const vector<set<size_t>>& vertex_set_list) {
  // Note: cliques CAN overlap, e.g.
  // consider two triangles {1,2,3} and {1,2,4}, with no edge between 3,4.
  for (const auto& other_set : vertex_set_list) {
    size_t vertex_count = 0;
    for (size_t i : vertices) {
      vertex_count += other_set.count(i);
    }
    if (vertex_count == vertices.size() && vertex_count == other_set.size()) {
      return true;
    }
  }
  return false;
}

// Returns the size of each clique (they should all be equal)
static size_t check_that_calculated_cliques_are_valid(
    const vector<set<size_t>>& calculated_clique_data,
    const AdjacencyData& cleaned_adjacency_data) {
  set<size_t> clique_sizes;
  for (const auto& calc_clique : calculated_clique_data) {
    clique_sizes.insert(calc_clique.size());
    REQUIRE(is_clique(calc_clique, cleaned_adjacency_data));
  }
  // There ARE some cliques, and they are all of the same size!
  REQUIRE(clique_sizes.size() == 1);
  //...and of nonzero size...
  const size_t max_clique_size_in_this_component = *clique_sizes.cbegin();
  REQUIRE(max_clique_size_in_this_component > 0);
  return max_clique_size_in_this_component;
}

// Returns true if the calculated clique list does include the expected clique.
// Also return true if the expected clique size is strictly smaller
// than the calculated clique size,
// so we do not epect to see the clique (as it's not of MAXIMUM possible size).
static bool expected_clique_is_present(
    const set<size_t>& expected_clique_vertices,
    const vector<set<size_t>>& calculated_clique_data,
    const set<size_t>& component, size_t max_clique_size_in_this_component) {
  size_t expected_vertices_present = 0;
  for (size_t vertex : expected_clique_vertices) {
    expected_vertices_present += component.count(vertex);
  }
  if (expected_vertices_present == 0) {
    return false;
  }
  REQUIRE(expected_vertices_present == expected_clique_vertices.size());

  // Now, we have an EXPECTED clique lying entirely within this component:
  // so EITHER it equals one of the calc cliques,
  // OR the calc clique is bigger.
  if (expected_vertices_present >= max_clique_size_in_this_component) {
    // Otherwise, purely by chance, i.e. adding extra edges, we have created
    // a strictly larger clique, so we don't expect to see our clique.
    REQUIRE(expected_vertices_present == max_clique_size_in_this_component);
    REQUIRE(set_is_present(expected_clique_vertices, calculated_clique_data));
  }
  return true;
}

// Just within a single component, check and compare
// the expected/calculated cliques.
static void test_cliques_in_single_component(
    const MaxCliqueTestData& test_data,
    const AdjacencyData& cleaned_adjacency_data, const set<size_t>& component,
    set<size_t>& clique_indices_seen) {
  const LargeCliquesResult calculated_clique_result(
      cleaned_adjacency_data, component, 1000);

  REQUIRE(calculated_clique_result.cliques_are_definitely_max_size);

  const size_t max_calc_clique_size_in_this_component =
      check_that_calculated_cliques_are_valid(
          calculated_clique_result.cliques, cleaned_adjacency_data);

  // Which EXPECTED cliques in this component are present
  // in the list of calculated cliques?
  for (size_t clique_index = 0; clique_index < test_data.cliques.size();
       ++clique_index) {
    const auto& expected_clique_vertices = test_data.cliques[clique_index];
    if (expected_clique_is_present(
            expected_clique_vertices, calculated_clique_result.cliques,
            component, max_calc_clique_size_in_this_component)) {
      clique_indices_seen.insert(clique_index);
    }
  }
}

// Returns the number of cliques seen, just as an extra check.
static size_t test_max_clique_generated_data(
    const MaxCliqueTestData& test_data) {
  // We'll check at the end that every expected clique DID occur.
  // The indices are for the vector of index sets.
  set<size_t> clique_indices_seen;

  const AdjacencyData cleaned_adjacency_data(test_data.raw_adjacency_data);

  const auto components =
      GraphRoutines::get_connected_components(cleaned_adjacency_data);

  for (const auto& component : components) {
    test_cliques_in_single_component(
        test_data, cleaned_adjacency_data, component, clique_indices_seen);
  }
  // The clique indices seen must be contiguous and complete,
  // i.e. all indices 0,1,...,N-1.
  // Otherwise, we're missing a clique.
  REQUIRE(clique_indices_seen.size() == test_data.cliques.size());
  if (!clique_indices_seen.empty()) {
    REQUIRE(*clique_indices_seen.crbegin() + 1 == test_data.cliques.size());
  }
  return clique_indices_seen.size();
}

SCENARIO("Correctly calculates max cliques") {
  RNG rng;
  MaxCliqueParameters parameters;
  size_t cliques_seen = 0;

  for (parameters.number_of_vertices = 10; parameters.number_of_vertices < 50;
       parameters.number_of_vertices += 20) {
    for (parameters.max_clique_size = 2; parameters.max_clique_size <= 5;
         ++parameters.max_clique_size) {
      for (parameters.approx_number_of_cliques = 1;
           parameters.approx_number_of_cliques < 5;
           parameters.approx_number_of_cliques += 2) {
        for (int counter = 0; counter < 5; ++counter) {
          INFO(
              "number_of_vertices = "
              << parameters.number_of_vertices << ", max_clique_size = "
              << parameters.max_clique_size << ", approx_number_of_cliques = "
              << parameters.approx_number_of_cliques
              << ", test counter = " << counter);

          const auto test_data = parameters.get_test_data(rng);
          cliques_seen += test_max_clique_generated_data(test_data);
        }
      }
    }
  }
  CHECK(cliques_seen == 160);
}

}  // namespace tests
}  // namespace graphs
}  // namespace tket
