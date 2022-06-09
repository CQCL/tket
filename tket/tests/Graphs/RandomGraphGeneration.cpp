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

#include "RandomGraphGeneration.hpp"

#include <algorithm>
#include <catch2/catch_test_macros.hpp>

#include "EdgeSequence.hpp"
#include "Graphs/AdjacencyData.hpp"
#include "Utils/RNG.hpp"

using std::vector;

namespace tket {
namespace graphs {
namespace tests {

bool RandomFibrousGraphParameters::add_edges(EdgeSequence& edge_sequence) {
  const size_t number_of_vertices =
      edge_sequence.adjacency_data.get_number_of_vertices();

  REQUIRE(number_of_vertices >= 2);
  REQUIRE(number_of_vertices <= 10000);

  const size_t max_vertex = number_of_vertices - 1;
  auto& rng = edge_sequence.rng;

  for (size_t i = 0; i < number_of_strands; ++i) {
    const size_t start_v = rng.get_size_t(max_vertex);
    size_t current_v = start_v;

    for (int counter = 0; counter < 1000; ++counter) {
      auto new_v = current_v;
      bool should_break = false;

      if (rng.check_percentage(percentage_for_each_strand_to_grow)) {
        new_v = rng.get_size_t(max_vertex);
      } else {
        should_break = true;
        // The strand ends now; but should we close it up to make a cycle?
        if (rng.check_percentage(percentage_for_strand_to_become_a_cycle)) {
          new_v = start_v;
        }
      }
      if (current_v != new_v) {
        edge_sequence.add_edge(current_v, new_v);
        break;
      }
      if (should_break) {
        break;
      }
    }
  }
  // Nothing to stop you adding more edges later.
  return true;
}

bool RandomTreeParameters::add_edges(EdgeSequence& edge_sequence) {
  // Easy to prove that all trees can be 2-coloured.
  max_chromatic_number = 2;

  const size_t number_of_vertices =
      edge_sequence.adjacency_data.get_number_of_vertices();

  REQUIRE(number_of_vertices >= 2);
  REQUIRE(number_of_vertices <= 10000);

  const size_t max_vertex = number_of_vertices - 1;
  auto& rng = edge_sequence.rng;
  const auto root_vertex = rng.get_size_t(max_vertex);

  // TODO: be fancy and use a single partitioned vector
  // for existing/unused vertices.
  vector<size_t> unused_vertices;
  unused_vertices.reserve(max_vertex);

  for (size_t ii = 0; ii < number_of_vertices; ++ii) {
    if (ii != root_vertex) {
      unused_vertices.push_back(ii);
    }
  }
  vector<size_t> existing_vertices;
  existing_vertices.push_back(root_vertex);

  for (size_t counter1 = 0; counter1 < approx_number_of_spawns; ++counter1) {
    const auto parent = rng.get_element(existing_vertices);

    for (size_t counter2 = 0; counter2 < approx_number_of_children_per_node;
         ++counter2) {
      if (unused_vertices.empty()) {
        return false;
      }

      const auto child_vertex_index =
          rng.get_size_t(unused_vertices.size() - 1);
      auto& vertex_entry = unused_vertices[child_vertex_index];

      edge_sequence.add_edge(parent, vertex_entry);

      existing_vertices.push_back(vertex_entry);

      // Remove the child vertex; simply swap with the last one
      vertex_entry = unused_vertices.back();
      unused_vertices.pop_back();
    }
  }
  return false;
}

bool RandomDenseGraphParameters::add_edges(EdgeSequence& edge_sequence) {
  const size_t number_of_vertices =
      edge_sequence.adjacency_data.get_number_of_vertices();

  REQUIRE(number_of_vertices >= 2);
  REQUIRE(number_of_vertices <= 10000);

  const size_t max_vertex = number_of_vertices - 1;
  auto& rng = edge_sequence.rng;

  size_t approx_max_total_attempts = 100000;
  approx_max_total_attempts = std::min(
      approx_max_total_attempts, number_of_vertices * number_of_vertices);

  size_t max_attempts = 10000;
  max_attempts =
      std::min(max_attempts, max_number_of_consecutive_add_edge_attempts);

  for (; approx_max_total_attempts > 0; --approx_max_total_attempts) {
    bool edge_added = false;
    for (size_t counter = 0; counter < max_attempts; ++counter) {
      const auto ii = rng.get_size_t(max_vertex);
      const auto jj = rng.get_size_t(max_vertex);
      if (ii != jj && edge_sequence.add_edge(ii, jj)) {
        edge_added = true;
        break;
      }
    }
    if (!edge_added) {
      return false;
    }
  }
  return false;
}

namespace {

// For initially randomly colouring the vertices.
// Once the vertex colours are known, it's trivial to add edges.
struct VerticesPartition {
  // Each inner vector element[i] lists all vertices with colour i.
  vector<vector<size_t>> single_colour_vertex_sets;

  vector<size_t> known_colouring;

  VerticesPartition(
      EdgeSequence& edge_sequence, size_t max_number_of_colours_to_use) {
    size_t max_colour_used = 0;

    // In each pair: first is the vertex, second is the colour.
    vector<std::pair<size_t, size_t>> colours(
        edge_sequence.adjacency_data.get_number_of_vertices());
    for (size_t ii = 0; ii < colours.size(); ++ii) {
      const size_t colour = ii % max_number_of_colours_to_use;
      max_colour_used = std::max(max_colour_used, colour);
      colours[ii].first = ii;
      colours[ii].second = colour;
    }
    edge_sequence.rng.do_shuffle(colours);
    single_colour_vertex_sets.resize(max_colour_used + 1);
    known_colouring.resize(colours.size());

    for (const auto& entry : colours) {
      const auto& vertex = entry.first;
      const auto& col = entry.second;
      known_colouring[vertex] = col;
      single_colour_vertex_sets[col].push_back(vertex);
    }
  }
};

}  // unnamed namespace

bool RandomColouredDenseGraphParameters::add_edges(
    EdgeSequence& edge_sequence) {
  const VerticesPartition partition(
      edge_sequence, max_number_of_colours_to_use);

  max_chromatic_number = partition.single_colour_vertex_sets.size();
  const size_t max_colour = max_chromatic_number - 1;

  known_colouring = GraphColouringResult(partition.known_colouring);

  // If there are N vertices, there are O(N^2) possible edges,
  // so if you still haven't found one after ~N^2 steps, time to give up.
  size_t max_attempts = edge_sequence.adjacency_data.get_number_of_vertices();
  max_attempts *= max_attempts;

  max_attempts =
      std::min(max_attempts, max_number_of_consecutive_add_edge_attempts);

  auto& rng = edge_sequence.rng;

  for (;;) {
    bool edge_added = false;
    for (size_t counter = 0; counter < max_attempts; ++counter) {
      const auto colour1 = rng.get_size_t(max_colour);
      const auto colour2 = rng.get_size_t(max_colour);

      if (colour1 == colour2) {
        continue;
      }
      const auto& vertex1 =
          rng.get_element(partition.single_colour_vertex_sets[colour1]);

      const auto& vertex2 =
          rng.get_element(partition.single_colour_vertex_sets[colour2]);

      // Guaranteed to be different vertices of different colours;
      // but was the edge seen before?
      edge_added = edge_sequence.add_edge(vertex1, vertex2);
      if (edge_added) {
        break;
      }
    }
    if (!edge_added) {
      return false;
    }
  }
  return false;
}

// Note: every graph that can be generated
// by this can, I believe, also be generated by
// RandomColouredDenseGraphParameters, by chance.
// However, the algorithms and probability distributions
// are so different that, in practice, this might as well be regarded
// as a totally different class of random graphs, so it's worth doing.
// PHILOSOPHY: we generate random test data in the hope of catching
// our subtle mistakes, due to a lack of understanding.
// Almost no human brain could fully understand, say, the set of
// all graphs with 10 vertices. We shouldn't be shy about
// running many different tests, even if (in theory) they are duplicates,
// but (in practice) they are not.
bool RandomColouredKPartiteGraphParameters::add_edges(
    EdgeSequence& edge_sequence) {
  const auto number_of_vertices =
      number_of_vertex_sets * number_of_vertices_in_each_set;

  edge_sequence.clear();
  edge_sequence.adjacency_data.clear(number_of_vertices);

  // Of course this MIGHT be suboptimal, if enough edges are missing.
  known_colouring.number_of_colours = number_of_vertex_sets;
  known_colouring.colours.resize(number_of_vertices);

  if (number_of_vertices == 0) {
    return false;
  }
  auto& rng = edge_sequence.rng;
  const auto labels = rng.get_permutation(number_of_vertices);

  // k is the number of colours.
  const auto& k = number_of_vertex_sets;
  const auto& m = number_of_vertices_in_each_set;

  for (size_t c0 = 0; c0 < k; ++c0) {
    for (size_t c1 = c0 + 1; c1 < k; ++c1) {
      for (size_t r0 = 0; r0 < m; ++r0) {
        for (size_t r1 = 0; r1 < m; ++r1) {
          size_t n0 = m * c0 + r0;
          size_t n1 = m * c1 + r1;
          auto l0 = labels[n0];
          auto l1 = labels[n1];

          if (rng.check_percentage(percentage_of_added_edges)) {
            INFO(
                "n0=" << n0 << ", n1=" << n1 << ", l0=" << l0 << ", l1=" << l1);
            REQUIRE(edge_sequence.add_edge(l0, l1));
          }
        }
      }
    }
  }
  for (size_t c = 0; c < k; ++c) {
    for (size_t r = 0; r < m; ++r) {
      known_colouring.colours[labels[m * c + r]] = c;
    }
  }
  return false;
}

bool EdgelessGraph::add_edges(EdgeSequence& edge_sequence) {
  max_chromatic_number = 1;
  known_colouring.number_of_colours = 1;
  known_colouring.colours.assign(
      edge_sequence.adjacency_data.get_number_of_vertices(), 0);

  return false;
}

bool CompleteGraph::add_edges(EdgeSequence& edge_sequence) {
  max_chromatic_number = edge_sequence.adjacency_data.get_number_of_vertices();

  known_colouring.number_of_colours = max_chromatic_number;
  known_colouring.colours.resize(max_chromatic_number);

  for (size_t ii = 0; ii < max_chromatic_number; ++ii) {
    known_colouring.colours[ii] = ii;

    for (size_t jj = ii + 1; jj < max_chromatic_number; ++jj) {
      REQUIRE(edge_sequence.add_edge(ii, jj));
    }
  }
  REQUIRE(
      edge_sequence.edges.size() ==
      (max_chromatic_number * (max_chromatic_number - 1)) / 2);
  return false;
}

}  // namespace tests
}  // namespace graphs
}  // namespace tket
