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

#include "EdgeSequenceColouringParameters.hpp"

#include <catch2/catch_test_macros.hpp>

#include "EdgeSequence.hpp"
#include "GraphTestingRoutines.hpp"
#include "Graphs/AdjacencyData.hpp"
#include "Graphs/GraphColouring.hpp"
#include "RandomGraphGeneration.hpp"

using std::size_t;
using std::string;
using std::vector;

namespace tket {
namespace graphs {
namespace tests {

size_t EdgeSequenceColouringParameters::test_colourings(
    RandomGraphParameters& parameters, EdgeSequence& edge_sequence) {
  const auto number_of_vertices =
      edge_sequence.adjacency_data.get_number_of_vertices();

  if (step_skip_size == 0 || max_number_of_colourings == 0 ||
      max_number_of_edges == 0 || number_of_vertices < 3) {
    return 0;
  }

  edge_sequence.edges.clear();
  parameters.add_edges(edge_sequence);

  size_t number_of_colourings = 0;
  size_t previous_number_of_colours = 0;

  const size_t known_max_chromatic_number =
      std::min(parameters.max_chromatic_number, number_of_vertices);

  // TODO: if the tests start getting slow, consider
  // not doing all this string manipulation until an error actually occurs.

  string colouring_str = "\nInitial known colouring: ";
  if (parameters.known_colouring.colours.empty()) {
    colouring_str += "[none]";
  } else {
    colouring_str += parameters.known_colouring.to_string();
  }
  AdjacencyData growing_adjacency_data(number_of_vertices);
  string prev_graph;
  string prev_colouring;

  for (size_t edge_index = 0; edge_index < max_number_of_edges; ++edge_index) {
    if (edge_index >= edge_sequence.edges.size() ||
        number_of_colourings >= max_number_of_colourings) {
      return number_of_colourings;
    }
    const auto vertex_pair = edge_sequence.edges[edge_index];
    growing_adjacency_data.add_edge(vertex_pair.first, vertex_pair.second);

    if (!parameters.known_colouring.colours.empty()) {
      GraphTestingRoutines::require_valid_suboptimal_colouring(
          parameters.known_colouring, growing_adjacency_data);
    }
    if (edge_index % step_skip_size != 0) {
      // No colouring here.
      continue;
    }
    // We are going to colour!
    ++number_of_colourings;

    const auto calculated_colouring =
        GraphColouringRoutines::get_colouring(growing_adjacency_data);

    ++total_number_of_colourings;

    INFO(
        "current edge index="
        << edge_index << ", calculated colouring: "
        << calculated_colouring.to_string() << colouring_str
        << "\nfor random graph: " << growing_adjacency_data.to_string()
        << "\n\nPREV graph was: " << prev_graph << "\nwith colouring "
        << prev_colouring);

    REQUIRE(
        calculated_colouring.number_of_colours <= known_max_chromatic_number);

    REQUIRE(
        calculated_colouring.number_of_colours >= previous_number_of_colours);

    if (!parameters.known_colouring.colours.empty()) {
      // We have a known colouring valid for ALL graphs in this sequence.
      REQUIRE(
          calculated_colouring.number_of_colours <=
          parameters.known_colouring.number_of_colours);
    }
    previous_number_of_colours = calculated_colouring.number_of_colours;

    // Check that the colourings are valid for ALL previous graphs,
    // as well as the current one.
    GraphTestingRoutines::require_valid_colouring(
        calculated_colouring, edge_sequence.edges, edge_index + 1);

    prev_graph = growing_adjacency_data.to_string();
    prev_colouring = calculated_colouring.to_string();
  }
  return number_of_colourings;
}

}  // namespace tests
}  // namespace graphs
}  // namespace tket
