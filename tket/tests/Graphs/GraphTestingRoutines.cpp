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

#include "GraphTestingRoutines.hpp"

#include <catch2/catch_test_macros.hpp>
#include <set>
#include <sstream>
#include <stdexcept>

#include "Graphs/AdjacencyData.hpp"
#include "Graphs/GraphColouring.hpp"

using std::exception;
using std::runtime_error;

using std::set;
using std::size_t;
using std::stringstream;
using std::vector;

namespace tket {
namespace graphs {
namespace tests {

void GraphTestingRoutines::require_valid_colouring_without_graph(
    const GraphColouringResult& colouring_result, bool require_no_colour_gaps) {
  const set<size_t> colours_seen(
      colouring_result.colours.cbegin(), colouring_result.colours.cend());

  try {
    if (require_no_colour_gaps && !colours_seen.empty()) {
      const auto min_col = *colours_seen.cbegin();
      const auto max_col = *colours_seen.crbegin();

      if (min_col != 0 || max_col + 1 != colours_seen.size()) {
        stringstream ss;
        ss << "The colours used should be an interval {0,1,2,...,m}"
              ", but we got min_col = "
           << min_col << ", max_col =" << max_col;
        throw runtime_error(ss.str());
      }
    }
    if (colours_seen.size() != colouring_result.number_of_colours) {
      throw runtime_error("number of colours mismatch");
    }
  } catch (const exception& e) {
    stringstream ss;
    ss << "total number of colours " << colours_seen.size() << ": "
       << colouring_result.to_string() << ": " << e.what();
    INFO("Invalid colouring object without graph: " << ss.str());
    REQUIRE(false);
  }
}

void GraphTestingRoutines::require_valid_colouring(
    const GraphColouringResult& colouring_result,
    const vector<std::pair<size_t, size_t>>& edges,
    size_t number_of_edges_to_check) {
  require_valid_colouring_without_graph(colouring_result);

  number_of_edges_to_check = std::min(number_of_edges_to_check, edges.size());

  for (size_t edge_index = 0; edge_index < number_of_edges_to_check;
       ++edge_index) {
    const auto& edge = edges[edge_index];
    try {
      if (edge.first == edge.second ||
          edge.first >= colouring_result.colours.size() ||
          edge.second >= colouring_result.colours.size()) {
        throw runtime_error("invalid vertex indices");
      }
      if (colouring_result.colours[edge.first] ==
          colouring_result.colours[edge.second]) {
        stringstream ss;
        ss << "adjacent vertices have the same colour "
           << colouring_result.colours[edge.first];
        throw runtime_error(ss.str());
      }
    } catch (const exception& e) {
      stringstream ss;
      ss << "checking colouring " << colouring_result.to_string()
         << " on edge index " << edge_index << " out of "
         << number_of_edges_to_check << " edges to check (" << edges.size()
         << " total). There is an edge between vertices " << edge.first << ", "
         << edge.second << ": " << e.what();
      INFO("GraphTestingRoutines: Invalid colouring for graph: " << ss.str());
      REQUIRE(false);
    }
  }
}

void GraphTestingRoutines::require_valid_suboptimal_colouring(
    const GraphColouringResult& colouring_result,
    const AdjacencyData& graph_data) {
  require_valid_colouring_without_graph(colouring_result);
  const auto number_of_vertices = colouring_result.colours.size();
  try {
    if (number_of_vertices != graph_data.get_number_of_vertices()) {
      throw runtime_error("Mismatching number of vertices");
    }
    for (size_t vertex = 0; vertex < number_of_vertices; ++vertex) {
      const auto& neighbours = graph_data.get_neighbours(vertex);
      for (size_t other_vertex : neighbours) {
        if (colouring_result.colours[vertex] ==
            colouring_result.colours[other_vertex]) {
          stringstream ss;
          ss << "adjacent vertices " << vertex << ", " << other_vertex
             << " have same colour " << colouring_result.colours[vertex];
          throw runtime_error(ss.str());
        }
      }
    }
  } catch (const exception& e) {
    stringstream ss;
    ss << "graph colouring: " << colouring_result.to_string() << "\nfor graph "
       << graph_data.to_string() << "\nwas invalid: " << e.what();
    INFO("GraphTestingRoutines: " << ss.str());
    REQUIRE(false);
  }
}

void GraphTestingRoutines::calculate_and_check_optimal_colouring(
    const GraphColouringResult& known_colouring,
    const AdjacencyData& adjacency_data) {
  require_valid_suboptimal_colouring(known_colouring, adjacency_data);

  const auto calculated_colouring =
      GraphColouringRoutines::get_colouring(adjacency_data);

  require_valid_suboptimal_colouring(calculated_colouring, adjacency_data);

  if (known_colouring.number_of_colours <
      calculated_colouring.number_of_colours) {
    INFO(
        "Graph: " << adjacency_data.to_string()
                  << "\nhas known valid colouring: "
                  << known_colouring.to_string()
                  << "\nwe recalculated the new colouring "
                  << calculated_colouring.to_string()
                  << "\nwhich cannot be optimal (it uses more colours)");

    REQUIRE(false);
  }
}

}  // namespace tests
}  // namespace graphs
}  // namespace tket
