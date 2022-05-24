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

#include "GraphColouring.hpp"

#include <algorithm>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

#include "AdjacencyData.hpp"
#include "BruteForceColouring.hpp"
#include "ColouringPriority.hpp"
#include "GraphRoutines.hpp"
#include "LargeCliquesResult.hpp"
#include "Utils/Assert.hpp"

using std::exception;
using std::map;
using std::runtime_error;
using std::set;
using std::string;
using std::stringstream;
using std::vector;

namespace tket {
namespace graphs {

string GraphColouringResult::to_string() const {
  stringstream ss;
  ss << "\nColouring: " << colours.size() << " vertices, " << number_of_colours
     << " colours : [ ";

  for (auto col : colours) {
    ss << col << ", ";
  }
  ss << "]";
  return ss.str();
}

GraphColouringResult::GraphColouringResult() {}

GraphColouringResult::GraphColouringResult(const vector<std::size_t>& _colours)
    : colours(_colours) {
  if (!colours.empty()) {
    number_of_colours = 1 + *std::max_element(colours.cbegin(), colours.cend());
  }
}

// Also updates the number of colours used in "result".
static void colour_single_component(
    const AdjacencyData& adjacency_data,
    const vector<set<std::size_t>>& connected_components,
    const vector<set<std::size_t>>& cliques, std::size_t component_index,
    GraphColouringResult& result) {
  result.number_of_colours =
      std::max(result.number_of_colours, cliques[component_index].size());

  const ColouringPriority colouring_priority(
      adjacency_data, connected_components[component_index],
      cliques[component_index]);

  const BruteForceColouring brute_force_colouring(
      colouring_priority, result.number_of_colours);

  const auto& partial_colour_map = brute_force_colouring.get_colours();

  for (const auto& entry : partial_colour_map) {
    const auto& vertex = entry.first;
    const auto& colour = entry.second;
    result.number_of_colours = std::max(result.number_of_colours, colour + 1);

    // GCOVR_EXCL_START
    try {
      if (vertex >= result.colours.size()) {
        throw runtime_error("illegal vertex index");
      }
      auto& colour_to_assign = result.colours[vertex];
      if (colour_to_assign < result.colours.size()) {
        stringstream ss;
        ss << "colour already assigned! Existing colour " << colour_to_assign;
        throw runtime_error(ss.str());
      }
      colour_to_assign = colour;
    } catch (const exception& e) {
      TKET_ASSERT(
          AssertMessage() << "colouring single component " << component_index
                          << " returned vertex " << vertex << " with colour "
                          << colour << " : " << e.what());
    }
    // GCOVR_EXCL_STOP
  }
}

// Check that everything was coloured,
// and we do have the correct number of colours
// (we might not have used all the colours we were allowed).
// (Don't bother trying to remove colour gaps, there shouldn't be any).
static void check_final_colouring(GraphColouringResult& result) {
  result.number_of_colours = 0;
  for (std::size_t i = 0; i < result.colours.size(); ++i) {
    const auto colour = result.colours[i];
    // GCOVR_EXCL_START
    if (colour >= result.colours.size()) {
      stringstream ss;
      ss << "vertex " << i << " has unassigned or illegal colour " << colour;
      throw runtime_error(ss.str());
    }
    // GCOVR_EXCL_STOP
    result.number_of_colours = std::max(result.number_of_colours, colour + 1);
  }
}

GraphColouringResult GraphColouringRoutines::get_colouring(
    const AdjacencyData& adjacency_data) {
  const auto connected_components =
      GraphRoutines::get_connected_components(adjacency_data);
  vector<set<std::size_t>> cliques(connected_components.size());
  vector<std::size_t> component_indices(connected_components.size());

  try {
    for (std::size_t i = 0; i < connected_components.size(); ++i) {
      const LargeCliquesResult cliques_in_this_component(
          adjacency_data, connected_components[i]);

      // GCOVR_EXCL_START
      if (cliques_in_this_component.cliques.empty()) {
        stringstream ss;
        ss << "component " << i << " has " << connected_components[i].size()
           << " vertices, but couldn't find a clique!";
        throw runtime_error(ss.str());
      }
      // GCOVR_EXCL_STOP
      cliques[i] = cliques_in_this_component.cliques[0];
      component_indices[i] = i;
    }

    // Now, we might as well start with the component containing the LARGEST
    // clique first (since, colouring becomes easier with more colours, and once
    // we've coloured one component, no point in trying to colour the others
    // with fewer colours).
    std::sort(
        component_indices.begin(), component_indices.end(),
        [&cliques](std::size_t lhs, std::size_t rhs) {
          return cliques[lhs].size() > cliques[rhs].size();
        });

    GraphColouringResult result;
    result.colours.assign(
        adjacency_data.get_number_of_vertices(),
        std::numeric_limits<std::size_t>::max());

    for (auto component_index : component_indices) {
      colour_single_component(
          adjacency_data, connected_components, cliques, component_index,
          result);
    }
    check_final_colouring(result);
    return result;
  } catch (const exception& e) {
    // GCOVR_EXCL_START
    TKET_ASSERT(
        AssertMessage() << "We had " << connected_components.size()
                        << " connected components, "
                        << adjacency_data.get_number_of_vertices()
                        << " vertices in total: " << e.what());
    // Some compilers error with "non-void function does not
    // return a value in all control paths..."
    return GraphColouringResult();
    // GCOVR_EXCL_STOP
  }
}

}  // namespace graphs
}  // namespace tket
