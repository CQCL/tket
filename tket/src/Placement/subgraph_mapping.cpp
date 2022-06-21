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

#include <algorithm>
#include <boost/algorithm/minmax_element.hpp>
#include <chrono>
#include <ctime>
#include <sstream>

#include "Architecture/Architecture.hpp"
#include "Graphs/Utils.hpp"
#include "Placement.hpp"
#include "Placement/Placement.hpp"
#include "Utils/Assert.hpp"
#include "Utils/GraphHeaders.hpp"
#include "Utils/TketLog.hpp"

namespace tket {

std::vector<qubit_bimap_t> monomorphism_edge_break(
    const Architecture& arc, const QubitGraph& q_graph, unsigned max_matches,
    unsigned timeout) {
  if (q_graph.n_nodes() > arc.n_nodes()) {
    throw ArchitectureInvalidity(
        "Interaction graph too large for architecture");
  }

  Architecture::UndirectedConnGraph undirected_target =
      arc.get_undirected_connectivity();
  QubitGraph::UndirectedConnGraph undirected_pattern =
      q_graph.get_undirected_connectivity();

  std::chrono::time_point<std::chrono::steady_clock> end_time =
      std::chrono::steady_clock::now() + std::chrono::milliseconds(timeout);

  while (true) {
    long search_timeout = std::chrono::duration_cast<std::chrono::milliseconds>(
                              (end_time - std::chrono::steady_clock::now()) / 2)
                              .count();
    if (search_timeout <= 0) search_timeout = 1;

    std::vector<qubit_bimap_t> all_maps = get_unweighted_subgraph_monomorphisms(
        undirected_pattern, undirected_target, max_matches, timeout);
    std::sort(all_maps.begin(), all_maps.end());

    if (std::chrono::steady_clock::now() >= end_time) {
      std::stringstream ss;
      ss << "subgraph monomorphism reached " << timeout
         << " millisecond timeout before reaching set max matches "
         << max_matches << ", instead finding " << all_maps.size()
         << " matches. "
            "Please change PlacementConfig.timeout to allow more matches.";
      tket_log()->warn(ss.str());
      if (all_maps.empty()) {
        throw std::runtime_error("No mappings found before timeout.");
      }
      return all_maps;
    }
    if (!all_maps.empty()) {
      return all_maps;
    }
    const unsigned current_number_of_edges =
        boost::num_edges(undirected_pattern);
    // It MUST have found a solution, if no pattern edges!
    TKET_ASSERT(current_number_of_edges > 0);
    using edge_t = graphs::utils::edge<QubitGraph::UndirectedConnGraph>;
    auto e_its = boost::edges(undirected_pattern);
    auto max_e_it = boost::first_max_element(
        e_its.first, e_its.second,
        [&undirected_pattern](const edge_t& lhs, const edge_t& rhs) {
          return undirected_pattern[lhs].weight <
                 undirected_pattern[rhs].weight;
        });
    graphs::utils::remove_edge(*max_e_it, undirected_pattern, true);
    TKET_ASSERT(boost::num_edges(undirected_pattern) < current_number_of_edges);
  }
}

}  // namespace tket
