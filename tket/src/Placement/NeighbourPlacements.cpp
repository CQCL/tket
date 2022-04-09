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

#include "NeighbourPlacements.hpp"

#include <cstdlib>
#include <sstream>
#include <string>

#include "TokenSwapping/SwapListOptimiser.hpp"
#include "Utils/TketLog.hpp"

namespace tket {

NeighbourPlacements::NeighbourPlacements(
    const Architecture& arc, const qubit_mapping_t& init_map)
    : arc_(arc), init_map_(init_map), u_to_node_(), rng_() {
  auto nodes = arc_.get_all_nodes_vec();
  for (unsigned i = 0; i < nodes.size(); ++i) {
    u_to_node_.left.insert({i, nodes[i]});
  }
}

NeighbourPlacements::ResultVec NeighbourPlacements::get(
    unsigned dist, unsigned n, bool optimise, unsigned seed,
    unsigned max_tries) {
  rng_.set_seed(seed);

  // define a comparison function for placements
  std::vector<Qubit> keys;
  for (auto [k, v] : init_map_) {
    keys.push_back(k);
  }
  auto map_compare = [&keys](
                         const qubit_mapping_t& a, const qubit_mapping_t& b) {
    for (auto k : keys) {
      if (a.at(k) < b.at(k)) {
        return true;
      } else if (a.at(k) > b.at(k)) {
        return false;
      }
    }
    return false;
  };
  // set of all generated placement maps
  std::set<qubit_mapping_t, decltype(map_compare)> placements(map_compare);

  ResultVec resvec;
  for (unsigned i = 0; i < n; ++i) {
    unsigned n_unsuccessful = 0;
    while (n_unsuccessful < max_tries) {
      Result res = gen_result(dist, optimise, max_tries);
      if (!placements.contains(res.map)) {
        resvec.push_back(res);
        placements.insert(res.map);
        break;
      }
      ++n_unsuccessful;
    }
    if (n_unsuccessful == max_tries) {
      std::stringstream ss;
      ss << "Could not generate " << n << " distinct placements";
      tket_log()->warn(ss.str());
    }
  }
  return resvec;
}

NeighbourPlacements::Result NeighbourPlacements::gen_result(
    unsigned dist, bool optimise, unsigned max_tries) {
  SwapList swaps;
  tsa_internal::SwapListOptimiser optimiser;

  // it might be impossible to find `dist` non-trivial swaps
  unsigned n_unsuccessful = 0;

  while (swaps.size() < dist && n_unsuccessful < max_tries) {
    Swap new_swap = gen_swap();

    if (optimise) {
      SwapList swaps_candidate = swaps;
      swaps_candidate.push_back(new_swap);
      optimiser.full_optimise(swaps_candidate);
      if (swaps_candidate.size() > swaps.size()) {
        swaps = std::move(swaps_candidate);
        n_unsuccessful = 0;
      } else {
        ++n_unsuccessful;
      }
    } else {
      swaps.push_back(new_swap);
    }
  }

  if (n_unsuccessful == max_tries) {
    std::stringstream ss;
    ss << "Unable to generate " << dist << " swaps for given architecture";
    tket_log()->warn(ss.str());
  }

  return convert_to_res(swaps.to_vector());
}

Swap NeighbourPlacements::gen_swap() {
  auto edges = arc_.get_all_edges_vec();
  unsigned m = edges.size();
  auto [n1, n2] = edges[rng_.get_size_t(m - 1)];
  Swap new_swap{u_to_node_.right.at(n1), u_to_node_.right.at(n2)};
  return new_swap;
}

NeighbourPlacements::Result NeighbourPlacements::convert_to_res(
    const SwapVec& swaps) {
  NodeSwapVec node_swaps;
  for (auto [u1, u2] : swaps) {
    node_swaps.push_back({u_to_node_.left.at(u1), u_to_node_.left.at(u2)});
  }

  qubit_bimap_t qubit_to_node;
  qubit_to_node.left.insert(init_map_.begin(), init_map_.end());
  for (auto [n1, n2] : node_swaps) {
    const Qubit q1 = qubit_to_node.right.at(n1);
    const Qubit q2 = qubit_to_node.right.at(n2);
    qubit_to_node.left.erase(q1);
    qubit_to_node.left.erase(q2);
    qubit_to_node.left.insert({q1, n2});
    qubit_to_node.left.insert({q2, n1});
  }
  qubit_mapping_t map;
  for (auto [k, v] : qubit_to_node.left) {
    map.insert({k, v});
  }
  return {map, node_swaps};
}

}  // namespace tket
