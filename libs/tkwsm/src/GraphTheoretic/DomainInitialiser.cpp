// Copyright 2019-2023 Cambridge Quantum Computing
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

#include "tkwsm/GraphTheoretic/DomainInitialiser.hpp"

#include <algorithm>
#include <tkassert/Assert.hpp>

#include "tkwsm/Common/GeneralUtils.hpp"
#include "tkwsm/GraphTheoretic/FilterUtils.hpp"
#include "tkwsm/GraphTheoretic/NearNeighboursData.hpp"
#include "tkwsm/GraphTheoretic/NeighboursData.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

bool DomainInitialiser::full_initialisation(
    InitialDomains& initial_domains,
    const NeighboursData& pattern_neighbours_data,
    NearNeighboursData& pattern_near_neighbours_data,
    const NeighboursData& target_neighbours_data,
    NearNeighboursData& target_near_neighbours_data, unsigned max_path_length) {
  return degree_sequence_initialisation(
             initial_domains, pattern_neighbours_data,
             target_neighbours_data) &&

         distance_counts_reduction(
             initial_domains, pattern_near_neighbours_data,
             target_near_neighbours_data, max_path_length);
}

bool DomainInitialiser::degree_sequence_initialisation(
    InitialDomains& initial_domains,
    const NeighboursData& pattern_neighbours_data,
    const NeighboursData& target_neighbours_data) {
  TKET_ASSERT(initial_domains.empty());
  initial_domains.resize(
      pattern_neighbours_data.get_number_of_nonisolated_vertices());

  const auto number_of_tv =
      target_neighbours_data.get_number_of_nonisolated_vertices();

  // Now get the target degree sequences.
  // For each PV, its domain will be all those TV
  // with a compatible degree sequence.
  // Thus we need to calculate data for ALL TV right now,
  // but the PV degree sequences can be discarded after use.
  //
  // NOTE: if in future we RESTRICT the initial domains before reaching
  // this point, then maybe we should calculate TV data lazily,
  // in case some vertices are never used.
  //
  // NOTE: it could be more dynamic, i.e. if a TV is not
  // in ANY domain, then we could erase it from the target graph completely.
  // We'd then recursively recalculate and break down,
  // and maybe end up reducing the target graph quite a bit in some cases.
  // However, a lot of extra calculation and computation,
  // probably not worth it
  // (ESPECIALLY as, in many applications, i.e. quantum computer
  // target architectures, the target graph has a lot of
  // homogeneity: e.g., IBM heavy hexagons, square grids, etc.;
  // so it's unlikely that many target vertices can be removed in this way).

  // FIRST: the target vertex;
  // SECOND: its degree sequence.
  typedef std::pair<VertexWSM, std::vector<std::size_t>> DegSeqData;
  std::vector<DegSeqData> target_degree_sequences;

  target_degree_sequences.reserve(number_of_tv);
  for (unsigned tv = 0; tv < number_of_tv; ++tv) {
    target_degree_sequences.emplace_back(
        tv, target_neighbours_data.get_sorted_degree_sequence_expensive(tv));
  }

  // Sort by decreasing sequence length.
  std::sort(
      target_degree_sequences.begin(), target_degree_sequences.end(),
      [](const DegSeqData& lhs, const DegSeqData& rhs) -> bool {
        return lhs.second.size() > rhs.second.size() ||
               // This is not crucial, it just ensures exact behaviour for
               // nonstable sorts:
               (lhs.second.size() == rhs.second.size() &&
                lhs.first < rhs.first);
      });

  std::vector<std::size_t> pattern_sequence;
  for (unsigned pv = 0; pv < initial_domains.size(); ++pv) {
    initial_domains[pv].resize(number_of_tv);

    pattern_sequence =
        pattern_neighbours_data.get_sorted_degree_sequence_expensive(pv);

    // Now, which target vertices have compatible degree sequences?
    for (const DegSeqData& deg_seq_data : target_degree_sequences) {
      const auto& target_sequence = deg_seq_data.second;
      if (target_sequence.size() < pattern_sequence.size()) {
        break;
      }
      if (FilterUtils::compatible_sorted_degree_sequences(
              pattern_sequence, target_sequence)) {
        initial_domains[pv].set(deg_seq_data.first);
      }
    }
    if (initial_domains[pv].none()) {
      return false;
    }
  }
  return true;
}

bool DomainInitialiser::distance_counts_reduction(
    InitialDomains& initial_domains,
    NearNeighboursData& pattern_near_neighbours_data,
    NearNeighboursData& target_near_neighbours_data, unsigned max_path_length) {
  if (max_path_length <= 1) {
    // Neighbour counts are already included in degree sequences.
    // Don't waste time recalculating.
    return true;
  }
  // Easier to build the TV list up, then erase them in a separate pass.
  boost::dynamic_bitset<> tv_to_erase(
      target_near_neighbours_data.get_number_of_vertices());

  for (unsigned pv = 0; pv < initial_domains.size(); ++pv) {
    boost::dynamic_bitset<>& domain = initial_domains[pv];
    tv_to_erase.reset();

    // Find which TV to erase.
    for (auto tv = domain.find_first(); tv < domain.size();
         tv = domain.find_next(tv)) {
      // If f is a valid monomorphism then
      //    dist(PV, PV')=d   ==>   dist(f(PV), f(PV')) <= d.
      // Now TV=f(PV) is the value to be checked;
      // so for each d, just count how many PV' are at that exact distance d,
      // and check that there are at least that many TV' with
      // dist(TV, TV') <= d.
      //
      // Notice, that we may break off early; because everything
      // is calculated lazily, it may save time.
      for (unsigned distance = 2; distance <= max_path_length; ++distance) {
        const auto number_of_pv_at_distance_d =
            pattern_near_neighbours_data.get_n_vertices_at_exact_distance(
                pv, distance);

        if (number_of_pv_at_distance_d > 0) {
          const auto number_of_tv_at_distance_le_d =
              target_near_neighbours_data.get_n_vertices_up_to_distance(
                  tv, distance);
          if (number_of_pv_at_distance_d > number_of_tv_at_distance_le_d) {
            tv_to_erase.set(tv);
            break;
          }
        }
      }
    }
    domain -= tv_to_erase;
    if (domain.none()) {
      return false;
    }
  }
  return true;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
