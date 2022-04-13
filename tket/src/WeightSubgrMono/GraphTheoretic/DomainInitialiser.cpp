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

#include "WeightSubgrMono/GraphTheoretic/DomainInitialiser.hpp"

#include <algorithm>

#include "WeightSubgrMono/Common/GeneralUtils.hpp"
#include "WeightSubgrMono/GraphTheoretic/FilterUtils.hpp"
#include "WeightSubgrMono/GraphTheoretic/NearNeighboursData.hpp"
#include "WeightSubgrMono/GraphTheoretic/NeighboursData.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

bool DomainInitialiser::full_initialisation(
    PossibleAssignments& possible_assignments,
    const NeighboursData& pattern_neighbours_data,
    NearNeighboursData& pattern_near_neighbours_data,
    const NeighboursData& target_neighbours_data,
    NearNeighboursData& target_near_neighbours_data, unsigned max_path_length) {
  return degree_sequence_initialisation(
             possible_assignments, pattern_neighbours_data,
             target_neighbours_data) &&

         distance_counts_reduction(
             possible_assignments, pattern_neighbours_data,
             pattern_near_neighbours_data, target_neighbours_data,
             target_near_neighbours_data, max_path_length);
}

bool DomainInitialiser::degree_sequence_initialisation(
    PossibleAssignments& possible_assignments,
    const NeighboursData& pattern_neighbours_data,
    const NeighboursData& target_neighbours_data) {
  possible_assignments.clear();
  for (const auto& entry : pattern_neighbours_data.get_map()) {
    possible_assignments[entry.first];
  }

  // Get the target degree sequences.

  // FIRST: the target vertex;
  // SECOND: its degree sequence.
  typedef std::pair<VertexWSM, std::vector<std::size_t>> DegSeqData;
  std::vector<DegSeqData> target_degree_sequences;

  const auto& target_map = target_neighbours_data.get_map();
  target_degree_sequences.reserve(target_map.size());
  for (const auto& t_entry : target_map) {
    const VertexWSM tv = t_entry.first;
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
  for (auto& entry : possible_assignments) {
    const auto& pv = entry.first;
    auto& domain = entry.second;
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
        domain.insert(deg_seq_data.first);
      }
    }
    if (domain.empty()) {
      return false;
    }
  }
  return true;
}

namespace {
// We may not need all TVs, or all distances, so just fill in the map lazily.
class TCountsLazyCalculator {
 public:
  TCountsLazyCalculator(
      const NeighboursData& target_neighbours_data,
      NearNeighboursData& target_near_neighbours_data)
      : m_target_neighbours_data(target_neighbours_data),
        m_target_near_neighbours_data(target_near_neighbours_data) {}

  const std::vector<std::size_t>& operator()(
      VertexWSM tv, unsigned max_distance) {
    const auto citer = m_target_counts_map.find(tv);
    if (citer != m_target_counts_map.cend() &&
        (citer->second.size() >= max_distance ||
         // Once it ends with a zero, no point in extending
         (!citer->second.empty() && citer->second.back() == 0))) {
      return citer->second;
    }
    auto& t_counts = m_target_counts_map[tv];
    m_target_near_neighbours_data.fill_counts_vector(
        tv, max_distance, t_counts);
    return t_counts;
  }

 private:
  const NeighboursData& m_target_neighbours_data;
  NearNeighboursData& m_target_near_neighbours_data;
  std::map<VertexWSM, std::vector<std::size_t>> m_target_counts_map;
};
}  // namespace

static void get_tv_to_erase(
    unsigned p_distance, std::vector<VertexWSM>& tv_to_erase,
    const std::vector<std::size_t>& pattern_counts,
    const std::set<VertexWSM>& domain,
    TCountsLazyCalculator& t_counts_calculator) {
  tv_to_erase.clear();
  for (auto tv : domain) {
    bool tv_is_acceptable = false;
    // Grow the t-subgraph. If the embedding succeeds for some d
    // then we can move onto another TV to test.
    for (unsigned t_distance = 2; t_distance <= p_distance; ++t_distance) {
      if (NearNeighboursData::test_against_target(
              pattern_counts, t_counts_calculator(tv, t_distance))) {
        // It's pointless to increase the t-distance once it has passed
        // for some d; the target subgraph is big enough.
        tv_is_acceptable = true;
        break;
      }
    }
    if (!tv_is_acceptable) {
      // TV failed all the way up to distance d, so this p-subgraph
      // definitely cannot embed with TV as the root.
      tv_to_erase.push_back(tv);
      break;
    }
  }
}

bool DomainInitialiser::distance_counts_reduction(
    PossibleAssignments& possible_assignments,
    const NeighboursData& pattern_neighbours_data,
    NearNeighboursData& pattern_near_neighbours_data,
    const NeighboursData& target_neighbours_data,
    NearNeighboursData& target_near_neighbours_data, unsigned max_path_length) {
  if (max_path_length <= 1) {
    // Neighbour counts are already included in degree sequences.
    return true;
  }

  TCountsLazyCalculator t_counts_calculator(
      target_neighbours_data, target_near_neighbours_data);

  std::vector<VertexWSM> tv_to_erase;

  // We only consider PV one-by-one, so we don't need all PV data at once.
  std::vector<std::size_t> pattern_counts;

  // What's the best approach, to terminate quickly upon an impossible problem?
  // We assume that growing a subgraph
  // is slower than testing/copying the distance vectors.
  // Thus we gradually grow the p-subgraphs "in parallel",
  // and we'll know as soon as one fails to embed in any TV
  // that the full problem is insoluble.

  // KEY: pv  VALUE: the previous length of the counts list;
  // stop as soon as it stops growing.
  std::map<VertexWSM, unsigned> last_p_subgraph_count_length;

  std::set<VertexWSM> assigned_pv;

  for (unsigned p_distance = 2; p_distance <= max_path_length; ++p_distance) {
    for (auto& entry : possible_assignments) {
      const VertexWSM& pv = entry.first;
      auto& domain = entry.second;

      pattern_near_neighbours_data.fill_counts_vector(
          pv, p_distance, pattern_counts);

      // Automatically zero initially.
      auto& previous_length = last_p_subgraph_count_length[pv];
      if (previous_length == pattern_counts.size()) {
        // No further change for PV, so finish testing this PV.
        continue;
      }
      previous_length = pattern_counts.size();
      get_tv_to_erase(
          p_distance, tv_to_erase, pattern_counts, domain, t_counts_calculator);
      if (tv_to_erase.size() == domain.size()) {
        return false;
      }
      for (auto tv : tv_to_erase) {
        domain.erase(tv);
      }

      // Try a simple alldiff reduction also,
      // but not a full propagation.
      if (domain.size() == 1 && assigned_pv.count(pv) == 0) {
        assigned_pv.insert(pv);
        const VertexWSM tv = *domain.cbegin();
        for (auto& entry_other : possible_assignments) {
          if (pv == entry_other.first) {
            continue;
          }
          auto& domain_other = entry_other.second;
          if (domain_other.erase(tv) > 0 && domain_other.empty()) {
            return false;
          }
        }
      }
    }
  }
  return true;
}


}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
