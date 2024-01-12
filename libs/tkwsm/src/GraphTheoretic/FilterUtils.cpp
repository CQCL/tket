// Copyright 2019-2024 Cambridge Quantum Computing
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

#include "tkwsm/GraphTheoretic/FilterUtils.hpp"

#include <algorithm>
#include <tkassert/Assert.hpp>

namespace tket {
namespace WeightedSubgraphMonomorphism {

bool FilterUtils::compatible_sorted_degree_sequences(
    const std::vector<std::size_t>& pattern_v_deg_seq,
    const std::vector<std::size_t>& target_v_deg_seq) {
  // If a pattern vertex v_P has deg sequence  p(0) <= p(1) <= p(2) <= ... <=
  // p(n) in the pattern graph, and target vertex v_T has deg seq   t(0) <= t(1)
  // <= ... <= t(m) in the target graph, then:
  //
  //  v_P can be mapped to v_T (at least locally) <==>
  //
  //  there exists injective F : {0,1,...,n} -> {0,1,...,m} such that p(i) <=
  //  t(F(i)) for all i
  //     [obvious by defn: F IS the vertex mapping]
  //
  //   <==>
  //
  // (p(i)) can be placed inside a dominating subsequence of (t(j)),
  //  i.e. we can make F be INCREASING
  //
  //    [Pr.: <== is trivial; for ==>: can do bubble sort:
  //    let F be given, inj but maybe not increasing. If F(i+1) < F(i), then
  //    p(i) <= p(i+1) <= t(F(i+1)) <= t(F(i)).
  //    Thus both p(i), p(i+1) are dominated by t(F(i+1)), t(F(i))
  //    so we can just swap F(i+1), F(i) and F remains valid].
  //
  //  <==>
  //
  //  we can construct a solution greedily, i.e. DEFINE G(0) = min { j : p(0) <=
  //  t(j)}, and G(i+1) = min { j > G(i) : p(i+1) <= t(j)}, and just set F=G
  //  (assuming that G(0),...,G(n) exist).
  //
  //    [Pr.: <== triv: if the alg succeeds, G is strictly increasing by
  //    construction.
  //    ==>: let F incr be given. Start trying to construct G by the alg.
  //    We have G(0) <= F(0). If at stage k, G(k) has been constructed AND G(i)
  //    <= F(i) for all i <= k, and k<n, consider trying to find G(k+1). Let S =
  //    { j > G(k) : p(k+1) <= t(j)}. By definition, F(k+1) is in S, since
  //    F(k+1) > F(k) >= G(k). Thus S is nonempty, so a value G(k+1) > G(k)
  //    exists. By construction G(k+1) is the MIN possible solution, so G(k+1)
  //    <= F(k+1). Thus by induction, G will be fully constructed and G(i) <=
  //    F(i) for all i. Now, by defn p(k) <= t(G(k)) for all k, so G really does
  //    work.]
  //
  //
  if (pattern_v_deg_seq.size() > target_v_deg_seq.size()) {
    return false;
  }
  if (pattern_v_deg_seq.empty()) {
    return true;
  }

  // This will remain the iterator to the first elem
  // to which we can map.
  auto start_t_citer = target_v_deg_seq.cbegin();

  for (std::size_t next_p_index_to_dominate = 0;
       next_p_index_to_dominate < pattern_v_deg_seq.size();
       ++next_p_index_to_dominate) {
    const std::size_t& p_degree = pattern_v_deg_seq[next_p_index_to_dominate];
    const auto t_citer =
        std::lower_bound(start_t_citer, target_v_deg_seq.cend(), p_degree);
    if (t_citer == target_v_deg_seq.cend()) {
      // Every "t" degree is < the "p" degree.
      return false;
    }
    // t_citer exactly points towards G(k).
    // How many t-values are still available (including this G(k) we're about to
    // assign?)
    const std::size_t remaining_t_values = target_v_deg_seq.cend() - t_citer;
    const std::size_t remaining_p_values =
        pattern_v_deg_seq.size() - next_p_index_to_dominate;
    if (remaining_t_values < remaining_p_values) {
      return false;
    }
    // A slight optim is possible: if remaining_t_values == remaining_p_values
    // we can just go through linearly, no binary search.
    // We could also have been fancy and searched up to somewhere before
    // target_v_deg_seq.cend() earlier.
    start_t_citer = t_citer + 1;
  }
  return true;
}

bool FilterUtils::compatible_sorted_degree_counts(
    const DegreeCounts& degree_counts1, const DegreeCounts& degree_counts2) {
  if (degree_counts1.empty()) {
    return true;
  }
  auto counts_to_satisfy = degree_counts1.back();
  TKET_ASSERT(counts_to_satisfy.first >= 1);
  TKET_ASSERT(counts_to_satisfy.second >= 1);
  if (degree_counts2.empty()) {
    return false;
  }
  auto next_counts_sink = degree_counts2.back();
  unsigned index1 = degree_counts1.size() - 1;
  unsigned index2 = degree_counts2.size() - 1;

  // Break out when we've cleared all the pattern vertices to be matched.
  for (;;) {
    while (next_counts_sink.first < counts_to_satisfy.first) {
      if (index2 == 0) {
        // We can't move down any further.
        return false;
      }
      --index2;
      next_counts_sink = degree_counts2[index2];
      TKET_ASSERT(next_counts_sink.first >= 1);
      TKET_ASSERT(next_counts_sink.second >= 1);
    }
    if (counts_to_satisfy.second <= next_counts_sink.second) {
      next_counts_sink.second -= counts_to_satisfy.second;
      if (index1 == 0) {
        break;
      }
      --index1;
      counts_to_satisfy = degree_counts1[index1];
      TKET_ASSERT(counts_to_satisfy.first >= 1);
      TKET_ASSERT(counts_to_satisfy.second >= 1);
      continue;
    }
    counts_to_satisfy.second -= next_counts_sink.second;
    if (index2 == 0) {
      // We can't move down any further.
      return false;
    }
    --index2;
    next_counts_sink = degree_counts2[index2];
    TKET_ASSERT(next_counts_sink.first >= 1);
    TKET_ASSERT(next_counts_sink.second >= 1);
  }
  return true;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
