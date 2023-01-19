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

#pragma once
#include <cstdint>
#include <vector>

namespace tket {
namespace WeightedSubgraphMonomorphism {

/** Functions useful for filtering, i.e. deciding whether to erase target
 * vertices from domains, etc. etc., based on graph theoretic properties.
 */
struct FilterUtils {
  /** Given sorted degree sequences (i.e., the vertex degrees of
   * the neighbours) of vertices P_v and T_v, is it maybe possible to
   * map P_v to T_v, based upon the sequences?
   * @param pattern_v_deg_seq The increasing vertex degrees of the neighbours of
   * PV
   * @param target_v_deg_seq The increasing vertex degrees of the neighbours of
   * TV
   * @return True if the target sequence dominates the pattern sequence, i.e. in
   * the absence of other graph theoretic information, there is no obstruction
   * to mapping PV->TV.
   */
  static bool compatible_sorted_degree_sequences(
      const std::vector<std::size_t>& pattern_v_deg_seq,
      const std::vector<std::size_t>& target_v_deg_seq);

  // In each element, FIRST is a vertex degree;
  // SECOND is the number of vertices with that degree.
  // Should be sorted lexicographically.
  // We require always that the first, second are both >= 1.
  typedef std::vector<std::pair<std::size_t, std::size_t>> DegreeCounts;

  /** Really the same function as "compatible_sorted_degree_sequences",
   * just with the data in a different format.
   */
  static bool compatible_sorted_degree_counts(
      const DegreeCounts& degree_counts1, const DegreeCounts& degree_counts2);
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
