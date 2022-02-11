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

#include "HybridTsa00.hpp"

#include "TokenSwapping/DistanceFunctions.hpp"
#include "Utils/Assert.hpp"

using std::vector;

namespace tket {
namespace tsa_internal {

HybridTsa00::HybridTsa00() {
  m_name = "HybridTSA_00";
  m_trivial_tsa.set(TrivialTSA::Options::BREAK_AFTER_PROGRESS);
}

// GCOVR_EXCL_START
CyclesPartialTsa& HybridTsa00::get_cycles_tsa_for_testing() {
  return m_cycles_tsa;
}

TrivialTSA& HybridTsa00::get_trivial_tsa_for_testing() { return m_trivial_tsa; }
// GCOVR_EXCL_STOP

void HybridTsa00::append_partial_solution(
    SwapList& swaps, VertexMapping& vertex_mapping,
    DistancesInterface& distances, NeighboursInterface& neighbours,
    PathFinderInterface& path_finder) {
  const auto initial_L = get_total_home_distances(vertex_mapping, distances);
  for (size_t counter = initial_L + 1; counter > 0; --counter) {
    const auto swaps_before = swaps.size();
    m_cycles_tsa.append_partial_solution(
        swaps, vertex_mapping, distances, neighbours, path_finder);

    m_trivial_tsa.append_partial_solution(
        swaps, vertex_mapping, distances, neighbours, path_finder);

    if (swaps_before == swaps.size()) {
      TKET_ASSERT(all_tokens_home(vertex_mapping));
      return;
    }
  }
  TKET_ASSERT(!"hybrid TSA termination");
}

}  // namespace tsa_internal
}  // namespace tket
