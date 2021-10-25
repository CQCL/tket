#include "HybridTsa00.hpp"

#include "TSAUtils/DistanceFunctions.hpp"
#include "Utils/Assert.hpp"

;
using std::vector;

namespace tket {
namespace tsa_internal {

HybridTsa00::HybridTsa00() {
  m_name = "HybridTSA_00";
  m_trivial_tsa.set(TrivialTSA::Options::BREAK_AFTER_PROGRESS);
}

CyclesPartialTsa& HybridTsa00::get_cycles_tsa_for_testing() {
  return m_cycles_tsa;
}

TrivialTSA& HybridTsa00::get_trivial_tsa_for_testing() { return m_trivial_tsa; }

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
