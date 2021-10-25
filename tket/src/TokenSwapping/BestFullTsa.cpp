#include "BestFullTsa.hpp"

#include "DistancesFromArchitecture.hpp"
#include "NeighboursFromArchitecture.hpp"
#include "RiverFlowPathFinder.hpp"
#include "TableLookup/VertexMapResizing.hpp"

namespace tket {
namespace tsa_internal {

BestFullTsa::BestFullTsa() { m_name = "BestFullTsa"; }

HybridTsa00& BestFullTsa::get_hybrid_tsa_for_testing() { return m_hybrid_tsa; }

void BestFullTsa::append_partial_solution(
    SwapList& swaps, VertexMapping& vertex_mapping,
    const ArchitectureMapping& arch_mapping) {
  DistancesFromArchitecture distances(arch_mapping);
  NeighboursFromArchitecture neighbours(arch_mapping);
  RiverFlowPathFinder path_finder(distances, neighbours, m_rng);
  m_rng.set_seed();
  append_partial_solution(
      swaps, vertex_mapping, distances, neighbours, path_finder);
}

void BestFullTsa::append_partial_solution(
    SwapList& swaps, VertexMapping& vertex_mapping,
    DistancesInterface& distances, NeighboursInterface& neighbours,
    PathFinderInterface& path_finder) {
  auto vm_copy = vertex_mapping;

  m_hybrid_tsa.append_partial_solution(
      swaps, vm_copy, distances, neighbours, path_finder);

  // Still subject to experimentation, but this seems the best
  m_swap_list_optimiser.optimise_pass_with_zero_travel(swaps);
  m_swap_list_optimiser.optimise_pass_with_token_tracking(swaps);
  m_swap_list_optimiser.optimise_pass_remove_empty_swaps(swaps, vertex_mapping);
  m_swap_list_optimiser.full_optimise(swaps, vertex_mapping);

  VertexMapResizing map_resizing(neighbours);
  std::set<size_t> vertices_with_tokens_at_start;
  for (const auto& entry : vertex_mapping) {
    vertices_with_tokens_at_start.insert(entry.first);
  }
  m_table_optimiser.optimise(
      vertices_with_tokens_at_start, map_resizing, swaps,
      m_swap_list_optimiser);
}

}  // namespace tsa_internal
}  // namespace tket
