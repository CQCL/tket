#ifndef _TKET_TokenSwapping_HybridTsa00_H_
#define _TKET_TokenSwapping_HybridTsa00_H_

#include "CyclesPartialTsa.hpp"
#include "TrivialTSA.hpp"

namespace tket {
namespace tsa_internal {

/** A full end-to-end TSA, combining the partial cycles TSA
 *  (hopefully good) with the full "trivial" TSA (not so good).
 */
class HybridTsa00 : public PartialTsaInterface {
 public:
  HybridTsa00();

  /** For the current token configuration, calculate a sequence of swaps
   *  to move all tokens home, and append them to the given list.
   *  As this is a full TSA, it guarantees to find a solution.
   *  @param swaps The list of swaps to append to.
   *  @param vertex_mapping The current desired mapping.
   *  @param distances An object to calculate distances between vertices.
   *  @param neighbours An object to calculate adjacent vertices to any given
   * vertex.
   *  @param path_finder An object to calculate a shortest path between any
   *    pair of vertices.
   */
  virtual void append_partial_solution(
      SwapList& swaps, VertexMapping& vertex_mapping,
      DistancesInterface& distances, NeighboursInterface& neighbours,
      PathFinderInterface& path_finder) override;

  /** Only for experiments; will be removed again
   *  once the best parameter combinations are found!
   *  @return A reference to the internal TSA object, to change parameters.
   */
  CyclesPartialTsa& get_cycles_tsa_for_testing();

  /** Temporary; only for experiments!
   *  @return A reference to the internal TSA object, to change parameters.
   */
  TrivialTSA& get_trivial_tsa_for_testing();

 private:
  CyclesPartialTsa m_cycles_tsa;
  TrivialTSA m_trivial_tsa;
};

}  // namespace tsa_internal
}  // namespace tket
#endif
