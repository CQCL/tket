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

#pragma once
#include <optional>
#include <string>
#include <tkrng/RNG.hpp>

#include "tkwsm/GraphTheoretic/GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {
class NeighboursData;
namespace InitialPlacement {

/** For storing the current solution and making random jumps;
 * but note that we don't actually store the scalar product,
 * we only need to calculate the DIFFERENCE from doing a jump.
 */
class SolutionJumper {
 public:
  SolutionJumper(
      const NeighboursData& pattern_ndata, const NeighboursData& target_ndata,
      // Any target edge not mentioned in target_ndata is given this weight.
      WeightWSM implicit_target_weight);

  // Element[PV] gives the current TV for PV.
  // (We're using contiguous labels of course, like NeighboursData).
  const std::vector<unsigned>& get_assignments() const;

  const NeighboursData& get_pattern_ndata() const;

  const NeighboursData& get_target_ndata() const;

  /** This is for resetting to a completely new solution.
   * Element[PV] = TV, where the assignments are PV->TV.
   * They are stored internally, and of the correct size.
   * The CALLER is responsible for
   * overwriting this vector with valid data
   * (the caller must not resize, and must fill with a valid
   * list of TVs).
   */
  std::vector<unsigned>& get_assignments_to_overwrite();

  /** Immediately after calling get_assignments_to_overwrite(),
   * and overwriting the values, call this function to update.
   * This also returns the scalar product of the new solution.
   * (which the caller will surely want to know).
   */
  WeightWSM reset_and_get_new_scalar_product();

  /** Suppose that currently PV1 -> TV1, and we want to change to PV1->TV2.
   * If TV2 is unoccupied, nothing else needs to change.
   * If instead TV2 is occupied (by PV2), i.e. we already have PV2->TV2,
   * we would also change to PV2->TV1, i.e. a straight target vertex swap.
   * If the scalar product would decrease, perform the move, and return
   * the amount by which it decreases.
   */
  std::optional<WeightWSM> perform_move_and_get_scalar_product_decrease(
      unsigned pv1, unsigned tv2, WeightWSM minimum_decrease = 1);

  void check_validity() const;

 private:
  const NeighboursData& m_pattern_ndata;
  const NeighboursData& m_target_ndata;
  const WeightWSM m_implicit_target_weight;

  // Element[PV] = TV, where currently we have assigned PV->TV.
  std::vector<unsigned> m_assigned_target_vertices;

  // The reverse of m_assigned_target_vertices.
  // Element[TV] gives the current PV such that we've assigned PV->TV.
  // HOWEVER, some TV might not have any PV assigned to them,
  // so they have a dummy value, larger than any valid PV.
  std::vector<unsigned> m_source_pattern_vertices;

  // For each PV, element[PV] gives the sum of (p_weight)*(t_weight)
  // over all p-edges joined to PV, with the current assignments.
  // Will be set to a 0 if it needs recalculating
  // (i.e., one of the neighbouring PVs has changed).
  mutable std::vector<WeightWSM> m_scalar_product_contributions;

  WeightWSM get_current_scalar_product_contribution(unsigned pv) const;

  /* What happens when we currently have PV->TV, and we want to reassign
    PV->TV'? Assuming that TV, TV' are different (otherwise there is no
    change to make), there are THREE possible cases:

    Case A (the usual case, "independent"): there is currently another PV'
        with PV'->TV', but there is NO edge PV--PV'.
        Thus, our trial jump will actually consist of a swap:
        PV->TV' and PV'->TV, and then the old and new scalar product
        contributions for PV, PV' don't affect each other.

    Case B ("empty image"): currently, there is NO p-vertex assigned to TV'.
        This can be viewed as the same as case A, if we invent a new isolated
        dummy vertex PV' with PV'->TV' currently: because it is isolated,
        it will never have a nonzero scalar product contribution.

    Case C ("dependent"): there is currently another PV'
        with PV'->TV', but the edge PV--PV' DOES exist.
        This is like Case A, but we have to be careful: after the trial jump
        PV->TV' and PV'->TV, BOTH new contributions would include the SAME
        contribution from the edge PV--PV'.
        This is not a problem, we just have to take it into account.

        Notice that if we subtract the CURRENT contributions from PV, PV'
        and then add the NEW contributions, the contribution from edge PV--PV'
        (i.e., weight(PV, PV') * weight(TV, TV'), the same before AND after
        the jump) is subtracted twice,
        then added twice.
        So the net effect is the same as if PV--PV' did not exist,
        so we don't need to do anything other than calculate
        the NEW contributions correctly.
    */

  struct HypotheticalScalarProductContribution {
    // This would be the new contribution from PV,
    // if the jump were to take place. In case C,
    // this includes the contribution from
    // the edge PV--PV'.
    WeightWSM contribution;

    // IF we are in case C (the "dependent" case),
    // this gives PV' which is adjacent to the original PV,
    // and currently
    std::optional<unsigned> case_c_other_pv_opt;
  };

  HypotheticalScalarProductContribution
  get_hypothetical_scalar_product_contribution(unsigned pv, unsigned tv) const;

  // Case C is required not to occur.
  WeightWSM get_hypothetical_scalar_product_contribution_disallowing_case_c(
      unsigned pv, unsigned tv) const;

  // If some PV is about to be reassigned, mark all its neighbour PV
  // contributions as invalid; they will need to be recalculated (lazily).
  void invalidate_neighbour_contributions(unsigned pv);

  // Returns the edge weight in all case, i.e. returns
  // the implicit edge weight if no explicit edge exists.
  WeightWSM get_target_edge_weight(unsigned tv1, unsigned tv2) const;
};

}  // namespace InitialPlacement
}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
