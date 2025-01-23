// Copyright Quantinuum
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
#include <optional>
#include <tkrng/RNG.hpp>

#include "FastRandomBits.hpp"
#include "MonteCarloManager.hpp"
#include "SolutionJumper.hpp"
#include "tkwsm/GraphTheoretic/GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {
class NeighboursData;
namespace InitialPlacement {

/** Try to find a good WSM solution for an (implicit) complete target graph,
 * using randomness.
 */
class MonteCarloCompleteTargetSolution {
 public:
  /** This is all for the newly labelled vertices
   * (so, contiguous vertices {0,1,2,...,N}).
   * Upon construction, runs through to completion.
   */
  MonteCarloCompleteTargetSolution(
      const NeighboursData& pattern_ndata, const NeighboursData& target_ndata,
      // Any target edge not mentioned in target_ndata is given this weight.
      WeightWSM implicit_target_weight,
      // If left unspecified, i.e. set to zero,
      // will choose a reasonable default.
      unsigned max_iterations = 0);

  const std::vector<unsigned>& get_best_assignments() const;

  /** Return the scalar product of the best solution found so far. */
  WeightWSM get_best_scalar_product() const;

  /** The number of iterations used. */
  unsigned iterations() const;

 private:
  const WeightWSM m_implicit_target_weight;
  RNG m_rng;
  FastRandomBits m_fast_random_bits;
  unsigned m_number_of_random_bits;
  MonteCarloManager m_manager;

  unsigned m_iterations;
  unsigned m_max_iterations;

  SolutionJumper m_solution_jumper;

  // Used to generate random PV->TV mappings.
  // In element[PV]:
  // The FIRST is a collection of random bits, used to sort (to generate
  // random permutations quickly).
  // The SECOND is the actual target vertex TV to which PV is assigned.
  std::vector<std::pair<std::uint64_t, unsigned>> m_random_bits_and_tv;

  // Element[PV] is the TV it is assigned to.
  std::vector<unsigned> m_best_assignments;

  WeightWSM m_best_scalar_product;

  // The current assignments are stored within m_solution_jumper
  WeightWSM m_current_scalar_product;

  void reset_target_vertices();

  bool new_solution_is_record_breaker();
};

}  // namespace InitialPlacement
}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
