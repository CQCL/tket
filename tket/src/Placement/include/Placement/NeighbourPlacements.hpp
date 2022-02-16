#pragma once

#include <map>

#include "Placement.hpp"
#include "TokenSwapping/SwapFunctions.hpp"
#include "Utils/BiMapHeaders.hpp"

namespace tket {

/**
 * @brief Given a placement map generates `n` nearby placement maps.
 *
 * Based on an architecture and a placement map, generates random
 * placements that can be achieved with `m` swaps.
 *
 * Optionally uses token swapping optimisations to try to ensure
 * that the generated placements cannot be obtained in less than `m`
 * swaps, but this cannot be guaranteed.
 */
class NeighbourPlacements {
 public:
  using Swap = tsa_internal::Swap;
  using SwapVec = std::vector<Swap>;
  using SwapList = tsa_internal::SwapList;
  using NodeSwap = std::pair<Node, Node>;
  using NodeSwapVec = std::vector<NodeSwap>;
  struct Result {
    qubit_mapping_t map;
    NodeSwapVec swaps;
  };
  using ResultVec = std::vector<Result>;

  /**
   * @brief Construct a new Swap Placement object.
   *
   * @param arc The architecture defining the allowed swaps.
   * @param init_map The initial Qubit => Node map.
   */
  NeighbourPlacements(const Architecture& arc, const qubit_mapping_t& init_map);

  /**
   * @brief Generate `n` placement maps using up to `dist` swaps for each map
   *
   * @param dist The number of swaps allowed on the architecture.
   * @param n The number of placement maps to generate (default n=1).
   * @param optimise Simplify the generated swap sequences (default true).
   * @param seed Seed for random number generator (default seed=0).
   * @return ResultVec A vector of the generated maps and swaps
   */
  ResultVec get(
      unsigned dist, unsigned n = 1, bool optimise = true, unsigned seed = 0);
 private:
  // generate a single Result
  Result gen_result(unsigned dist, bool optimise = true);
  // generate a swap list with `dist` swaps
  SwapList gen_swap_list(unsigned dist);
  // apply swap list to init_map and return new placement map
  Result convert_to_res(const SwapVec& swaps);

  Architecture arc_;
  qubit_mapping_t init_map_;
  boost::bimap<unsigned, Node> u_to_node_;
};

}  // namespace tket