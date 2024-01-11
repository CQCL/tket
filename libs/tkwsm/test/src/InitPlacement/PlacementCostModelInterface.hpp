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

#include <string>
#include <tkwsm/GraphTheoretic/GeneralStructs.hpp>

namespace tket {
namespace WeightedSubgraphMonomorphism {
namespace InitialPlacement {
namespace tests {

/** Just for testing; we want reaonsably large graphs where
 * sensible token swapping is easy to work out.
 *
 * Thus we get a reasonable simplified model for the cost of applying
 * a sequence of gates, once we've got an initial qubit placement.
 *
 * However, this is only for very quick sanity tests; a full benchmark
 * should be carried out with many different graphs,
 * token swapping strategies, etc. (Even GIVEN an initial placement,
 * and GIVEN a sequence of gates to apply, there may be many
 * possible reorderings from parallel gates;
 * so there are probably many possible solutions, even for e.g. trees
 * where shortest paths are unique (although smallest WEIGHT paths are not).
 *
 * Initialising with the detailed graph edges, weights, etc. should be done
 * via a constructor.
 */
class PlacementCostModelInterface {
 public:
  /** Given the initial PV->TV assignments, set up all the tokens etc.
   * ready to do token swapping and cost calculations.
   * We don't impose restrictions on the PV or TV (apart from them being
   * "sensibly sized", as we use larger values for debug purposes).
   * The KEY is PV the VALUE is TV.
   */
  void initialise_with_qubit_placement(
      const std::vector<std::pair<VertexWSM, VertexWSM>>& placement);

  // Wherever the tokens corresponding to pv1, pv2 currently are,
  // move them to be adjacent,
  // so that a gate can be applied between them.
  // Return the total cost
  // of the gate AND the swaps (remembering that a single swap may count as
  // multiple primitive gates, according to
  // m_number_of_primitive_gates_in_swap).
  // We DON'T then move them back!
  WeightWSM do_token_swapping_and_apply_gate(VertexWSM pv1, VertexWSM pv2);

  // Call only after "initialise_with_qubit_placement".
  // This performs the given gates, in that order, and returns the total cost.
  WeightWSM get_cost(const std::vector<std::pair<VertexWSM, VertexWSM>>& gates);

  // Simply runs forever, printing out best results found as it goes.
  // A very crude test function to try to find good solutions.
  void try_random_placements(
      const std::vector<std::pair<VertexWSM, VertexWSM>>& gates);

  // KEY: the PV (which we can also think of as the token);
  // VALUE: the TV where it has currently been assigned to
  const std::map<VertexWSM, VertexWSM>& get_current_placement() const;

  // KEY: a TV.
  // VALUE: the PV (i.e., token - we think of the initial placement as defining
  //    the tokens - a token number is exactly the PV where it came from)
  //    of the token currently sitting at TV.
  //    For SPEED, note that some values might be "dummy" values,
  //    larger than any valid PV.
  //    This is because repeatedly inserting and erasing keys in a map
  //    can be relatively a lot slower than altering an existing value
  //    at an existing key.
  const std::map<VertexWSM, VertexWSM>& get_current_tokens() const;

  // Does a "require" test that the vertices and tokens are all correct.
  void require_valid() const;

  // For testing, calculate Path(v1,v2) and return it as a string.
  std::string get_path_str(VertexWSM v1, VertexWSM v2) const;

  // Calculate the edges and weights in standard format.
  virtual GraphEdgeWeights get_graph_data() const = 0;

  virtual ~PlacementCostModelInterface();

 protected:
  WeightWSM m_number_of_primitive_gates_in_swap = 3;
  VertexWSM m_invalid_token = 100000;

  // KEY: PV (i.e., a token)
  // VALUE: TV (i.e., a vertex in this graph).
  std::map<VertexWSM, VertexWSM> m_current_placement;

  // KEY: TV
  // VALUE: PV (or a dummy large value).
  std::map<VertexWSM, VertexWSM> m_current_tokens;

  // The path is a list of (vertex, edge weight pairs),
  // giving an actual path to be traversed.
  // In each element, the second value is the edge weight
  // to move from the PREVIOUS vertex in the path;
  // thus element[0] always starts with weight 0.
  using Path = std::vector<std::pair<VertexWSM, WeightWSM>>;

  // If it's a path from v1 to v2, change into a path from v2 to v1.
  static void reverse_path(Path& path);

  static WeightWSM get_total_weight(const Path& path);

  // For the path [v0 v1 v2 ... v(n)], get information about
  // the edge weights for  v(i)--v(i+1),  0 <= i <= n-1,
  // and also return the index i such that v(i)--v(i+1) has highest weight.
  //
  // Ironically, when applying a gate by swapping tokens along a path,
  // we should always let the actual gate have the HIGHEST edge weight
  //
  // (i.e., ** we deliberately choose the WORST possible fidelity
  // for our gate!!!!! **)
  //
  // This is because we MUST use each edge in the path exactly once
  // (either for a SWAP, or the actual gate);
  // but swaps cost more than primitive gates, so we should use
  // the worst fidelity for our primitive gate, to reduce the total cost.
  struct WeightsData {
    // The sum of the weights of those edges along which we apply a swap
    // (so, NOT including the edge weight along which we apply
    // the 2-qubit gate, and also NOT including the swap weight multiplier).
    WeightWSM sum_of_swap_edge_weights;
    WeightWSM highest_edge_weight;

    // Set to the value of i such that Weight(v[i], v[i+1]) is
    // (jointly) the maximum edge weight in the path.
    unsigned index_of_largest_edge_weight;

    explicit WeightsData(const Path& path);
  };

  // Actually move the tokens around along the path,
  // and the weights data (since, it depends on which edge along the path
  // is chosen for the gate).
  void enact_swaps(const Path& path, const WeightsData& weights_data);

  // Return some valid path, stored and cached internally somehow.
  // To bring tokens on v1,v2 together, We will do token swapping
  // along this path.
  virtual const Path& get_path_to_use(
      VertexWSM vertex1, VertexWSM vertex2) const = 0;
};

}  // namespace tests
}  // namespace InitialPlacement
}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
