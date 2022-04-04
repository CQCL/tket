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

#pragma once
#include <map>
#include <optional>
#include <set>
#include <string>

#include "../Common/GeneralUtils.hpp"
#include "SearchNode.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

/** Contains raw search node data, with functions to manipulate it. */
class SearchNodeWrapper {
 public:
  SearchNodeWrapper();

  explicit SearchNodeWrapper(SearchNode node);

  /** Convenient to have a single getter function.
   * @return A const reference to the stored node.
   */
  const SearchNode& get() const;

  /** Usually, avoid this; the whole purpose of a wrapper
   * is to restrict the operations done on the stored SearchNode.
   * However, the algorithms are still subject to change, so the interface
   * is not yet stable.
   * @return A reference to the stored node.
   */
  SearchNode& get_mutable();

  /** Add the given amount to the stored total p edge weights.
   * @param extra_weight The amount to add.
   * @return For chaining, a reference to this object.
   */
  SearchNodeWrapper& add_p_edge_weights(WeightWSM extra_weight);

  /** Add the given amount to the stored total weight (scalar product).
   * @param extra_weight The amount to add.
   * @return For chaining, a reference to this object.
   */
  SearchNodeWrapper& add_scalar_product(WeightWSM extra_weight);

  /** The caller must ensure that every element of the new domain
   * is within the existing domain (not checked) - so it's a subset.
   * The caller must also ensure that the existing domain DOES exist,
   * i.e. it is not already an assignment.
   *
   * It's ASSUMED that the caller already has checked
   * that the new domain is nonempty.
   *
   * Check if the new domain has size 1,
   * updating assignments and removing the domain if it does,
   * and overwriting the domain of PV with the new domain
   * if it has size >1.
   * @param new_domain The new list of target vertices which gives the domain of
   * PV.
   * @param pv The pattern vertex, whose domain is being changed.
   * @param assignments The assignments which will be updated if the new domain
   * has size 1. If so, the domain of PV is erased completely (not just
   * reduced).
   */
  void overwrite_domain(
      const std::vector<VertexWSM>& new_domain, VertexWSM pv,
      Assignments& assignments);

  /** Remove the given target vertex from that domain,
   * and return the new size. (So, 0 means a nogood).
   * Also return 0 if there is no domain for that pv.
   * However, note that we LEAVE domains of size 1 in place;
   * but we take care of new domains of size 1,
   * by adding to assignments and
   * m_node.chosen_assignments if necessary.
   * @param pv The pattern vertex.
   * @param target_vertex The target vertex, to be removed from the domain of
   * pv.
   * @param assignments ALL current PV->TV mappings (not just those made in this
   * node).
   * @return The new size of the domain.
   */
  std::size_t remove_element_from_domain(
      VertexWSM pv, VertexWSM target_vertex, Assignments& assignments);

  /** We have a container of target graph vertices,
   * and a pattern vertex assumed to have a domain.
   * Remove all the target vertices from that domain,
   * and return the new size. (So, 0 means a nogood).
   * Also return 0 if there is no domain for that pv.
   * However, note that we LEAVE domains of size 1 in place.
   * @param pv The pattern vertex.
   * @param target_vertices A container of target vertices, all to be removed
   * from the domain of pv.
   * @param assignments ALL current PV->TV mappings (not just those made in this
   * node).
   * @return The new size of the domain.
   */
  template <class Container>
  std::size_t remove_elements_from_domain(
      VertexWSM pv, const Container& target_vertices, Assignments& assignments);

 private:
  SearchNode m_node;

  // To save unnecessary relookup, pass in the domain
  // of pv.
  std::size_t remove_element_from_domain(
      VertexWSM pv, VertexWSM target_vertex, std::set<VertexWSM>& domain,
      // ALL known PV->TV mappings at this stage.
      Assignments& assignments);
};

template <class Container>
std::size_t SearchNodeWrapper::remove_elements_from_domain(
    VertexWSM pv, const Container& target_vertices, Assignments& assignments) {
  const auto citer = m_node.pattern_v_to_possible_target_v.find(pv);
  if (citer == m_node.pattern_v_to_possible_target_v.cend() ||
      citer->second.empty()) {
    return 0;
  }
  auto& domain = citer->second;
  for (VertexWSM tv : target_vertices) {
    const auto new_size =
        remove_element_from_domain(pv, tv, domain, assignments);
    if (new_size == 0) {
      return 0;
    }
  }
  return domain.size();
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
