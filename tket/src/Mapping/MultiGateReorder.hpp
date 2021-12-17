#ifndef _TKET_MultiGateReorder_H_
#define _TKET_MultiGateReorder_H_

#include "Mapping/MappingFrontier.hpp"
#include "Mapping/RoutingMethod.hpp"

namespace tket {

class MultiGateReorder {
 public:
  /**
   * Class Constructor
   * @param _architecture Architecture object added operations must respect
   * @param _mapping_frontier Contains Circuit object to be modified
   */
  MultiGateReorder(
      const ArchitecturePtr& _architecture,
      std::shared_ptr<MappingFrontier>& _mapping_frontier);

  /**
   * Try to commute any multi-qubit gates to the quantum frontier
   */
  void solve();

 private:
  // Check an edge is in the frontier
  bool edge_in_frontier(const Edge& edge);
  // Try to commute a multiqubit gate to the quantum frontier
  bool try_commute_multi_to_front(const Vertex& vert);
  // Architecture all new physical operations must respect
  ArchitecturePtr architecture_;
  // Contains circuit for finding SWAP from and non-routed/routed boundary
  std::shared_ptr<MappingFrontier> mapping_frontier_;
  EdgeVec u_frontier_edges_;
};

class MultiGateReorderRoutingMethod : public RoutingMethod {
 public:
  /**
   * Checking and Routing methods redefined using MultiGateReorder.
   *
   */
  MultiGateReorderRoutingMethod(){};

  /**
   * @return true if method can route subcircuit, false if not
   */
  bool check_method(
      const std::shared_ptr<MappingFrontier>& /*mapping_frontier*/,
      const ArchitecturePtr& /*architecture*/) const;

  /**
   * @param mapping_frontier Contains boundary of routed/unrouted circuit for
   * modifying
   * @param architecture Architecture providing physical constraints
   * @return Logical to Physical mapping at boundary due to modification.
   *
   */
  unit_map_t routing_method(
      std::shared_ptr<MappingFrontier>& mapping_frontier,
      const ArchitecturePtr& architecture) const;
};

}  // namespace tket

#endif