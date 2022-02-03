#ifndef _TKET_MultiGateReorder_H_
#define _TKET_MultiGateReorder_H_

#include "Mapping/MappingFrontier.hpp"
#include "Mapping/RoutingMethod.hpp"

namespace tket {

class BoxDecomposition {
 public:
  /**
   * Class Constructor
   * @param _architecture Architecture object added operations must respect
   * @param _mapping_frontier Contains Circuit object to be modified
   */
  BoxDecomposition(
      const ArchitecturePtr& _architecture,
      std::shared_ptr<MappingFrontier>& _mapping_frontier);

  /**
   * Decompose any boxes on the frontier
   */
  void solve();

 private:
  // Architecture all new physical operations must respect
  ArchitecturePtr architecture_;
  std::shared_ptr<MappingFrontier> mapping_frontier_;
};

class BoxDecompositionRoutingMethod : public RoutingMethod {
 public:
  /**
   * Decompose any boxes on the frontier
   */
  BoxDecompositionRoutingMethod();

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