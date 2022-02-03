#ifndef _TKET_LexiLabelling_H_
#define _TKET_LexiLabelling_H_

#include "Mapping/RoutingMethod.hpp"
#include "Mapping/LexiRoute.hpp"

namespace tket{
    
class LabellingRoutingMethod : public RoutingMethod {
 public:
  /**
   * Checking and Routing methods redefined for dynamically assigning qubits to
   * some Architecture.
   */
  LabellingRoutingMethod();

  /**
   * @return true if method can label unlabelled qubits
   */
  bool check_method(
      const std::shared_ptr<MappingFrontier>& mapping_frontier,
      const ArchitecturePtr& architecture) const override;

  /**
   * @param mapping_frontier Contains boundary of routed/unrouted circuit for
   * modifying
   * @param architecture Architecture providing physical constraints
   * @return Logical to Physical mapping at boundary due to modification.
   *
   */
  unit_map_t routing_method(
      std::shared_ptr<MappingFrontier>& mapping_frontier,
      const ArchitecturePtr& architecture) const override;
 
 private:
    std::vector<std::vector<Node>> interacting_nodes;
}
}  // namespace tket
#endif