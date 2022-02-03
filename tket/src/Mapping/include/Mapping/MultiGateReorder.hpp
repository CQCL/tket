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
   * @param max_depth Maximum number of layers of gates checked for commutation.
   * @param max_size Maximum number of gates checked for commutation.
   */
  void solve(unsigned max_depth, unsigned max_size);

 private:
  // Architecture all new physical operations must respect
  ArchitecturePtr architecture_;
  std::shared_ptr<MappingFrontier> mapping_frontier_;
  EdgeVec u_frontier_edges_;
};

class MultiGateReorderRoutingMethod : public RoutingMethod {
 public:
  /**
   * Checking and Routing methods redefined using MultiGateReorder.
   * @param _max_depth Maximum number of layers of gates checked for
   * commutation.
   * @param _max_size Maximum number of gates checked for commutation.
   */
  MultiGateReorderRoutingMethod(
      unsigned _max_depth = 10, unsigned _max_size = 10);

  /**
   * @return true if method can route subcircuit, false if not
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

  nlohmann::json serialize() const override;

  static MultiGateReorderRoutingMethod deserialize(const nlohmann::json& j);

  /**
   * @return Maximum number of layers of gates checked for commutation.
   */
  unsigned get_max_depth() const;

  /**
   * @return Maximum number of gates checked for commutation.
   */
  unsigned get_max_size() const;

 private:
  unsigned max_depth_;
  unsigned max_size_;
};

}  // namespace tket

#endif