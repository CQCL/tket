#ifndef _TKET_RoutingMethodCircuit_H_
#define _TKET_RoutingMethodCircuit_H_

#include "Mapping/RoutingMethod.hpp"

namespace tket {

class RoutingMethodCircuit : public RoutingMethod {
 public:
  virtual ~RoutingMethodCircuit() {}
  /**
   * RoutingMethodCircuit objects hold methods for partially routing subcircuits
   * in the incremental routing of full circuits.
   *
   * @param _route_subcircuit Function ptr for partial routing method
   * @param _check_subcircuit Function ptr for confirming if method sufficient
   * @param _max_size Max number of gates in partial routing circuit
   * @param _max_depth Max depth of partial routing circuit
   */
  RoutingMethodCircuit(
      const std::function<std::tuple<Circuit, unit_map_t, unit_map_t>(
          const Circuit&, const ArchitecturePtr&)>
          _route_subcircuit,
      const std::function<bool(const Circuit&, const ArchitecturePtr&)>
          _check_subcircuit,
      unsigned _max_size, unsigned _max_depth);

  /**
   * @param mapping_frontier Contains boundary of gates to be checked for method
   * @param architecture Architecture method would work with if permitted
   * @return true if method can route subcircuit, false if not
   */
  bool check_method(
      const std::shared_ptr<MappingFrontier>& mapping_frontier,
      const ArchitecturePtr& architecture) const;

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

 private:
  const std::function<std::tuple<Circuit, unit_map_t, unit_map_t>(
      const Circuit&, const ArchitecturePtr&)>
      route_subcircuit_;
  const std::function<bool(const Circuit&, const ArchitecturePtr&)>
      check_subcircuit_;
  unsigned max_size_, max_depth_;
};
}  // namespace tket

#endif
