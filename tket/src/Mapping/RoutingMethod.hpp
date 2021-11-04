#ifndef _TKET_RoutingMethod_H_
#define _TKET_RoutingMethod_H_

#include "Mapping/MappingFrontier.hpp"
#include "Utils/Json.hpp"

namespace tket {

class RoutingMethod {
 public:
  RoutingMethod(){};
  virtual ~RoutingMethod() {}
  /**
   * check_method returns true if held method can route given circuit.
   * This is completed by converting boundary subcircuit to a Circuit object
   * which is then passed to check_subcircuit_ as defined in constructor.
   *
   * Overloded parameter mapping_frontier contains boundary of gates to be
   * checked for method.
   * Overloaded parameter architecture is the architecture method works with
   * if permitted.
   * @return true if method can route subcircuit, false if not
   */
  virtual bool check_method(
      const std::shared_ptr<MappingFrontier>& /*mapping_frontier*/,
      const ArchitecturePtr& /*architecture*/) const {
    return false;
  }

  /**
   * routing_method modifies circuit held in mapping_frontier with gates for the
   * purpose of moving circuit closer to one physically permitted by given
   * architecture. Returns new initial mapping of qubits incase permutation via
   * swap network is then required, or new ancilla qubits are added.
   * This is completed by converting boundaty subcircuit in mapping frontier to
   * a Circuit object which is then passed to route_subcircuit_ as defined in
   * the constructor.
   *
   * Overloaded parameter mapping_frontier contains boundary of routed/unrouted
   * circuit for modifying.
   * Overloaded parameter architecture provides physical constraints
   * @return Logical to Physical mapping at boundary due to modification.
   *
   */
  virtual unit_map_t routing_method(
      std::shared_ptr<MappingFrontier>& /*mapping_frontier*/,
      const ArchitecturePtr& /*architecture*/) const {
    return {};
  }

  virtual nlohmann::json serialize() const {
    throw JsonError(
        "JSON serialization not implemented for given RoutingMethod.");
  }
};

typedef std::reference_wrapper<RoutingMethod> RoutingMethodWrapper;
}  // namespace tket

#endif