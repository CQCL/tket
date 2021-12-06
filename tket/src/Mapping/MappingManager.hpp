#ifndef _TKET_MappingManager_H_
#define _TKET_MappingManager_H_

#include "Architecture/Architecture.hpp"
#include "Circuit/Circuit.hpp"
#include "Mapping/RoutingMethod.hpp"
#include "Utils/UnitID.hpp"

namespace tket {

// list of error types to throw out
class MappingManagerError : public std::logic_error {
 public:
  explicit MappingManagerError(const std::string& message)
      : std::logic_error(message) {}
};

typedef ArchitecturePtr ArchitecturePtr;

class MappingManager {
 public:
  /* Mapping Manager Constructor */
  // MappingManager object defined by Architecture initialised with
  MappingManager(const ArchitecturePtr& _architecture);

  /**
   * route_circuit
   * Referenced Circuit modified such that all multi-qubit gates are permitted
   * by this->architecture_ RoutingIncompability thrown if Circuit has more
   * logical qubits than Architecture has physical qubits RoutingIncompability
   * thrown if Circuit has a gate of OpType not in Architecture's permitted
   * OpTypes
   *
   * @param circuit Circuit to be routed
   * @param routing_methods Ranked RoutingMethod objects to use for routing
   * segments.
   *
   * @return True if circuit is modified
   */
  bool route_circuit(
      Circuit& circuit,
      const std::vector<RoutingMethodPtr>& routing_methods) const;

 private:
  ArchitecturePtr architecture_;
};
}  // namespace tket

#endif