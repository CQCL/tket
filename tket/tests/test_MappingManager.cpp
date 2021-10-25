#include <catch2/catch.hpp>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "Mapping/MappingManager.hpp"

namespace tket {

SCENARIO("Test MappingManager::route_circuit") {
  Node node0("test_node", 0), node1("test_node", 1), node2("test_node", 2);
  Architecture arc({{node0, node1}, {node1, node2}});
  ArchitecturePtr shared_arc = std::make_shared<Architecture>(arc);
  MappingManager test_mm(shared_arc);
  RoutingMethod test_rm;
  std::vector<std::reference_wrapper<RoutingMethod>> test_vrm = {test_rm};
  GIVEN("More qubits than architecture has qubits.") {
    Circuit circ(5);
    REQUIRE_THROWS_AS(
        test_mm.route_circuit(circ, test_vrm), MappingManagerError);
  }
  GIVEN("Circuit unmodified.") {
    Circuit circ(2);
    REQUIRE(!test_mm.route_circuit(circ, test_vrm));
  }
  GIVEN("No method can route circuit.") {
    Circuit circ(3);
    std::vector<Qubit> qubits = circ.all_qubits();
    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[2]});
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], node0}, {qubits[1], node1}, {qubits[2], node2}};
    circ.rename_units(rename_map);
    REQUIRE_THROWS_AS(
        test_mm.route_circuit(circ, test_vrm), MappingManagerError);
  }
}
}  // namespace tket
