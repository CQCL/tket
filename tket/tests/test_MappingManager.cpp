#include <catch2/catch.hpp>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "Mapping/MappingManager.hpp"

namespace tket {

class TokenSwappingTester : public RoutingMethod {
 public:
  TokenSwappingTester(){};

  bool check_method(
      const std::shared_ptr<MappingFrontier>& /*mapping_frontier*/,
      const ArchitecturePtr& /*architecture*/) const {
    return true;
  }

  /**
   * @param mapping_frontier Contains boundary of routed/unrouted circuit for
   * modifying
   * @param architecture Architecture providing physical constraints
   * @return Logical to Physical mapping at boundary due to modification.
   *
   */
  unit_map_t routing_method(
      std::shared_ptr<MappingFrontier>& /*mapping_frontier*/,
      const ArchitecturePtr& /*architecture*/) const {
    Node node0("test_node", 0), node1("test_node", 1), node2("test_node", 2);
    return {{node0, node1}, {node1, node2}, {node2, node0}};
  }
};

SCENARIO("Test MappingManager::route_circuit") {
  Node node0("test_node", 0), node1("test_node", 1), node2("test_node", 2);
  Architecture arc({{node0, node1}, {node1, node2}});
  ArchitecturePtr shared_arc = std::make_shared<Architecture>(arc);
  MappingManager test_mm(shared_arc);
  std::vector<RoutingMethodPtr> test_vrm = {std::make_shared<RoutingMethod>()};
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
  GIVEN("Method that invokes a permutation from token swapping stage.") {
    Circuit circ(3);
    std::vector<Qubit> qubits = circ.all_qubits();
    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[2]});
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], node0}, {qubits[1], node1}, {qubits[2], node2}};
    circ.rename_units(rename_map);
    std::vector<RoutingMethodPtr> test_ts_rm = {
        std::make_shared<TokenSwappingTester>()};
    test_mm.route_circuit(circ, test_ts_rm);

    std::vector<Command> commands = circ.get_commands();
    REQUIRE(commands.size() == 3);
    Command c0 = commands[0];
    unit_vector_t uid_swap_12 = {node1, node2};
    REQUIRE(c0.get_args() == uid_swap_12);
    REQUIRE(*c0.get_op_ptr() == *get_op_ptr(OpType::SWAP));

    Command c1 = commands[1];
    unit_vector_t uid_swap_01 = {node0, node1};
    REQUIRE(c1.get_args() == uid_swap_01);
    REQUIRE(*c1.get_op_ptr() == *get_op_ptr(OpType::SWAP));

    Command c2 = commands[2];
    unit_vector_t uid_cx_10 = {node1, node0};
    REQUIRE(c2.get_args() == uid_cx_10);
    REQUIRE(*c2.get_op_ptr() == *get_op_ptr(OpType::CX));
  }
}
}  // namespace tket
