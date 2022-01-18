#include <catch2/catch.hpp>

#include "Mapping/MappingManager.hpp"
#include "Mapping/MultiGateReorder.hpp"
#include "Simulation/CircuitSimulator.hpp"
#include "Simulation/ComparisonFunctions.hpp"

namespace tket {
SCENARIO("Reorder circuits") {
  std::vector<Node> nodes = {
      Node("test_node", 0), Node("test_node", 1), Node("test_node", 2),
      Node("node_test", 3)};

  // n0 -- n1 -- n2 -- n3
  Architecture architecture(
      {{nodes[0], nodes[1]}, {nodes[1], nodes[2]}, {nodes[2], nodes[3]}});
  ArchitecturePtr shared_arc = std::make_shared<Architecture>(architecture);
  GIVEN("Simple CZ circuit.") {
    Circuit circ(4);
    std::vector<Qubit> qubits = circ.all_qubits();
    // Physically invalid operations
    circ.add_op<UnitID>(OpType::CZ, {qubits[0], qubits[2]});
    circ.add_op<UnitID>(OpType::CZ, {qubits[0], qubits[3]});
    // Physically valid operations
    circ.add_op<UnitID>(OpType::CZ, {qubits[0], qubits[1]});
    circ.add_op<UnitID>(OpType::CZ, {qubits[2], qubits[3]});
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]},
        {qubits[1], nodes[1]},
        {qubits[2], nodes[2]},
        {qubits[3], nodes[3]}};
    circ.rename_units(rename_map);
    Circuit circ_copy(circ);
    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);
    mf->advance_frontier_boundary(shared_arc);
    MultiGateReorder mr(shared_arc, mf);
    mr.solve(20, 20);
    std::vector<Command> commands = circ.get_commands();
    for (unsigned i = 0; i < 2; i++) {
      std::vector<Node> nodes;
      for (auto arg : commands[i].get_args()) {
        nodes.push_back(Node(arg));
      }
      REQUIRE(shared_arc->valid_operation(nodes));
    }
    const auto u = tket_sim::get_unitary(circ);
    const auto u1 = tket_sim::get_unitary(circ_copy);
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(
        u, u1, tket_sim::MatrixEquivalence::EQUAL));
  }

  GIVEN("Simple CZ circuit 2.") {
    Circuit circ(4);
    std::vector<Qubit> qubits = circ.all_qubits();
    // Physically valid operations
    circ.add_op<UnitID>(OpType::CZ, {qubits[1], qubits[0]});
    // Physically invalid operations
    circ.add_op<UnitID>(OpType::CZ, {qubits[0], qubits[2]});
    circ.add_op<UnitID>(OpType::CZ, {qubits[0], qubits[3]});
    // Physically valid operations
    circ.add_op<UnitID>(OpType::CZ, {qubits[0], qubits[1]});
    circ.add_op<UnitID>(OpType::CZ, {qubits[2], qubits[3]});
    // Physically invalid operations
    circ.add_op<UnitID>(OpType::CZ, {qubits[3], qubits[0]});
    // Physically valid operations
    circ.add_op<UnitID>(OpType::CZ, {qubits[3], qubits[2]});
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]},
        {qubits[1], nodes[1]},
        {qubits[2], nodes[2]},
        {qubits[3], nodes[3]}};
    circ.rename_units(rename_map);
    Circuit circ_copy(circ);
    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);
    mf->advance_frontier_boundary(shared_arc);
    MultiGateReorder mr(shared_arc, mf);
    mr.solve(20, 20);
    std::vector<Command> commands = circ.get_commands();
    for (unsigned i = 0; i < 4; i++) {
      std::vector<Node> nodes;
      for (auto arg : commands[i].get_args()) {
        nodes.push_back(Node(arg));
      }
      REQUIRE(shared_arc->valid_operation(nodes));
    }
    const auto u = tket_sim::get_unitary(circ);
    const auto u1 = tket_sim::get_unitary(circ_copy);
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(
        u, u1, tket_sim::MatrixEquivalence::EQUAL));
  }
  GIVEN("Simple CZ circuit with single_qs.") {
    Circuit circ(4, 1);
    std::vector<Qubit> qubits = circ.all_qubits();
    // Physically valid operations
    circ.add_op<UnitID>(OpType::CZ, {qubits[1], qubits[0]});
    // Physically invalid operations
    circ.add_op<UnitID>(OpType::CZ, {qubits[0], qubits[2]});
    circ.add_op<UnitID>(OpType::CZ, {qubits[0], qubits[3]});
    // Physically valid operations
    circ.add_op<UnitID>(OpType::Rz, 0.5, {qubits[0]});
    circ.add_op<UnitID>(OpType::Rz, 0.5, {qubits[2]});
    circ.add_op<UnitID>(OpType::Rz, 0.5, {qubits[3]});
    circ.add_op<UnitID>(OpType::CZ, {qubits[0], qubits[1]});
    circ.add_op<UnitID>(OpType::Measure, {qubits[2], Bit(0)});
    circ.add_op<UnitID>(OpType::CZ, {qubits[2], qubits[3]});
    // Physically invalid operations
    circ.add_op<UnitID>(OpType::CZ, {qubits[3], qubits[0]});
    // Physically valid operations
    circ.add_op<UnitID>(OpType::H, {qubits[3]});
    circ.add_op<UnitID>(OpType::CZ, {qubits[3], qubits[2]});
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]},
        {qubits[1], nodes[1]},
        {qubits[2], nodes[2]},
        {qubits[3], nodes[3]}};
    circ.rename_units(rename_map);
    Circuit circ_copy(circ);

    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);
    mf->advance_frontier_boundary(shared_arc);
    MultiGateReorder mr(shared_arc, mf);
    mr.solve(20, 20);
    std::vector<Command> commands = circ.get_commands();
    for (unsigned i = 0; i < 2; i++) {
      std::vector<Node> nodes;
      for (auto arg : commands[i].get_args()) {
        nodes.push_back(Node(arg));
      }
      REQUIRE(shared_arc->valid_operation(nodes));
    }
    const auto u = tket_sim::get_unitary(circ);
    const auto u1 = tket_sim::get_unitary(circ_copy);
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(
        u, u1, tket_sim::MatrixEquivalence::EQUAL));
  }

  GIVEN("Circuit with multi qubit gates.") {
    Circuit circ(4, 1);
    std::vector<Qubit> qubits = circ.all_qubits();
    // Physically valid operations
    circ.add_op<UnitID>(OpType::CZ, {qubits[1], qubits[0]});
    // Physically invalid operations
    circ.add_op<UnitID>(OpType::CZ, {qubits[0], qubits[2]});
    // Physically valid operations
    circ.add_op<UnitID>(OpType::CCX, {qubits[1], qubits[2], qubits[3]});
    circ.add_op<UnitID>(OpType::Rx, 0.5, {qubits[3]});
    circ.add_op<UnitID>(OpType::CX, {qubits[2], qubits[3]});
    circ.add_op<UnitID>(OpType::Rz, 0.5, {qubits[0]});
    circ.add_op<UnitID>(OpType::CRz, 0.5, {qubits[0], qubits[1]});
    circ.add_op<UnitID>(OpType::ZZPhase, 0.2, {qubits[0], qubits[1]});
    // Physically invalid operations
    circ.add_op<UnitID>(OpType::CZ, {qubits[3], qubits[0]});
    // Physically valid operations
    circ.add_op<UnitID>(OpType::H, {qubits[3]});
    circ.add_op<UnitID>(OpType::CZ, {qubits[3], qubits[2]});

    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]},
        {qubits[1], nodes[1]},
        {qubits[2], nodes[2]},
        {qubits[3], nodes[3]}};
    circ.rename_units(rename_map);
    Circuit circ_copy(circ);

    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);
    mf->advance_frontier_boundary(shared_arc);
    MultiGateReorder mr(shared_arc, mf);
    mr.solve(20, 20);
    std::vector<Command> commands = circ.get_commands();
    for (unsigned i = 0; i < 6; i++) {
      std::vector<Node> nodes;
      for (auto arg : commands[i].get_args()) {
        nodes.push_back(Node(arg));
      }
      REQUIRE(shared_arc->valid_operation(nodes));
    }
    const auto u = tket_sim::get_unitary(circ);
    const auto u1 = tket_sim::get_unitary(circ_copy);
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(
        u, u1, tket_sim::MatrixEquivalence::EQUAL));
  }
}

SCENARIO("Reorder circuits with limited search space") {
  std::vector<Node> nodes = {
      Node("test_node", 0), Node("test_node", 1), Node("test_node", 2),
      Node("node_test", 3)};

  // n0 -- n1 -- n2 -- n3
  Architecture architecture(
      {{nodes[0], nodes[1]}, {nodes[1], nodes[2]}, {nodes[2], nodes[3]}});
  ArchitecturePtr shared_arc = std::make_shared<Architecture>(architecture);
  GIVEN("Simple CZ circuit.") {
    Circuit circ(4);
    std::vector<Qubit> qubits = circ.all_qubits();
    // Physically invalid operations
    circ.add_op<UnitID>(OpType::CZ, {qubits[0], qubits[2]});
    circ.add_op<UnitID>(OpType::CZ, {qubits[0], qubits[3]});
    // Physically valid operations
    circ.add_op<UnitID>(OpType::CZ, {qubits[0], qubits[1]});
    circ.add_op<UnitID>(OpType::CZ, {qubits[2], qubits[1]});
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]},
        {qubits[1], nodes[1]},
        {qubits[2], nodes[2]},
        {qubits[3], nodes[3]}};
    circ.rename_units(rename_map);
    Circuit circ_copy(circ);
    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);
    mf->advance_frontier_boundary(shared_arc);
    MultiGateReorder mr(shared_arc, mf);
    mr.solve(3, 3);
    // Check only the first valid CZ get commuted to the front
    std::vector<Command> commands = circ.get_commands();
    REQUIRE(shared_arc->valid_operation(
        {Node(commands[0].get_args()[0]), Node(commands[0].get_args()[1])}));
    REQUIRE(!shared_arc->valid_operation(
        {Node(commands[1].get_args()[0]), Node(commands[1].get_args()[1])}));
    const auto u = tket_sim::get_unitary(circ);
    const auto u1 = tket_sim::get_unitary(circ_copy);
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(
        u, u1, tket_sim::MatrixEquivalence::EQUAL));
  }
}

SCENARIO("Test MultiGateReorderRoutingMethod") {
  std::vector<Node> nodes = {
      Node("test_node", 0), Node("test_node", 1), Node("test_node", 2),
      Node("node_test", 3)};

  // n0 -- n1 -- n2 -- n3
  Architecture architecture(
      {{nodes[0], nodes[1]}, {nodes[1], nodes[2]}, {nodes[2], nodes[3]}});
  ArchitecturePtr shared_arc = std::make_shared<Architecture>(architecture);
  GIVEN("Simple CZ circuit.") {
    Circuit circ(4);
    std::vector<Qubit> qubits = circ.all_qubits();
    // Physically valid operations
    circ.add_op<UnitID>(OpType::CZ, {qubits[0], qubits[1]});
    circ.add_op<UnitID>(OpType::CZ, {qubits[2], qubits[3]});
    // Physically invalid operations
    circ.add_op<UnitID>(OpType::CZ, {qubits[0], qubits[2]});
    circ.add_op<UnitID>(OpType::CZ, {qubits[0], qubits[3]});
    // Physically valid operations
    circ.add_op<UnitID>(OpType::CZ, {qubits[0], qubits[1]});
    circ.add_op<UnitID>(OpType::CZ, {qubits[2], qubits[3]});
    circ.add_op<UnitID>(OpType::CZ, {qubits[2], qubits[3]});
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]},
        {qubits[1], nodes[1]},
        {qubits[2], nodes[2]},
        {qubits[3], nodes[3]}};
    circ.rename_units(rename_map);
    Circuit circ_copy(circ);
    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);
    mf->advance_frontier_boundary(shared_arc);
    MultiGateReorderRoutingMethod mrrm;
    REQUIRE(mrrm.check_method(mf, shared_arc));

    unit_map_t init_map = mrrm.routing_method(mf, shared_arc);
    REQUIRE(init_map.size() == 0);
    std::vector<Command> commands = circ.get_commands();
    for (unsigned i = 0; i < 5; i++) {
      std::vector<Node> nodes;
      for (auto arg : commands[i].get_args()) {
        nodes.push_back(Node(arg));
      }
      REQUIRE(shared_arc->valid_operation(nodes));
    }
    const auto u = tket_sim::get_unitary(circ);
    const auto u1 = tket_sim::get_unitary(circ_copy);
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(
        u, u1, tket_sim::MatrixEquivalence::EQUAL));

    // Test with limits
    Circuit circ2(circ_copy);

    std::shared_ptr<MappingFrontier> mf2 =
        std::make_shared<MappingFrontier>(circ2);
    mf2->advance_frontier_boundary(shared_arc);
    MultiGateReorderRoutingMethod mrrm2(4, 4);
    REQUIRE(mrrm2.check_method(mf2, shared_arc));

    unit_map_t init_map2 = mrrm2.routing_method(mf2, shared_arc);
    REQUIRE(init_map2.size() == 0);
    std::vector<Command> commands2 = circ2.get_commands();
    for (unsigned i = 0; i < 4; i++) {
      std::vector<Node> nodes;
      for (auto arg : commands2[i].get_args()) {
        nodes.push_back(Node(arg));
      }
      REQUIRE(shared_arc->valid_operation(nodes));
    }
    std::vector<Node> nodes;
    for (auto arg : commands2[4].get_args()) {
      nodes.push_back(Node(arg));
    }
    REQUIRE(!shared_arc->valid_operation(nodes));
    const auto u2 = tket_sim::get_unitary(circ2);
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(
        u2, u1, tket_sim::MatrixEquivalence::EQUAL));
  }
}
}  // namespace tket