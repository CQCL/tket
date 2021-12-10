#include <catch2/catch.hpp>
#include "Transformations/Transform.hpp"

namespace tket {
    SCENARIO("Test Transform::reorder_cz") {
    std::vector<Node> nodes = {Node("test_node", 0), Node("test_node", 1),
                             Node("test_node", 2), Node("node_test", 3)};
    
    // n0 -- n1 -- n2 -- n3
    Architecture architecture(
      {{nodes[0], nodes[1]},
       {nodes[1], nodes[2]},
       {nodes[2], nodes[3]}});
    ArchitecturePtr shared_arc = std::make_shared<Architecture>(architecture);
    Circuit circ(4);
    std::vector<Qubit> qubits = circ.all_qubits();
    // Physically invalid operations
    circ.add_op<UnitID>(OpType::CZ, {qubits[0], qubits[2]});
    circ.add_op<UnitID>(OpType::CZ, {qubits[0], qubits[3]});
    // Physically valid operations
    circ.add_op<UnitID>(OpType::CZ, {qubits[0], qubits[1]});
    circ.add_op<UnitID>(OpType::CZ, {qubits[2], qubits[3]});
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]}, {qubits[1], nodes[1]}, {qubits[2], nodes[2]},
        {qubits[3], nodes[3]}};
    circ.rename_units(rename_map);
    std::cout<<circ;
    Transform::reorder_cz(shared_arc).apply(circ);
    std::cout<<circ;
    std::vector<Command> commands = circ.get_commands();
    Command c0 = commands[0];
    Command c1 = commands[1];
    REQUIRE(shared_arc->valid_operation({Node(c0.get_args()[0]), Node(c0.get_args()[1])}));
    REQUIRE(shared_arc->valid_operation({Node(c1.get_args()[0]), Node(c1.get_args()[1])}));
    
}
}