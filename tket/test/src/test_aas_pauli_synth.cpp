#include <catch2/catch_test_macros.hpp>
#include <tket/Transformations/Decomposition.hpp>

#include "CircuitsForTesting.hpp"
#include "testutil.hpp"
#include "tket/ArchAwareSynth/SteinerForest.hpp"
#include "tket/Architecture/Architecture.hpp"
#include "tket/Circuit/Boxes.hpp"
#include "tket/Circuit/PauliExpBoxes.hpp"
#include "tket/Converters/Converters.hpp"
#include "tket/Converters/PauliGadget.hpp"
#include "tket/Diagonalisation/Diagonalisation.hpp"
#include "tket/Gate/SymTable.hpp"
#include "tket/PauliGraph/ConjugatePauliFunctions.hpp"
#include "tket/PauliGraph/PauliGraph.hpp"
#include "tket/Simulation/CircuitSimulator.hpp"
#include "tket/Transformations/PauliOptimisation.hpp"
#include "tket/Transformations/Rebase.hpp"

namespace tket {
namespace test_aas_pauli_synth {

static void add_ops_to_prepend_1(Circuit& circ) {
  circ.add_op<unsigned>(OpType::Rx, 1.511, {2});
  circ.add_op<unsigned>(OpType::Rz, 0.745, {2});
}
static void add_ops_to_prepend_2(Circuit& circ) {
  circ.add_op<unsigned>(OpType::Rx, 0.849, {3});
  circ.add_op<unsigned>(OpType::Rz, 0.102, {3});
}

SCENARIO("Test AAS pauli synth") {
  GIVEN("4 qb 3 Pauli Gadget circuit") {
    auto prepend = CircuitsForTesting::get_prepend_circuit(3);
    add_ops_to_prepend_1(prepend);
    Circuit circ(4);
    PauliExpBox peb({Pauli::Z, Pauli::Z, Pauli::Z, Pauli::Z}, 0.333);
    circ.add_box(peb, {0, 1, 2, 3});
    PauliExpBox peb2({Pauli::X, Pauli::Z, Pauli::X, Pauli::I}, 0.233);
    circ.add_box(peb2, {0, 1, 2, 3});
    PauliExpBox peb3({Pauli::X, Pauli::X, Pauli::X, Pauli::X}, 0.174);
    circ.add_box(peb3, {0, 1, 2, 3});
    Circuit test_circ = prepend >> circ;
    // 2. Defined a grid architecture
    std::vector<Node> nodes = {
        Node("a", 0), Node("a", 1), Node("a", 2), Node("a", 3)};

    Architecture arch(
        {{nodes[0], nodes[1]},
         {nodes[2], nodes[3]},
         {nodes[0], nodes[3]},
         {nodes[1], nodes[3]}});
    // 3. manually place the qubits
    std::vector<Qubit> qubits = test_circ.all_qubits();
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]},
        {qubits[1], nodes[1]},
        {qubits[2], nodes[2]},
        {qubits[3], nodes[3]}};
    test_circ.rename_units(rename_map);
    // 4. synthesis
    PauliGraph pg = circuit_to_pauli_graph(test_circ);
    Circuit out_circ = pauli_graph_to_circuit_lazy_aas(pg, arch);
    // 5. check correctness
    REQUIRE(test_statevector_comparison(test_circ, out_circ, true));
    for (const auto& cmd : out_circ) {
      std::vector<Node> nodes;
      for (auto arg : cmd.get_args()) nodes.push_back(Node(arg));
      REQUIRE(arch.valid_operation(nodes));
    }
  }
  GIVEN("5 qubit 7 Pauli Gadget circuit") {
    // 1. prepare the test_circ
    auto prepend = CircuitsForTesting::get_prepend_circuit(5);
    add_ops_to_prepend_1(prepend);
    add_ops_to_prepend_2(prepend);
    prepend.add_op<unsigned>(OpType::Rx, 0.466, {4});
    prepend.add_op<unsigned>(OpType::Rz, 1.303, {4});
    Circuit circ(5);
    PauliExpBox peb0(
        {Pauli::I, Pauli::X, Pauli::Z, Pauli::I, Pauli::Z}, 0.3112);
    circ.add_box(peb0, {0, 1, 2, 3, 4});
    PauliExpBox peb1({Pauli::I, Pauli::Y, Pauli::I, Pauli::Z, Pauli::Y}, 1.178);
    circ.add_box(peb1, {0, 1, 2, 3, 4});
    PauliExpBox peb2(
        {Pauli::X, Pauli::X, Pauli::I, Pauli::Y, Pauli::I}, -0.911);
    circ.add_box(peb2, {0, 1, 2, 3, 4});
    PauliExpBox peb3(
        {Pauli::Y, Pauli::Y, Pauli::X, Pauli::I, Pauli::I}, 0.7122);
    circ.add_box(peb3, {0, 1, 2, 3, 4});
    PauliExpBox peb4({Pauli::Z, Pauli::I, Pauli::Y, Pauli::X, Pauli::X}, 1.102);
    circ.add_box(peb4, {0, 1, 2, 3, 4});
    PauliExpBox peb5({Pauli::Z, Pauli::X, Pauli::I, Pauli::Z, Pauli::Z}, 0.151);
    circ.add_box(peb5, {0, 1, 2, 3, 4});
    PauliExpBox peb6({Pauli::Z, Pauli::Y, Pauli::Z, Pauli::I, Pauli::Y}, 1.223);
    circ.add_box(peb6, {0, 1, 2, 3, 4});
    Circuit test_circ = prepend >> circ;

    // 2. Defined a line architecture
    std::vector<Node> nodes = {
        Node("a", 0), Node("a", 1), Node("a", 2), Node("a", 3), Node("a", 4)};

    Architecture arch(
        {{nodes[0], nodes[1]},
         {nodes[1], nodes[2]},
         {nodes[2], nodes[3]},
         {nodes[3], nodes[4]}});
    // 3. manually place the qubits
    std::vector<Qubit> qubits = test_circ.all_qubits();
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]},
        {qubits[1], nodes[1]},
        {qubits[2], nodes[2]},
        {qubits[3], nodes[3]},
        {qubits[4], nodes[4]}};
    test_circ.rename_units(rename_map);
    // 4. synthesis
    PauliGraph pg = circuit_to_pauli_graph(test_circ);
    Circuit out_circ = pauli_graph_to_circuit_lazy_aas(pg, arch);
    // 5. check correctness
    REQUIRE(test_statevector_comparison(test_circ, out_circ, true));
    for (const auto& cmd : out_circ) {
      std::vector<Node> nodes;
      for (auto arg : cmd.get_args()) nodes.push_back(Node(arg));
      REQUIRE(arch.valid_operation(nodes));
    }
  }
}

}  // namespace test_aas_pauli_synth
}  // namespace tket
