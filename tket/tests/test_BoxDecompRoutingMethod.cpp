// Copyright 2019-2022 Cambridge Quantum Computing
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
#include <catch2/catch_test_macros.hpp>

#include "Mapping/BoxDecomposition.hpp"
#include "Mapping/LexiRoute.hpp"
#include "Mapping/MappingManager.hpp"
#include "Predicates/Predicates.hpp"
#include "Simulation/CircuitSimulator.hpp"
#include "Simulation/ComparisonFunctions.hpp"

namespace tket {
SCENARIO("Decompose boxes") {
  std::vector<Node> nodes = {
      Node("test_node", 0), Node("test_node", 1), Node("test_node", 2),
      Node("node_test", 3)};

  // n0 -- n1 -- n2 -- n3
  Architecture architecture(
      {{nodes[0], nodes[1]}, {nodes[1], nodes[2]}, {nodes[2], nodes[3]}});
  ArchitecturePtr shared_arc = std::make_shared<Architecture>(architecture);

  Eigen::Matrix4cd m;
  m << 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0;
  Unitary2qBox ubox(m);

  GIVEN("A box") {
    Circuit circ(4);
    std::vector<Qubit> qubits = circ.all_qubits();

    circ.add_box(ubox, {0, 2});
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]},
        {qubits[1], nodes[1]},
        {qubits[2], nodes[2]},
        {qubits[3], nodes[3]}};
    circ.rename_units(rename_map);
    Circuit circ_copy(circ);
    MappingFrontier_ptr mf = std::make_shared<MappingFrontier>(circ);
    BoxDecomposition bd(shared_arc, mf);
    bd.solve();
    const auto u = tket_sim::get_unitary(circ);
    const auto u1 = tket_sim::get_unitary(circ_copy);
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(
        u, u1, tket_sim::MatrixEquivalence::EQUAL));
    std::vector<Command> commands = mf->circuit_.get_commands();
    for (Command c : commands) {
      REQUIRE(!c.get_op_ptr()->get_desc().is_box());
    }
  }

  GIVEN("A conditional box") {
    Circuit circ(4, 1);
    std::vector<Qubit> qubits = circ.all_qubits();
    Conditional cond(std::make_shared<Unitary2qBox>(ubox), 1, 1);
    circ.add_op<UnitID>(
        std::make_shared<Conditional>(cond), {Bit(0), Qubit(0), Qubit(1)});
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]},
        {qubits[1], nodes[1]},
        {qubits[2], nodes[2]},
        {qubits[3], nodes[3]}};
    circ.rename_units(rename_map);
    MappingFrontier_ptr mf = std::make_shared<MappingFrontier>(circ);
    BoxDecomposition bd(shared_arc, mf);
    bd.solve();
    std::vector<Command> commands = mf->circuit_.get_commands();
    for (Command c : commands) {
      Op_ptr op = c.get_op_ptr();
      REQUIRE(
          !(op->get_desc().is_box() || (op->get_type() == OpType::Conditional &&
                                        static_cast<const Conditional &>(*op)
                                            .get_op()
                                            ->get_desc()
                                            .is_box())));
    }
  }

  GIVEN("Test BoxDecompositionRoutingMethod") {
    Circuit circ(4, 1);
    std::vector<Qubit> qubits = circ.all_qubits();
    circ.add_box(ubox, {0, 3});
    circ.add_op<UnitID>(OpType::CZ, {qubits[0], qubits[1]});
    circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[3]});
    circ.add_box(ubox, {1, 3});
    circ.add_box(ubox, {0, 1});
    circ.add_op<UnitID>(OpType::X, {qubits[1]});
    circ.add_op<unsigned>(OpType::Measure, {0, 0});
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]},
        {qubits[1], nodes[1]},
        {qubits[2], nodes[2]},
        {qubits[3], nodes[3]}};
    circ.rename_units(rename_map);
    MappingFrontier_ptr mf = std::make_shared<MappingFrontier>(circ);
    MappingManager mm(shared_arc);
    std::vector<RoutingMethodPtr> vrm = {

        std::make_shared<LexiRouteRoutingMethod>(10),
        std::make_shared<BoxDecompositionRoutingMethod>()};
    bool res = mm.route_circuit(circ, vrm);
    REQUIRE(res);
    PredicatePtr routed_correctly =
        std::make_shared<ConnectivityPredicate>(architecture);
    REQUIRE(routed_correctly->verify(circ));
    std::vector<Command> commands = mf->circuit_.get_commands();
    for (Command c : commands) {
      REQUIRE(!c.get_op_ptr()->get_desc().is_box());
    }
  }
}

SCENARIO("Test JSON serialisation for BoxDecompositionRoutingMethod") {
  GIVEN("BoxDecompositionRoutingMethod") {
    nlohmann::json j_rm;
    j_rm["name"] = "BoxDecompositionRoutingMethod";
    BoxDecompositionRoutingMethod rm_loaded =
        BoxDecompositionRoutingMethod::deserialize(j_rm);
    nlohmann::json j_rm_serialised = rm_loaded.serialize();
    REQUIRE(j_rm == j_rm_serialised);
  }

  GIVEN("BoxDecompositionRoutingMethod vector") {
    nlohmann::json j_rms = {
        {{"name", "BoxDecompositionRoutingMethod"}},
        {
            {"name", "LexiRouteRoutingMethod"},
            {"depth", 3},
        }};
    std::vector<RoutingMethodPtr> rms =
        j_rms.get<std::vector<RoutingMethodPtr>>();
    nlohmann::json j_rms_serialised = rms;
    REQUIRE(j_rms == j_rms_serialised);
  }
}

}  // namespace tket
