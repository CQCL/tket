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

#include <catch2/catch.hpp>

#include "Mapping/LexiRoute.hpp"
#include "Mapping/LexiLabelling.hpp"
#include "Mapping/MappingManager.hpp"
#include "Mapping/Verification.hpp"
#include "Predicates/CompilationUnit.hpp"
#include "Predicates/CompilerPass.hpp"
#include "Predicates/PassGenerators.hpp"
#include "Predicates/PassLibrary.hpp"
#include "Transformations/Decomposition.hpp"
#include "testutil.hpp"
namespace tket {
SCENARIO("Test LexiRoute::solve") {
  std::vector<Node> nodes = {Node("test_node", 0), Node("test_node", 1),
                             Node("test_node", 2), Node("node_test", 3),
                             Node("node_test", 4), Node("node_test", 5),
                             Node("test_node", 6), Node("node_test", 7)};
  // n0 -- n1 -- n2 -- n3 -- n4
  //             |     |
  //             n5    n7
  //             |
  //             n6
  Architecture architecture(
      {{nodes[0], nodes[1]},
       {nodes[1], nodes[2]},
       {nodes[2], nodes[3]},
       {nodes[3], nodes[4]},
       {nodes[2], nodes[5]},
       {nodes[5], nodes[6]},
       {nodes[3], nodes[7]}});
  ArchitecturePtr shared_arc = std::make_shared<Architecture>(architecture);
  GIVEN("Single best solution, all qubits labelled.") {
    Circuit circ(6);
    std::vector<Qubit> qubits = circ.all_qubits();
    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[2]});
    circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[3]});
    circ.add_op<UnitID>(OpType::CX, {qubits[4], qubits[5]});
    // n0 -- n1 -- n2 -- n3 -- n4
    //             |     |
    //             n5    n7
    //             |
    //             n6
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]}, {qubits[1], nodes[1]}, {qubits[2], nodes[2]},
        {qubits[3], nodes[3]}, {qubits[4], nodes[6]}, {qubits[5], nodes[5]}};
    circ.rename_units(rename_map);
    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);
    LexiRoute lr(shared_arc, mf);

    lr.solve(4);
    std::vector<Command> commands = mf->circuit_.get_commands();
    REQUIRE(commands.size() == 4);
    Command swap_c = commands[1];
    unit_vector_t uids = {nodes[1], nodes[2]};
    REQUIRE(swap_c.get_args() == uids);
    REQUIRE(*swap_c.get_op_ptr() == *get_op_ptr(OpType::SWAP));
  }
  GIVEN("Single best solution, one qubit unlabelled.") {
    Circuit circ(6);
    std::vector<Qubit> qubits = circ.all_qubits();
    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[2]});
    circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[3]});
    circ.add_op<UnitID>(OpType::CX, {qubits[4], qubits[5]});
    // n0 -- n1 -- n2 -- n3 -- n4
    //             |     |
    //             n5    n7
    //             |
    //             n6
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]},
        {qubits[1], nodes[1]},
        {qubits[2], nodes[2]},
        {qubits[3], nodes[3]},
        {qubits[5], nodes[5]}};
    circ.rename_units(rename_map);
    std::shared_ptr<MappingFrontier> mf0 =
        std::make_shared<MappingFrontier>(circ);
    LexiRoute lr(shared_arc, mf0);

    lr.solve(4);

    REQUIRE(mf0->circuit_.n_gates() == 3);

    rename_map = {{qubits[4], nodes[6]}};
    mf0->circuit_.rename_units(rename_map);

    std::shared_ptr<MappingFrontier> mf1 =
        std::make_shared<MappingFrontier>(circ);
    LexiRoute lr1(shared_arc, mf1);
    lr1.solve(4);
    std::vector<Command> commands = mf1->circuit_.get_commands();
    Command swap_c = commands[1];
    unit_vector_t uids = {nodes[1], nodes[2]};
    REQUIRE(swap_c.get_args() == uids);
    REQUIRE(*swap_c.get_op_ptr() == *get_op_ptr(OpType::SWAP));
  }
  GIVEN("Single best solution, one stage of look-ahead required.") {
    Circuit circ(8);
    std::vector<Qubit> qubits = circ.all_qubits();
    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[4]});
    circ.add_op<UnitID>(OpType::CX, {qubits[6], qubits[7]});
    circ.add_op<UnitID>(OpType::CX, {qubits[2], qubits[7]});
    //                   n7
    //                   |
    // n0 -- n1 -- n2 -- n3 -- n4
    //             |
    //             n5
    //             |
    //             n6
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]}, {qubits[1], nodes[1]}, {qubits[2], nodes[2]},
        {qubits[3], nodes[3]}, {qubits[4], nodes[4]}, {qubits[5], nodes[5]},
        {qubits[6], nodes[6]}, {qubits[7], nodes[7]}};
    circ.rename_units(rename_map);
    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);
    LexiRoute lr(shared_arc, mf);

    lr.solve(4);
    std::vector<Command> commands = mf->circuit_.get_commands();
    REQUIRE(commands.size() == 4);
    Command swap_c = commands[0];
    unit_vector_t uids = {nodes[7], nodes[3]};
    REQUIRE(swap_c.get_args() == uids);
    REQUIRE(*swap_c.get_op_ptr() == *get_op_ptr(OpType::SWAP));

    Command changed_c = commands[3];
    uids = {nodes[2], nodes[3]};
    REQUIRE(changed_c.get_args() == uids);
  }
  GIVEN("All unlabelled, labelling can give complete solution.") {
    Circuit circ(5);
    std::vector<Qubit> qubits = circ.all_qubits();
    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[1]});
    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[2]});
    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[3]});
    circ.add_op<UnitID>(OpType::CX, {qubits[3], qubits[4]});

    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);
    LexiRoute lr0(shared_arc, mf);
    lr0.solve(20);
    std::vector<Command> commands = mf->circuit_.get_commands();
    REQUIRE(commands.size() == 4);
    Command c = commands[0];
    unit_vector_t uids = {nodes[2], nodes[1]};
    REQUIRE(c.get_args() == uids);
    mf->advance_frontier_boundary(shared_arc);

    LexiRoute lr1(shared_arc, mf);
    lr1.solve(20);
    uids = {nodes[2], nodes[3]};
    REQUIRE(mf->circuit_.get_commands()[1].get_args() == uids);
    mf->advance_frontier_boundary(shared_arc);

    LexiRoute lr2(shared_arc, mf);
    lr2.solve(20);
    uids = {nodes[2], nodes[5]};
    REQUIRE(mf->circuit_.get_commands()[2].get_args() == uids);
    mf->advance_frontier_boundary(shared_arc);

    LexiRoute lr3(shared_arc, mf);
    lr3.solve(20);
    uids = {nodes[5], nodes[6]};
    REQUIRE(mf->circuit_.get_commands()[3].get_args() == uids);
  }
  GIVEN("Bridge preferred, CX.") {
    Circuit circ(5);
    std::vector<Qubit> qubits = circ.all_qubits();
    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[1]});
    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[2]});
    circ.add_op<UnitID>(OpType::CX, {qubits[3], qubits[1]});
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[1]},
        {qubits[1], nodes[3]},
        {qubits[2], nodes[0]},
        {qubits[3], nodes[7]},
        {qubits[4], nodes[2]}};
    circ.rename_units(rename_map);
    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);

    mf->advance_frontier_boundary(shared_arc);
    LexiRoute lr(shared_arc, mf);
    lr.solve(4);
    Command bridge_c = mf->circuit_.get_commands()[0];
    unit_vector_t uids = {nodes[1], nodes[2], nodes[3]};
    REQUIRE(bridge_c.get_args() == uids);
    REQUIRE(*bridge_c.get_op_ptr() == *get_op_ptr(OpType::BRIDGE));
  }
  GIVEN("Bridge preferred, CZ.") {
    Circuit circ(5);
    std::vector<Qubit> qubits = circ.all_qubits();
    circ.add_op<UnitID>(OpType::CZ, {qubits[0], qubits[1]});
    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[2]});
    circ.add_op<UnitID>(OpType::CX, {qubits[3], qubits[1]});
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[1]},
        {qubits[1], nodes[3]},
        {qubits[2], nodes[0]},
        {qubits[3], nodes[7]},
        {qubits[4], nodes[2]}};
    circ.rename_units(rename_map);
    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);

    mf->advance_frontier_boundary(shared_arc);
    LexiRoute lr(shared_arc, mf);
    lr.solve(4);
    REQUIRE(mf->circuit_.get_commands().size() == 4);
  }
  GIVEN("Bridge preferred, conditional CX.") {
    Circuit circ(5, 1);
    std::vector<Qubit> qubits = circ.all_qubits();
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {0}, 1);
    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[2]});
    circ.add_op<UnitID>(OpType::CX, {qubits[3], qubits[1]});
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[1]},
        {qubits[1], nodes[3]},
        {qubits[2], nodes[0]},
        {qubits[3], nodes[7]},
        {qubits[4], nodes[2]}};
    circ.rename_units(rename_map);
    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);

    mf->advance_frontier_boundary(shared_arc);
    LexiRoute lr(shared_arc, mf);
    lr.solve(4);
    REQUIRE(mf->circuit_.get_commands().size() == 4);
  }
  GIVEN("Bridge preferred, conditional CZ.") {
    Circuit circ(5, 1);
    std::vector<Qubit> qubits = circ.all_qubits();
    circ.add_conditional_gate<unsigned>(OpType::CZ, {}, {0, 1}, {0}, 1);
    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[2]});
    circ.add_op<UnitID>(OpType::CX, {qubits[3], qubits[1]});
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[1]},
        {qubits[1], nodes[3]},
        {qubits[2], nodes[0]},
        {qubits[3], nodes[7]},
        {qubits[4], nodes[2]}};
    circ.rename_units(rename_map);
    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);

    mf->advance_frontier_boundary(shared_arc);
    LexiRoute lr(shared_arc, mf);
    lr.solve(4);
    REQUIRE(mf->circuit_.get_commands().size() == 4);
  }

  GIVEN("Ancilla assignment, one valid node.") {
    Circuit circ(3);
    std::vector<Qubit> qubits = circ.all_qubits();
    circ.add_op<UnitID>(OpType::CZ, {qubits[0], qubits[1]});
    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[2]});

    std::vector<Node> nodes = {
        Node("test_node", 0), Node("test_node", 1), Node("test_node", 2),
        Node("node_test", 3), Node("node_test", 4)};
    // just a ring

    Architecture architecture(
        {{nodes[0], nodes[1]},
         {nodes[1], nodes[2]},
         {nodes[2], nodes[3]},
         {nodes[3], nodes[4]},
         {nodes[4], nodes[0]}});
    ArchitecturePtr shared_arc = std::make_shared<Architecture>(architecture);

    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[2]}, {qubits[1], nodes[4]}};
    circ.rename_units(rename_map);

    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);
    mf->advance_frontier_boundary(shared_arc);
    LexiRoute lr0(shared_arc, mf);
    lr0.solve(20);
    REQUIRE(circ.all_qubits()[1] == nodes[4]);

    mf->advance_frontier_boundary(shared_arc);
    LexiRoute lr1(shared_arc, mf);
    lr1.solve(20);
    REQUIRE(circ.all_qubits()[0] == nodes[3]);
  }

  GIVEN("Ancilla assignment, multiple valid Node.") {
    Circuit circ(3);
    std::vector<Qubit> qubits = circ.all_qubits();
    circ.add_op<UnitID>(OpType::CZ, {qubits[0], qubits[1]});
    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[2]});

    std::vector<Node> nodes = {Node("test_node", 0), Node("test_node", 1),
                               Node("test_node", 2), Node("node_test", 3),
                               Node("node_test", 4), Node("node_test", 5),
                               Node("node_test", 6)};
    // A ring, but with two identical length paths where ancilla could be
    // assigned
    Architecture architecture(
        {{nodes[0], nodes[1]},
         {nodes[1], nodes[2]},
         {nodes[2], nodes[3]},
         {nodes[2], nodes[5]},
         {nodes[3], nodes[6]},
         {nodes[5], nodes[6]},
         {nodes[3], nodes[4]},
         {nodes[5], nodes[4]},
         {nodes[4], nodes[0]}});
    ArchitecturePtr shared_arc = std::make_shared<Architecture>(architecture);

    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[2]}, {qubits[1], nodes[4]}};
    circ.rename_units(rename_map);

    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);
    mf->advance_frontier_boundary(shared_arc);
    LexiRoute lr0(shared_arc, mf);
    lr0.solve(20);

    mf->advance_frontier_boundary(shared_arc);
    LexiRoute lr1(shared_arc, mf);
    lr1.solve(20);

    REQUIRE(circ.all_qubits()[1] == nodes[5]);
  }
  GIVEN("Ancilla assignment, one valid Node, with merge.") {
    Circuit circ(4);
    std::vector<Qubit> qubits = circ.all_qubits();
    circ.add_op<UnitID>(OpType::CZ, {qubits[0], qubits[1]});
    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[2]});
    circ.add_op<UnitID>(OpType::H, {qubits[3]});

    std::vector<Node> nodes = {
        Node("test_node", 0), Node("test_node", 1), Node("test_node", 2),
        Node("node_test", 3), Node("node_test", 4)};
    // just a ring

    Architecture architecture(
        {{nodes[0], nodes[1]},
         {nodes[1], nodes[2]},
         {nodes[2], nodes[3]},
         {nodes[3], nodes[4]},
         {nodes[4], nodes[0]}});
    ArchitecturePtr shared_arc = std::make_shared<Architecture>(architecture);

    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[2]}, {qubits[1], nodes[4]}, {qubits[3], nodes[3]}};
    circ.rename_units(rename_map);

    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);
    mf->ancilla_nodes_.insert(nodes[3]);
    mf->advance_frontier_boundary(shared_arc);
    LexiRoute lr0(shared_arc, mf);
    lr0.solve(20);

    REQUIRE(circ.all_qubits()[1] == nodes[4]);
    REQUIRE(circ.all_qubits()[0] == nodes[3]);
  }
  GIVEN(
      "Single best solution, with measurements and classically controlled "
      "gates.") {
    Circuit circ(6, 1);
    std::vector<Qubit> qubits = circ.all_qubits();
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 2}, {0}, 1);
    circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[3]});
    circ.add_conditional_gate<unsigned>(OpType::X, {}, {0}, {0}, 1);
    circ.add_op<UnitID>(OpType::Measure, {qubits[1], Bit(0)});
    circ.add_op<UnitID>(OpType::CX, {qubits[4], qubits[5]});
    circ.add_op<UnitID>(OpType::Measure, {qubits[3], Bit(0)});
    // n0 -- n1 -- n2 -- n3 -- n4
    //             |     |
    //             n5    n7
    //             |
    //             n6
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]}, {qubits[1], nodes[1]}, {qubits[2], nodes[2]},
        {qubits[3], nodes[3]}, {qubits[4], nodes[6]}, {qubits[5], nodes[5]}};
    circ.rename_units(rename_map);
    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);
    LexiRoute lr(shared_arc, mf);

    lr.solve(4);
    std::vector<Command> commands = mf->circuit_.get_commands();
    REQUIRE(commands.size() == 7);
    Command swap_c = commands[1];
    unit_vector_t uids = {nodes[1], nodes[2]};
    REQUIRE(swap_c.get_args() == uids);
    REQUIRE(*swap_c.get_op_ptr() == *get_op_ptr(OpType::SWAP));
  }
  GIVEN(
      "Labelling is required, but there are no free remaining qubits, for one "
      "updated label, order 0.") {
    Circuit circ(9);
    std::vector<Qubit> qubits = circ.all_qubits();
    circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[8]});
    // n0 -- n1 -- n2 -- n3 -- n4
    //             |     |
    //             n5    n7
    //             |
    //             n6
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]}, {qubits[1], nodes[1]}, {qubits[2], nodes[2]},
        {qubits[3], nodes[3]}, {qubits[4], nodes[4]}, {qubits[5], nodes[5]},
        {qubits[6], nodes[6]}, {qubits[7], nodes[7]}};
    circ.rename_units(rename_map);
    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);
    LexiRoute lr(shared_arc, mf);
    REQUIRE_THROWS_AS(lr.solve(1), LexiRouteError);
  }
  GIVEN(
      "Labelling is required, but there are no free remaining qubits, for one "
      "updated label, order 1.") {
    Circuit circ(9);
    std::vector<Qubit> qubits = circ.all_qubits();
    circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[8]});
    // n0 -- n1 -- n2 -- n3 -- n4
    //             |     |
    //             n5    n7
    //             |
    //             n6
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]}, {qubits[8], nodes[1]}, {qubits[2], nodes[2]},
        {qubits[3], nodes[3]}, {qubits[4], nodes[4]}, {qubits[5], nodes[5]},
        {qubits[6], nodes[6]}, {qubits[7], nodes[7]}};
    circ.rename_units(rename_map);
    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);
    LexiRoute lr(shared_arc, mf);
    REQUIRE_THROWS_AS(lr.solve(1), LexiRouteError);
  }
  GIVEN(
      "Labelling is required, but there are no free remaining qubits, for two "
      "updated labels.") {
    Circuit circ(10);
    std::vector<Qubit> qubits = circ.all_qubits();
    circ.add_op<UnitID>(OpType::CX, {qubits[9], qubits[8]});
    // n0 -- n1 -- n2 -- n3 -- n4
    //             |     |
    //             n5    n7
    //             |
    //             n6
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]}, {qubits[1], nodes[1]}, {qubits[2], nodes[2]},
        {qubits[3], nodes[3]}, {qubits[4], nodes[4]}, {qubits[5], nodes[5]},
        {qubits[6], nodes[6]}, {qubits[7], nodes[7]}};
    circ.rename_units(rename_map);
    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);
    LexiRoute lr(shared_arc, mf);
    REQUIRE_THROWS_AS(lr.solve(1), LexiRouteError);
  }
}

SCENARIO("Test LexiRouteRoutingMethod") {
  std::vector<Node> nodes = {
      Node("test_node", 0), Node("test_node", 1), Node("test_node", 2),
      Node("node_test", 3), Node("node_test", 4), Node("node_test", 5),
      Node("test_node", 6), Node("node_test", 7), Node("node_test", 8),
      Node("node_test", 9), Node("node_test", 10)};
  //       n9 -- n8 -- n10
  //             |     |
  // n0 -- n1 -- n2 -- n3 -- n4
  //             |     |
  //             n5    n7
  //             |
  //             n6
  Architecture architecture(
      {{nodes[0], nodes[1]},
       {nodes[1], nodes[2]},
       {nodes[2], nodes[3]},
       {nodes[3], nodes[4]},
       {nodes[2], nodes[5]},
       {nodes[5], nodes[6]},
       {nodes[3], nodes[7]},
       {nodes[2], nodes[8]},
       {nodes[8], nodes[9]},
       {nodes[8], nodes[10]},
       {nodes[3], nodes[10]}});
  ArchitecturePtr shared_arc = std::make_shared<Architecture>(architecture);
  GIVEN("Circuit with all qubits, labelled, stage 0.") {
    Circuit circ(11);
    std::vector<Qubit> qubits = circ.all_qubits();
    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[4]});
    circ.add_op<UnitID>(OpType::CX, {qubits[6], qubits[7]});
    circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[10]});
    circ.add_op<UnitID>(OpType::CX, {qubits[8], qubits[5]});
    circ.add_op<UnitID>(OpType::CX, {qubits[3], qubits[9]});

    circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[5]});
    circ.add_op<UnitID>(OpType::CX, {qubits[3], qubits[9]});
    circ.add_op<UnitID>(OpType::CX, {qubits[10], qubits[0]});
    circ.add_op<UnitID>(OpType::CX, {qubits[6], qubits[0]});

    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]}, {qubits[1], nodes[1]},  {qubits[2], nodes[2]},
        {qubits[3], nodes[3]}, {qubits[4], nodes[4]},  {qubits[5], nodes[5]},
        {qubits[6], nodes[6]}, {qubits[7], nodes[7]},  {qubits[8], nodes[8]},
        {qubits[9], nodes[9]}, {qubits[10], nodes[10]}};
    circ.rename_units(rename_map);

    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);
    LexiRouteRoutingMethod lrrm(100);
    REQUIRE(lrrm.check_method(mf, shared_arc));

    unit_map_t init_map = lrrm.routing_method(mf, shared_arc);
    REQUIRE(init_map.size() == 0);

    std::vector<Command> commands = mf->circuit_.get_commands();
    REQUIRE(commands.size() == 9);
    Command bridge_c = commands[2];
    unit_vector_t uids = {nodes[8], nodes[2], nodes[5]};
    REQUIRE(bridge_c.get_args() == uids);
    REQUIRE(*bridge_c.get_op_ptr() == *get_op_ptr(OpType::BRIDGE));
  }
  GIVEN("Circuit with all qubits, labelled, stage 1.") {
    Circuit circ(11);
    std::vector<Qubit> qubits = circ.all_qubits();
    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[4]});
    circ.add_op<UnitID>(OpType::CX, {qubits[6], qubits[7]});
    circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[10]});
    circ.add_op<UnitID>(OpType::CX, {qubits[8], qubits[5]});
    circ.add_op<UnitID>(OpType::CX, {qubits[3], qubits[9]});
    //       n9 -- n8 -- n3
    //             |     |
    // n0 -- n1 -- n2 -- n10 -- n4
    //             |     |
    //             n6    n7
    //             |
    //             n5
    circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[5]});
    circ.add_op<UnitID>(OpType::CX, {qubits[3], qubits[9]});
    circ.add_op<UnitID>(OpType::CX, {qubits[10], qubits[0]});
    circ.add_op<UnitID>(OpType::CX, {qubits[6], qubits[0]});

    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]}, {qubits[1], nodes[1]},  {qubits[2], nodes[2]},
        {qubits[3], nodes[3]}, {qubits[4], nodes[4]},  {qubits[5], nodes[6]},
        {qubits[6], nodes[5]}, {qubits[7], nodes[7]},  {qubits[8], nodes[8]},
        {qubits[9], nodes[9]}, {qubits[10], nodes[10]}};
    circ.rename_units(rename_map);

    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);
    LexiRouteRoutingMethod lrrm(100);
    unit_map_t init_map = lrrm.routing_method(mf, shared_arc);
    REQUIRE(init_map.size() == 0);
    std::vector<Command> commands = mf->circuit_.get_commands();
    REQUIRE(commands.size() == 10);
    Command swap_c = commands[0];
    unit_vector_t uids = {nodes[3], nodes[10]};
    REQUIRE(swap_c.get_args() == uids);
    REQUIRE(*swap_c.get_op_ptr() == *get_op_ptr(OpType::SWAP));
  }
}
SCENARIO("Test MappingManager::route_circuit with lc_route_subcircuit") {
  GIVEN("11 Node Architecture, 11 Qubit circuit, multiple SWAP required.") {
    std::vector<Node> nodes = {
        Node("test_node", 0), Node("test_node", 1), Node("test_node", 2),
        Node("node_test", 3), Node("node_test", 4), Node("node_test", 5),
        Node("test_node", 6), Node("node_test", 7), Node("node_test", 8),
        Node("node_test", 9), Node("node_test", 10)};
    //       n9 -- n8 -- n10
    //             |     |
    // n0 -- n1 -- n2 -- n3 -- n4
    //             |     |
    //             n5    n7
    //             |
    //             n6
    Architecture architecture(
        {{nodes[0], nodes[1]},
         {nodes[1], nodes[2]},
         {nodes[2], nodes[3]},
         {nodes[3], nodes[4]},
         {nodes[2], nodes[5]},
         {nodes[5], nodes[6]},
         {nodes[3], nodes[7]},
         {nodes[2], nodes[8]},
         {nodes[8], nodes[9]},
         {nodes[8], nodes[10]},
         {nodes[3], nodes[10]}});
    ArchitecturePtr shared_arc = std::make_shared<Architecture>(architecture);
    Circuit circ(11);
    std::vector<Qubit> qubits = circ.all_qubits();
    for (unsigned i = 0; i < 10; i++) {
      circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[4]});
      circ.add_op<UnitID>(OpType::CX, {qubits[6], qubits[7]});
      circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[10]});
      circ.add_op<UnitID>(OpType::CX, {qubits[8], qubits[5]});
      circ.add_op<UnitID>(OpType::CX, {qubits[3], qubits[9]});
      circ.add_op<UnitID>(OpType::CX, {qubits[2], qubits[8]});

      circ.add_op<UnitID>(OpType::CX, {qubits[1], qubits[5]});
      circ.add_op<UnitID>(OpType::CX, {qubits[3], qubits[9]});
      circ.add_op<UnitID>(OpType::CX, {qubits[10], qubits[0]});
      circ.add_op<UnitID>(OpType::CX, {qubits[6], qubits[0]});
    }

    Circuit copy_circ(circ);
    // transform stuff
    PassPtr dec = gen_decompose_routing_gates_to_cxs_pass(architecture, false);

    MappingManager mm(shared_arc);
    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(copy_circ);

    std::vector<RoutingMethodPtr> vrm = {
        std::make_shared<LexiRouteRoutingMethod>(100)};
    REQUIRE(vrm[0]->check_method(mf, shared_arc));

    bool res = mm.route_circuit(circ, vrm);

    PredicatePtr routed_correctly =
        std::make_shared<ConnectivityPredicate>(architecture);
    PredicatePtrMap preds{CompilationUnit::make_type_pair(routed_correctly)};
    CompilationUnit cu0(circ, preds);
    dec->apply(cu0);
    REQUIRE(res);
    REQUIRE(cu0.check_all_predicates());
  }
  GIVEN("Square Grid Architecture, large number of gates.") {
    SquareGrid sg(5, 10);
    ArchitecturePtr shared_arc = std::make_shared<Architecture>(sg);
    Circuit circ(35);
    std::vector<Qubit> qubits = circ.all_qubits();
    for (unsigned i = 0; i < qubits.size() - 1; i++) {
      circ.add_op<UnitID>(OpType::CX, {qubits[i], qubits[i + 1]});
    }
    for (unsigned i = 0; i < qubits.size() - 2; i++) {
      circ.add_op<UnitID>(OpType::CZ, {qubits[i], qubits[i + 2]});
    }
    // transform stuff
    PassPtr dec = gen_decompose_routing_gates_to_cxs_pass(sg, false);

    MappingManager mm(shared_arc);
    std::vector<RoutingMethodPtr> vrm = {
        std::make_shared<LexiRouteRoutingMethod>(100)};
    bool res = mm.route_circuit(circ, vrm);

    PredicatePtr routed_correctly = std::make_shared<ConnectivityPredicate>(sg);
    PredicatePtrMap preds{CompilationUnit::make_type_pair(routed_correctly)};
    CompilationUnit cu(circ, preds);
    dec->apply(cu);
    REQUIRE(res);
    REQUIRE(cu.check_all_predicates());
    REQUIRE(circ.n_gates() == 88);
  }
}

SCENARIO(
    "Check that an already solved routing problem will not add unecessary "
    "swaps") {
  GIVEN("A solved problem") {
    Circuit test_circuit;
    test_circuit.add_blank_wires(4);
    add_2qb_gates(test_circuit, OpType::CX, {{0, 1}, {1, 2}, {2, 3}, {3, 0}});

    // Ring of size 4
    RingArch arc(4);
    MappingManager mm(std::make_shared<Architecture>(arc));
    REQUIRE(mm.route_circuit(
        test_circuit, {std::make_shared<LexiLabellingMethod>(),
                       std::make_shared<LexiRouteRoutingMethod>()}));
    REQUIRE(test_circuit.n_gates() == 4);
  }
  GIVEN("A solved problem supplied with map and custom architecture") {
    Circuit test_circuit;
    test_circuit.add_blank_wires(4);
    add_2qb_gates(test_circuit, OpType::CX, {{0, 1}, {1, 2}, {2, 3}, {3, 0}});

    Architecture test_arc({{0, 1}, {1, 2}, {2, 3}, {3, 0}});
    Placement test_p(test_arc);

    qubit_mapping_t map_;
    for (unsigned nn = 0; nn <= 3; ++nn) {
      map_[Qubit(nn)] = Node(nn);
    }
    test_p.place_with_map(test_circuit, map_);
    qubit_vector_t all_qs_post_place = test_circuit.all_qubits();

    MappingManager mm(std::make_shared<Architecture>(test_arc));
    REQUIRE(!mm.route_circuit(
        test_circuit, {std::make_shared<LexiLabellingMethod>(),
                       std::make_shared<LexiRouteRoutingMethod>()}));

    qubit_vector_t all_qs_post_solve = test_circuit.all_qubits();
    REQUIRE(all_qs_post_place == all_qs_post_solve);
    REQUIRE(test_circuit.n_gates() == 4);
  }
}

SCENARIO("Empty Circuit test") {
  GIVEN("An Empty Circuit") {
    Circuit circ;
    circ.add_blank_wires(4);
    Architecture arc({{0, 1}, {1, 2}, {2, 3}});
    MappingManager mm(std::make_shared<Architecture>(arc));
    REQUIRE(!mm.route_circuit(
        circ, {std::make_shared<LexiLabellingMethod>(),
               std::make_shared<LexiRouteRoutingMethod>()}));
    REQUIRE(circ.n_gates() == 0);
  }
}

SCENARIO("Routing on circuit with no multi-qubit gates") {
  GIVEN("A circuit with no multi-qubit gates") {
    Circuit circ;
    circ.add_blank_wires(4);
    add_1qb_gates(circ, OpType::X, {0, 2});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::Y, {1});

    unsigned orig_vertices = circ.n_vertices();
    Architecture arc({{0, 1}, {1, 2}, {2, 3}});
    MappingManager mm(std::make_shared<Architecture>(arc));
    REQUIRE(!mm.route_circuit(
        circ, {std::make_shared<LexiLabellingMethod>(),
               std::make_shared<LexiRouteRoutingMethod>()}));
    REQUIRE(orig_vertices - 8 == circ.n_gates());
  }
}

SCENARIO("Test routing on a directed architecture with bidirectional edges") {
  GIVEN("A simple two-qubit circuit") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    Architecture arc({{0, 1}, {1, 0}});
    Architecture arc2(std::vector<std::pair<unsigned, unsigned>>{{0, 1}});

    // routing ignored bi directional edge and solves correctly
    MappingManager mm(std::make_shared<Architecture>(arc));
    REQUIRE(mm.route_circuit(
        circ, {std::make_shared<LexiLabellingMethod>(),
               std::make_shared<LexiRouteRoutingMethod>()}));
    REQUIRE(circ.n_gates() == 2);
    CHECK(respects_connectivity_constraints(circ, arc, false));
  }
}

SCENARIO(
    "Test routing on a directed architecture doesn't throw an error if "
    "non-cx optype is presented") {
  GIVEN(
      "A simple two-qubit circuit with non-cx multi-qubit gates and a "
      "directed architecture") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CU1, 0.5, {1, 0});
    circ.add_op<unsigned>(OpType::CU1, 0.5, {0, 1});
    circ.add_op<unsigned>(OpType::CY, {1, 0});
    circ.add_op<unsigned>(OpType::CY, {0, 1});
    circ.add_op<unsigned>(OpType::CZ, {1, 0});
    circ.add_op<unsigned>(OpType::CZ, {0, 1});
    circ.add_op<unsigned>(OpType::CRz, 0.5, {1, 0});
    circ.add_op<unsigned>(OpType::CRz, 0.5, {0, 1});

    Architecture arc(std::vector<std::pair<unsigned, unsigned>>{{0, 1}});
    MappingManager mm(std::make_shared<Architecture>(arc));
    REQUIRE(mm.route_circuit(
        circ, {std::make_shared<LexiLabellingMethod>(),
               std::make_shared<LexiRouteRoutingMethod>()}));
    REQUIRE(circ.n_gates() == 8);
  }
}

SCENARIO("Dense CX circuits route succesfully") {
  GIVEN(
      "Complex CX circuits for large directed architecture based off "
      "IBMTokyo") {
    Circuit circ(17);
    for (unsigned x = 0; x < 17; ++x) {
      for (unsigned y = 0; y + 1 < x; ++y) {
        if (x % 2) {  // swap the way directed chain runs each time
          add_2qb_gates(circ, OpType::CX, {{x, y}, {y + 1, y}});
        } else {
          add_2qb_gates(circ, OpType::CX, {{y, x}, {y, y + 1}});
        }
      }
    }
    Architecture arc(
        {{0, 1},   {1, 2},   {2, 3},   {3, 4},   {0, 5},   {1, 6},   {1, 7},
         {2, 6},   {2, 7},   {3, 8},   {3, 9},   {4, 8},   {4, 9},   {5, 6},
         {5, 10},  {5, 11},  {6, 10},  {6, 11},  {6, 7},   {7, 12},  {7, 13},
         {7, 8},   {8, 12},  {8, 13},  {8, 9},   {10, 11}, {11, 16}, {11, 17},
         {11, 12}, {12, 16}, {12, 17}, {12, 13}, {13, 18}, {13, 19}, {13, 14},
         {14, 18}, {14, 19}, {15, 16}, {16, 17}, {17, 18}, {18, 19}});
    MappingManager mm(std::make_shared<Architecture>(arc));
    REQUIRE(mm.route_circuit(
        circ, {std::make_shared<LexiLabellingMethod>(),
               std::make_shared<LexiRouteRoutingMethod>()}));
    (Transforms::decompose_SWAP_to_CX() >> Transforms::decompose_BRIDGE_to_CX())
        .apply(circ);

    Transforms::decompose_CX_directed(arc).apply(circ);
    REQUIRE(respects_connectivity_constraints(circ, arc, true));
  }
}

SCENARIO(
    "Dense CX circuits route succesfully on undirected Ring with "
    "placement.") {
  GIVEN("Complex CX circuits, big ring") {
    Circuit circ(29);
    for (unsigned x = 0; x < 29; ++x) {
      for (unsigned y = 0; y + 1 < x; ++y) {
        if (x % 2) {
          add_2qb_gates(circ, OpType::CX, {{x, y}, {y + 1, y}});
        } else {
          add_2qb_gates(circ, OpType::CX, {{y, x}, {y, y + 1}});
        }
      }
    }
    RingArch arc(29);
    MappingManager mm(std::make_shared<Architecture>(arc));
    REQUIRE(mm.route_circuit(
        circ, {std::make_shared<LexiLabellingMethod>(),
               std::make_shared<LexiRouteRoutingMethod>()}));
    Transforms::decompose_SWAP_to_CX().apply(circ);
    REQUIRE(respects_connectivity_constraints(circ, arc, false, true));
  }
}

SCENARIO(
    "Dense CX circuits route succesfully on smart placement unfriendly "
    "architecture.") {
  GIVEN("Complex CX circuits, big ring") {
    Circuit circ(13);
    for (unsigned x = 0; x < 13; ++x) {
      for (unsigned y = 0; y + 1 < x; ++y) {
        if (x % 2) {
          add_2qb_gates(circ, OpType::CX, {{x, y}, {y + 1, y}});
        } else {
          add_2qb_gates(circ, OpType::CX, {{y, x}, {y, y + 1}});
        }
      }
    }
    Architecture arc(
        {{0, 1},
         {2, 0},
         {2, 4},
         {6, 4},
         {8, 6},
         {8, 10},
         {12, 10},
         {3, 1},
         {3, 5},
         {7, 5},
         {7, 9},
         {11, 9},
         {11, 13},
         {12, 13},
         {6, 7}});
    MappingManager mm(std::make_shared<Architecture>(arc));
    REQUIRE(mm.route_circuit(
        circ, {std::make_shared<LexiLabellingMethod>(),
               std::make_shared<LexiRouteRoutingMethod>()}));
    REQUIRE(respects_connectivity_constraints(circ, arc, false, true));
  }
}

SCENARIO("Empty circuits, with and without blank wires") {
  GIVEN("An empty circuit with some qubits") {
    Circuit circ(6);
    RingArch arc(6);
    MappingManager mm(std::make_shared<Architecture>(arc));
    REQUIRE(!mm.route_circuit(
        circ, {std::make_shared<LexiLabellingMethod>(),
               std::make_shared<LexiRouteRoutingMethod>()}));
    REQUIRE(circ.depth() == 0);
    REQUIRE(circ.n_gates() == 0);
    REQUIRE(circ.n_qubits() == 6);
    REQUIRE(!respects_connectivity_constraints(circ, arc, true));
  }
  GIVEN("An empty circuit with no qubits") {
    Circuit circ(0);
    RingArch arc(6);
    MappingManager mm(std::make_shared<Architecture>(arc));
    REQUIRE(!mm.route_circuit(
        circ, {std::make_shared<LexiLabellingMethod>(),
               std::make_shared<LexiRouteRoutingMethod>()}));
    REQUIRE(circ.depth() == 0);
    REQUIRE(circ.n_gates() == 0);
    REQUIRE(circ.n_qubits() == 0);
  }
}

}  // namespace tket