#include <catch2/catch.hpp>

#include "Mapping/LexiLabelling.hpp"
#include "Mapping/LexiRoute.hpp"
#include "Mapping/MappingManager.hpp"
#include "Predicates/CompilationUnit.hpp"
#include "Predicates/CompilerPass.hpp"
#include "Predicates/PassGenerators.hpp"
#include "Predicates/PassLibrary.hpp"

namespace tket {

SCENARIO("Test LexiRoute::solve and LexiRoute::solve_labelling") {
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
    lr.solve_labelling();
    // lr.solve(4);

    REQUIRE(mf0->circuit_.n_gates() == 3);

    rename_map = {{qubits[4], nodes[6]}};
    mf0->circuit_.rename_units(rename_map);

    std::shared_ptr<MappingFrontier> mf1 =
        std::make_shared<MappingFrontier>(circ);
    LexiRoute lr1(shared_arc, mf1);
    // lr1.solve_labelling();
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
    lr0.solve_labelling();
    // lr0.solve(20);
    std::vector<Command> commands = mf->circuit_.get_commands();
    REQUIRE(commands.size() == 4);
    Command c = commands[0];
    unit_vector_t uids = {nodes[2], nodes[1]};
    REQUIRE(c.get_args() == uids);
    mf->advance_frontier_boundary(shared_arc);

    LexiRoute lr1(shared_arc, mf);
    lr1.solve_labelling();
    uids = {nodes[2], nodes[3]};
    REQUIRE(mf->circuit_.get_commands()[1].get_args() == uids);
    mf->advance_frontier_boundary(shared_arc);

    LexiRoute lr2(shared_arc, mf);
    lr2.solve_labelling();
    uids = {nodes[2], nodes[5]};
    REQUIRE(mf->circuit_.get_commands()[2].get_args() == uids);
    mf->advance_frontier_boundary(shared_arc);

    LexiRoute lr3(shared_arc, mf);
    lr3.solve_labelling();
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
    lr1.solve_labelling();
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
    lr0.solve_labelling();

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
    lr0.solve_labelling();

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
      "Labelling is required, but there are no free remaining qubits, for"
      "one updated label, order 0.") {
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
    REQUIRE_THROWS_AS(lr.solve_labelling(), LexiRouteError);
  }
  GIVEN(
      "Labelling is required, but there are no free remaining qubits, for "
      " one updated label, order 1.") {
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
    REQUIRE_THROWS_AS(lr.solve_labelling(), LexiRouteError);
  }
  GIVEN(
      "Labelling is required, but there are no free remaining qubits, for"
      "two updated labels.") {
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
    REQUIRE_THROWS_AS(lr.solve_labelling(), LexiRouteError);
  }
}

SCENARIO("Test LexiLabellingMethod") {
  std::vector<Node> nodes = {
      Node("test_node", 0), Node("test_node", 1), Node("test_node", 2),
      Node("node_test", 3), Node("node_test", 4)};

  // straight line
  Architecture architecture(
      {{nodes[0], nodes[1]},
       {nodes[1], nodes[2]},
       {nodes[2], nodes[3]},
       {nodes[3], nodes[4]}});
  ArchitecturePtr shared_arc = std::make_shared<Architecture>(architecture);
  GIVEN("No qubit to label, empty frontier, check_method.") {
    Circuit circ(5);
    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);
    LexiLabellingMethod lrm;
    REQUIRE(!lrm.check_method(mf, shared_arc));
  }
  GIVEN("No qubit to label, partially filled frontier, check_method.") {
    Circuit circ(5);
    std::vector<Qubit> qubits = circ.all_qubits();
    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[4]});
    circ.add_op<UnitID>(OpType::CZ, {qubits[1], qubits[2]});
    circ.add_op<UnitID>(OpType::ZZPhase, 0.3, {qubits[3], qubits[0]});
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]},
        {qubits[1], nodes[1]},
        {qubits[2], nodes[2]},
        {qubits[3], nodes[3]},
        {qubits[4], nodes[4]}};
    circ.rename_units(rename_map);
    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);
    LexiLabellingMethod lrm;
    REQUIRE(!lrm.check_method(mf, shared_arc));
  }
  GIVEN("Qubit to label, but casually restricted, check_method.") {
    Circuit circ(5);
    std::vector<Qubit> qubits = circ.all_qubits();
    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[4]});
    circ.add_op<UnitID>(OpType::CZ, {qubits[1], qubits[2]});
    circ.add_op<UnitID>(OpType::ZZPhase, 0.3, {qubits[3], qubits[0]});
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]},
        {qubits[1], nodes[1]},
        {qubits[2], nodes[2]},
        {qubits[4], nodes[4]}};
    circ.rename_units(rename_map);
    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);
    LexiLabellingMethod lrm;
    REQUIRE(!lrm.check_method(mf, shared_arc));
  }
  GIVEN(
      "Two Qubit to label in future slice, causally restricted, "
      "check_method.") {
    Circuit circ(5);
    std::vector<Qubit> qubits = circ.all_qubits();
    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[1]});
    circ.add_op<UnitID>(OpType::CZ, {qubits[1], qubits[2]});
    circ.add_op<UnitID>(OpType::CZ, {qubits[2], qubits[3]});
    circ.add_op<UnitID>(OpType::ZZPhase, 0.3, {qubits[3], qubits[4]});
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]}, {qubits[1], nodes[1]}, {qubits[2], nodes[2]}};
    circ.rename_units(rename_map);
    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);
    LexiLabellingMethod lrm;
    REQUIRE(!lrm.check_method(mf, shared_arc));
  }
  GIVEN("Three Qubit Gate, all labelled, first slice, check_method.") {
    Circuit circ(5);
    std::vector<Qubit> qubits = circ.all_qubits();
    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[4]});
    circ.add_op<UnitID>(OpType::CCX, {qubits[1], qubits[2], qubits[3]});
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]},
        {qubits[1], nodes[1]},
        {qubits[2], nodes[2]},
        {qubits[3], nodes[3]},
        {qubits[4], nodes[4]}};
    circ.rename_units(rename_map);
    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);
    LexiLabellingMethod lrm;
    REQUIRE(!lrm.check_method(mf, shared_arc));
  }
  GIVEN("One unlabelled qubit, one slice, check and route.") {
    Circuit circ(5);
    std::vector<Qubit> qubits = circ.all_qubits();
    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[1]});
    circ.add_op<UnitID>(OpType::CX, {qubits[2], qubits[3]});
    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]}, {qubits[1], nodes[1]}, {qubits[2], nodes[2]}};
    circ.rename_units(rename_map);
    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);
    VertPort pre_label =
        mf->quantum_boundary->get<TagKey>().find(qubits[3])->second;
    LexiLabellingMethod lrm;
    REQUIRE(lrm.check_method(mf, shared_arc));
    lrm.routing_method(mf, shared_arc);
    REQUIRE(
        mf->quantum_boundary->get<TagKey>().find(qubits[3]) ==
        mf->quantum_boundary->get<TagKey>().end());
    VertPort post_label =
        mf->quantum_boundary->get<TagKey>().find(nodes[3])->second;
    REQUIRE(pre_label == post_label);
  }
  GIVEN(
      "One unlabelled qubit, two slices, lookahead for better solution, check"
      " and route.") {
    Circuit circ(5);
    std::vector<Qubit> qubits = circ.all_qubits();
    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[1]});
    circ.add_op<UnitID>(OpType::ZZPhase, 0.8, {qubits[2], qubits[3]});
    circ.add_op<UnitID>(OpType::CZ, {qubits[2], qubits[0]});

    std::map<UnitID, UnitID> rename_map = {
        {qubits[0], nodes[0]}, {qubits[1], nodes[1]}, {qubits[3], nodes[3]}};
    circ.rename_units(rename_map);
    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);
    VertPort pre_label =
        mf->quantum_boundary->get<TagKey>().find(qubits[2])->second;
    LexiLabellingMethod lrm;
    REQUIRE(lrm.check_method(mf, shared_arc));
    lrm.routing_method(mf, shared_arc);
    REQUIRE(
        mf->quantum_boundary->get<TagKey>().find(qubits[2]) ==
        mf->quantum_boundary->get<TagKey>().end());
    VertPort post_label =
        mf->quantum_boundary->get<TagKey>().find(nodes[2])->second;
    REQUIRE(pre_label == post_label);
  }
  GIVEN("Two unlabelled qubits, one slice, check and route.") {
    Circuit circ(5);
    std::vector<Qubit> qubits = circ.all_qubits();
    circ.add_op<UnitID>(OpType::CX, {qubits[0], qubits[1]});
    circ.add_op<UnitID>(OpType::ZZPhase, 0.8, {qubits[2], qubits[3]});

    std::map<UnitID, UnitID> rename_map = {
        {qubits[2], nodes[2]}, {qubits[1], nodes[1]}};
    circ.rename_units(rename_map);
    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);
    VertPort pre_label_0 =
        mf->quantum_boundary->get<TagKey>().find(qubits[0])->second;
    VertPort pre_label_3 =
        mf->quantum_boundary->get<TagKey>().find(qubits[3])->second;
    LexiLabellingMethod lrm;
    REQUIRE(lrm.check_method(mf, shared_arc));
    lrm.routing_method(mf, shared_arc);
    REQUIRE(
        mf->quantum_boundary->get<TagKey>().find(qubits[0]) ==
        mf->quantum_boundary->get<TagKey>().end());
    REQUIRE(
        mf->quantum_boundary->get<TagKey>().find(qubits[3]) ==
        mf->quantum_boundary->get<TagKey>().end());
    VertPort post_label_0 =
        mf->quantum_boundary->get<TagKey>().find(nodes[0])->second;
    REQUIRE(pre_label_0 == post_label_0);
    VertPort post_label_3 =
        mf->quantum_boundary->get<TagKey>().find(nodes[3])->second;
    REQUIRE(pre_label_3 == post_label_3);
  }
  GIVEN("Two unlabelled qubits, two slices, lookahead, check and route.") {
    Circuit circ(5);
    std::vector<Qubit> qubits = circ.all_qubits();
    circ.add_op<UnitID>(OpType::CX, {qubits[2], qubits[1]});
    circ.add_op<UnitID>(OpType::ZZPhase, 0.8, {qubits[4], qubits[3]});
    circ.add_op<UnitID>(OpType::CX, {qubits[2], qubits[4]});

    std::map<UnitID, UnitID> rename_map = {
        {qubits[4], nodes[4]}, {qubits[1], nodes[1]}};
    circ.rename_units(rename_map);
    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);
    VertPort pre_label_0 =
        mf->quantum_boundary->get<TagKey>().find(qubits[2])->second;
    VertPort pre_label_3 =
        mf->quantum_boundary->get<TagKey>().find(qubits[3])->second;
    LexiLabellingMethod lrm;
    REQUIRE(lrm.check_method(mf, shared_arc));
    lrm.routing_method(mf, shared_arc);
    REQUIRE(
        mf->quantum_boundary->get<TagKey>().find(qubits[2]) ==
        mf->quantum_boundary->get<TagKey>().end());
    REQUIRE(
        mf->quantum_boundary->get<TagKey>().find(qubits[3]) ==
        mf->quantum_boundary->get<TagKey>().end());
    VertPort post_label_0 =
        mf->quantum_boundary->get<TagKey>().find(nodes[0])->second;
    REQUIRE(pre_label_0 == post_label_0);
    VertPort post_label_3 =
        mf->quantum_boundary->get<TagKey>().find(nodes[3])->second;
    REQUIRE(pre_label_3 == post_label_3);
  }
  GIVEN(
      "Two unlabelled qubits, two slices, lookahead unrouted, check and "
      "route.") {
    Circuit circ(5);
    std::vector<Qubit> qubits = circ.all_qubits();
    circ.add_op<UnitID>(OpType::CX, {qubits[2], qubits[1]});
    circ.add_op<UnitID>(OpType::ZZPhase, 0.8, {qubits[4], qubits[3]});
    circ.add_op<UnitID>(OpType::CX, {qubits[2], qubits[0]});

    std::map<UnitID, UnitID> rename_map = {
        {qubits[4], nodes[4]}, {qubits[1], nodes[1]}};
    circ.rename_units(rename_map);
    std::shared_ptr<MappingFrontier> mf =
        std::make_shared<MappingFrontier>(circ);
    VertPort pre_label_0 =
        mf->quantum_boundary->get<TagKey>().find(qubits[2])->second;
    VertPort pre_label_3 =
        mf->quantum_boundary->get<TagKey>().find(qubits[3])->second;
    LexiLabellingMethod lrm;
    REQUIRE(lrm.check_method(mf, shared_arc));
    lrm.routing_method(mf, shared_arc);
    REQUIRE(
        mf->quantum_boundary->get<TagKey>().find(qubits[2]) ==
        mf->quantum_boundary->get<TagKey>().end());
    REQUIRE(
        mf->quantum_boundary->get<TagKey>().find(qubits[3]) ==
        mf->quantum_boundary->get<TagKey>().end());
    VertPort post_label_0 =
        mf->quantum_boundary->get<TagKey>().find(nodes[0])->second;
    REQUIRE(pre_label_0 == post_label_0);
    VertPort post_label_3 =
        mf->quantum_boundary->get<TagKey>().find(nodes[3])->second;
    REQUIRE(pre_label_3 == post_label_3);
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
SCENARIO("Test MappingManager with LexiRouteRoutingMethod and LexiLabelling") {
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
    for (unsigned i = 0; i < 11; i++) {
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

    LexiLabellingMethod lrm;
    std::vector<RoutingMethodPtr> vrm = {
        std::make_shared<LexiLabellingMethod>(lrm),
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
    LexiLabellingMethod lrm;
    std::vector<RoutingMethodPtr> vrm = {
        std::make_shared<LexiLabellingMethod>(lrm),
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
}  // namespace tket