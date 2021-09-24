// Copyright 2019-2021 Cambridge Quantum Computing
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

#include "Program/Program.hpp"
#include "testutil.hpp"

namespace tket {

SCENARIO("Basic Program construction") {
  GIVEN("An empty Program") {
    Program p;
    p.check_valid();
    REQUIRE(p.get_n_vertices() == 2);
  }
  GIVEN("A single block via add_block") {
    Program p(2, 2);
    Circuit c(2, 2);
    c.add_op(OpType::X, uvec{0});
    c.add_measure(0, 0);
    p.add_block(c);
    p.check_valid();
    REQUIRE(p.get_n_vertices() == 3);
  }
  GIVEN("A single block via add_op") {
    Program p(2, 2);
    p.add_op(OpType::X, uvec{0});
    p.add_op(OpType::Measure, uvec{1, 1});
    p.check_valid();
    REQUIRE(p.get_n_vertices() == 3);
    REQUIRE(p.bit_readout().size() == 2);
    REQUIRE(p.qubit_readout().size() == 1);
    REQUIRE(p.qubit_readout().at(Qubit(1)) == 1);
  }
  GIVEN("A straight-line block sequence") {
    Program p(2, 2);
    Circuit c(2, 2);
    c.add_op(OpType::X, uvec{0});
    c.add_measure(0, 0);
    p.add_block(c);
    p.add_block(c);
    p.add_block(c);
    p.check_valid();
    REQUIRE(p.get_n_vertices() == 5);
  }
  GIVEN("Conditional execution") {
    Program p(2, 2);
    p.add_op(OpType::X, uvec{0});
    p.add_op(OpType::Measure, uvec{0, 0});
    Program body(2, 2);
    body.add_op(OpType::X, uvec{1});
    p.append_if(Bit(0), body);
    p.add_op(OpType::Z, uvec{0});
    p.check_valid();
    REQUIRE(p.get_n_vertices() == 6);
  }
  GIVEN("If-else") {
    Program p(2, 2);
    p.add_op(OpType::X, uvec{0});
    p.add_op(OpType::Measure, uvec{0, 0});
    Program ifbody(2, 2);
    ifbody.add_op(OpType::X, uvec{1});
    Program elsebody(2, 2);
    elsebody.add_op(OpType::Y, uvec{1});
    p.append_if_else(Bit(0), ifbody, elsebody);
    p.add_op(OpType::Z, uvec{0});
    p.check_valid();
    REQUIRE(p.get_n_vertices() == 8);
  }
  GIVEN("While loop") {
    Program p(2, 2);
    p.add_op(OpType::X, uvec{0});
    p.add_op(OpType::Measure, uvec{0, 0});
    Program whilebody(2, 2);
    whilebody.add_op(OpType::H, uvec{0});
    whilebody.add_op(OpType::Measure, uvec{0, 0});
    p.append_while(Bit(0), whilebody);
    p.add_op(OpType::Y, uvec{1});
    p.check_valid();
    REQUIRE(p.get_n_vertices() == 7);
  }
  GIVEN("Appending interesting programs") {
    Program p(2, 2);
    Program ifbody(2, 2);
    ifbody.add_op(OpType::X, uvec{1});
    p.append_if(Bit(0), ifbody);
    Program p2(2, 2);
    Program whilebody(2, 2);
    whilebody.add_op(OpType::H, uvec{0});
    whilebody.add_op(OpType::Measure, uvec{0, 0});
    p2.append_while(Bit(0), whilebody);
    p.append(p2);
    p.check_valid();
    REQUIRE(p.get_n_vertices() == 7);
  }
}

SCENARIO("Check producing Graphviz does not fail") {
  GIVEN("If-else") {
    Program p(2, 2);
    p.add_op(OpType::X, uvec{0});
    p.add_op(OpType::Measure, uvec{0, 0});
    Program ifbody(2, 2);
    ifbody.add_op(OpType::X, uvec{1});
    Program elsebody(2, 2);
    elsebody.add_op(OpType::Y, uvec{1});
    p.append_if_else(Bit(0), ifbody, elsebody);
    p.add_op(OpType::Z, uvec{0});
    REQUIRE_NOTHROW(p.to_graphviz_file("ifelse.dot"));
    remove("ifelse.dot");
  }
}

SCENARIO("Block iteration") {
  GIVEN("An empty Program") {
    Program p;
    REQUIRE(p.block_begin() == p.block_end());
  }
  GIVEN("A single block via add_block") {
    Program p(2, 2);
    Circuit c(2, 2);
    c.add_op(OpType::X, uvec{0});
    c.add_measure(0, 0);
    FGVert block = p.add_block(c);
    Program::BlockIterator blit = p.block_begin();
    CHECK(*blit == block);
    ++blit;
    CHECK(blit == p.block_end());
  }
  GIVEN("A single block via add_op") {
    Program p(2, 2);
    p.add_op(OpType::X, uvec{0});
    p.add_op(OpType::Measure, uvec{0, 0});
    Program::BlockIterator blit = p.block_begin();
    CHECK(blit != p.block_end());
    ++blit;
    CHECK(blit == p.block_end());
  }
  GIVEN("A straight-line block sequence") {
    Program p(2, 2);
    Circuit c(2, 2);
    c.add_op(OpType::X, uvec{0});
    c.add_measure(0, 0);
    FGVert block0 = p.add_block(c);
    FGVert block1 = p.add_block(c);
    FGVert block2 = p.add_block(c);
    Program::BlockIterator blit = p.block_begin();
    CHECK(*blit == block0);
    ++blit;
    CHECK(*blit == block1);
    ++blit;
    CHECK(*blit == block2);
    ++blit;
    CHECK(blit == p.block_end());
  }
  GIVEN("Conditional execution") {
    Program p(2, 2);
    p.add_op(OpType::X, uvec{0});
    p.add_op(OpType::Measure, uvec{0, 0});
    Program body(2, 2);
    body.add_op(OpType::X, uvec{1});
    p.append_if(Bit(0), body);
    p.add_op(OpType::Z, uvec{0});
    Program::BlockIterator blit = p.block_begin();
    // Initial gates
    CHECK(blit.get_circuit_ref().n_gates() == 2);
    ++blit;
    // Block is just a branch
    CHECK(blit.get_condition());
    ++blit;
    // Block AFTER the if
    CHECK(blit.get_circuit_ref().n_gates() == 1);
    CHECK(p.n_in_edges(*blit) == 2);
    ++blit;
    // Body of the if
    CHECK(blit.get_circuit_ref().n_gates() == 1);
    CHECK(p.n_in_edges(*blit) == 1);
    ++blit;
    // Exit
    CHECK(blit == p.block_end());
  }
  GIVEN("If-else") {
    Program p(2, 2);
    p.add_op(OpType::X, uvec{0});
    p.add_op(OpType::Measure, uvec{0, 0});
    Program ifbody(2, 2);
    ifbody.add_op(OpType::X, uvec{1});
    Program elsebody(1, 1);
    elsebody.add_op(OpType::Y, uvec{0});
    p.append_if_else(Bit(0), ifbody, elsebody);
    p.add_op(OpType::Z, uvec{0});
    Program::BlockIterator blit = p.block_begin();
    // Initial gates
    CHECK(blit.get_circuit_ref().n_gates() == 2);
    ++blit;
    // Block is just a branch
    CHECK(blit.get_condition());
    ++blit;
    // Else block
    CHECK(blit.get_circuit_ref().n_qubits() == 1);
    ++blit;
    // Block after if-else
    CHECK(blit.get_circuit_ref().n_gates() == 1);
    CHECK(p.n_in_edges(*blit) == 2);
    ++blit;
    // Body of the if
    CHECK(blit.get_circuit_ref().n_gates() == 1);
    CHECK(p.n_in_edges(*blit) == 1);
    ++blit;
    // Exit of if
    CHECK(blit.get_circuit_ref().n_gates() == 0);
    CHECK(p.n_in_edges(*blit) == 1);
    ++blit;
    // Exit
    CHECK(blit == p.block_end());
  }
  GIVEN("While loop") {
    Program p(2, 2);
    p.add_op(OpType::X, uvec{0});
    p.add_op(OpType::Y, uvec{1});
    p.add_op(OpType::Measure, uvec{0, 0});
    Program whilebody(2, 2);
    whilebody.add_op(OpType::H, uvec{0});
    whilebody.add_op(OpType::Measure, uvec{0, 0});
    p.append_while(Bit(0), whilebody);
    p.add_op(OpType::Y, uvec{1});
    Program::BlockIterator blit = p.block_begin();
    // Initial gates
    CHECK(blit.get_circuit_ref().n_gates() == 3);
    ++blit;
    // Empty block
    CHECK(blit.get_circuit_ref().n_gates() == 0);
    CHECK(!blit.get_condition());
    ++blit;
    // Condition
    CHECK(blit.get_condition());
    ++blit;
    // Final command
    CHECK(blit.get_circuit_ref().n_gates() == 1);
    ++blit;
    // While body
    CHECK(blit.get_circuit_ref().n_gates() == 2);
    ++blit;
    // Exit
    CHECK(blit == p.block_end());
  }
  GIVEN("Appending interesting programs") {
    Bit b0(0);
    Bit b1(1);
    Program p(2, 2);
    Program ifbody(2, 2);
    ifbody.add_op(OpType::X, uvec{1});
    p.append_if(b1, ifbody);
    Program p2(2, 2);
    Program whilebody(2, 2);
    whilebody.add_op(OpType::H, uvec{0});
    whilebody.add_op(OpType::Measure, uvec{0, 0});
    p2.append_while(b0, whilebody);
    p.append(p2);
    Program::BlockIterator blit = p.block_begin();
    // Initial empty block for if condition
    CHECK(blit.get_circuit_ref().n_gates() == 0);
    CHECK(blit.get_condition() == b1);
    ++blit;
    // Empty block after if
    CHECK(blit.get_circuit_ref().n_gates() == 0);
    CHECK_FALSE(blit.get_condition());
    ++blit;
    // While condition
    CHECK(blit.get_circuit_ref().n_gates() == 0);
    CHECK(blit.get_condition() == b0);
    ++blit;
    // While body
    CHECK(blit.get_circuit_ref().n_gates() == 2);
    ++blit;
    // If body
    CHECK(blit.get_circuit_ref().n_gates() == 1);
    ++blit;
    // Exit
    CHECK(blit == p.block_end());
  }
}

static void test_command_types_sequence(
    const Program& p, const std::vector<OpType>& expected_types) {
  Program::CommandIterator cit = p.begin();
  for (auto type : expected_types) {
    CHECK(cit->get_op_ptr()->get_type() == type);
    ++cit;
  }
  CHECK(cit == p.end());
}

SCENARIO("Command iteration") {
  GIVEN("An empty Program") {
    Program p;
    test_command_types_sequence(p, {OpType::Stop});
  }
  GIVEN("A single block via add_op") {
    Program p(2, 2);
    p.add_op(OpType::X, uvec{0});
    p.add_op(OpType::Measure, uvec{0, 0});
    test_command_types_sequence(p, {OpType::X, OpType::Measure, OpType::Stop});
  }
  GIVEN("A straight-line block sequence") {
    unsigned N = 3;
    Program p(2, 2);
    Circuit c(2, 2);
    c.add_op(OpType::X, uvec{0});
    c.add_measure(0, 0);
    for (unsigned i = 0; i < N; ++i) {
      p.add_block(c);
    }
    Program::CommandIterator cit = p.begin();
    for (unsigned i = 0; i < N; ++i) {
      CHECK(cit->get_op_ptr()->get_type() == OpType::X);
      ++cit;
      CHECK(cit->get_op_ptr()->get_type() == OpType::Measure);
      ++cit;
    }
    CHECK(cit->get_op_ptr()->get_type() == OpType::Stop);
    ++cit;
    CHECK(cit == p.end());
  }
  GIVEN("Conditional execution") {
    Program p(2, 2);
    p.add_op(OpType::X, uvec{0});
    p.add_op(OpType::Measure, uvec{0, 0});
    Program body(2, 2);
    body.add_op(OpType::X, uvec{1});
    p.append_if(Bit(0), body);
    p.add_op(OpType::Z, uvec{0});

    test_command_types_sequence(
        p, {// Initial gates
            OpType::X, OpType::Measure,
            // Block is just a branch
            OpType::Branch,
            // Block AFTER the if
            OpType::Label, OpType::Z, OpType::Goto,
            // Body of the if
            OpType::Label, OpType::X, OpType::Goto,
            // Exit
            OpType::Label, OpType::Stop});
  }
  GIVEN("If-else") {
    Program p(2, 2);
    p.add_op(OpType::X, uvec{0});
    p.add_op(OpType::Measure, uvec{0, 0});
    Program ifbody(2, 2);
    ifbody.add_op(OpType::X, uvec{1});
    Program elsebody(1, 1);
    elsebody.add_op(OpType::Y, uvec{0});
    p.append_if_else(Bit(0), ifbody, elsebody);
    p.add_op(OpType::Z, uvec{0});

    test_command_types_sequence(
        p, {// Initial gates
            OpType::X, OpType::Measure,
            // Block is just a branch
            OpType::Branch,
            // Else block
            OpType::Y,
            // Block after if-else
            OpType::Label, OpType::Z, OpType::Goto,
            // Body of the if
            OpType::Label, OpType::X, OpType::Goto,
            // Exit
            OpType::Label, OpType::Stop});
  }
  GIVEN("While loop") {
    Program p(2, 2);
    p.add_op(OpType::X, uvec{0});
    p.add_op(OpType::Y, uvec{1});
    p.add_op(OpType::Measure, uvec{0, 0});
    Program whilebody(2, 2);
    whilebody.add_op(OpType::H, uvec{0});
    whilebody.add_op(OpType::Measure, uvec{0, 0});
    p.append_while(Bit(0), whilebody);
    p.add_op(OpType::Y, uvec{1});

    test_command_types_sequence(
        p, {// Initial gates
            OpType::X, OpType::Y, OpType::Measure,
            // Condition
            OpType::Label, OpType::Branch,
            // Final command
            OpType::Y, OpType::Goto,
            // While body
            OpType::Label, OpType::H, OpType::Measure, OpType::Goto,
            // Exit
            OpType::Label, OpType::Stop});
  }
  GIVEN("Appending interesting programs") {
    Bit b0(0);
    Bit b1(1);
    Program p(2, 2);
    Program ifbody(2, 2);
    ifbody.add_op(OpType::X, uvec{1});
    p.append_if(b1, ifbody);
    Program p2(2, 2);
    Program whilebody(2, 2);
    whilebody.add_op(OpType::H, uvec{0});
    whilebody.add_op(OpType::Measure, uvec{0, 0});
    p2.append_while(b0, whilebody);
    p.append(p2);

    test_command_types_sequence(
        p, {// Initial empty block for if condition
            OpType::Branch,
            // Empty block after if
            OpType::Label,
            // While condition
            OpType::Label, OpType::Branch, OpType::Goto,
            // While body
            OpType::Label, OpType::H, OpType::Measure, OpType::Goto,
            // If body
            OpType::Label, OpType::X, OpType::Goto,
            // Exit
            OpType::Label, OpType::Stop});
  }
}

}  // namespace tket
