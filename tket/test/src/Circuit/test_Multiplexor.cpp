// Copyright 2019-2023 Cambridge Quantum Computing
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

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include "../testutil.hpp"
#include "Circuit/Boxes.hpp"
#include "Circuit/CircUtils.hpp"
#include "Circuit/Circuit.hpp"
#include "Circuit/Multiplexor.hpp"
#include "Eigen/src/Core/Matrix.h"
#include "Gate/SymTable.hpp"
#include "Simulation/CircuitSimulator.hpp"
#include "Utils/HelperFunctions.hpp"

namespace tket {
namespace test_Multiplexor {

static bool check_multiplexor(
    const ctrl_op_map_t &op_map, const Circuit &circ) {
  // Assume ops in op_map have the same number of q wires and have no c wires.
  // Also assumes op_map is not empty
  auto first_op = op_map.begin();
  unsigned n_ctrl_bits = (unsigned)first_op->first.size();
  unsigned total_ops = 1 << n_ctrl_bits;
  unsigned n_targets = first_op->second->n_qubits();
  std::vector<unsigned> target_qubits(n_targets);
  std::iota(std::begin(target_qubits), std::end(target_qubits), 0);
  unsigned block_size = 1 << n_targets;
  Eigen::MatrixXcd correct_u = Eigen::MatrixXcd::Identity(
      total_ops * block_size, total_ops * block_size);
  for (unsigned i = 0; i < total_ops; i++) {
    // Find the binary rep for i
    std::vector<bool> bin = dec_to_bin(i, n_ctrl_bits);
    auto it = op_map.find(bin);
    if (it != op_map.end()) {
      // get the matrix for op
      Op_ptr op = it->second;
      Circuit c(n_targets);
      c.add_op<unsigned>(op, target_qubits);
      c.decompose_boxes_recursively();
      Eigen::MatrixXcd block_m = tket_sim::get_unitary(c);
      correct_u.block(i * block_size, i * block_size, block_size, block_size) =
          block_m;
    }
  }
  Circuit circ_copy(circ);
  circ_copy.decompose_boxes_recursively();
  Eigen::MatrixXcd circ_u = tket_sim::get_unitary(circ_copy);
  return (correct_u - circ_u).cwiseAbs().sum() < ERR_EPS;
}

SCENARIO("MultiplexorBox decomposition", "[boxes]") {
  GIVEN("Simple MultiplexorBox construction") {
    Circuit c0(2);
    c0.add_op<unsigned>(OpType::H, {0});
    CircBox cbox(c0);
    Op_ptr op0 = std::make_shared<CircBox>(cbox);
    ctrl_op_map_t op_map = {
        {{1, 1}, op0},
        {{0, 1}, get_op_ptr(OpType::CX)},
        {{1, 0}, get_op_ptr(OpType::TK2, std::vector<Expr>{0.2, 0.4, 0.4})}};
    MultiplexorBox multiplexor(op_map);
    std::shared_ptr<Circuit> c = multiplexor.to_circuit();
    std::vector<Command> cmds = c->get_commands();
    // 4 X gates and 3 QControlBoxes
    REQUIRE(cmds.size() == 7);
    for (auto cmd : cmds) {
      REQUIRE(
          (cmd.get_op_ptr()->get_type() == OpType::QControlBox ||
           cmd.get_op_ptr()->get_type() == OpType::X));
    }
    REQUIRE(check_multiplexor(op_map, *c));
  }
  GIVEN("MultiplexorBox with one control") {
    ctrl_op_map_t op_map = {{{1}, get_op_ptr(OpType::H)}};
    MultiplexorBox multiplexor(op_map);
    std::shared_ptr<Circuit> c = multiplexor.to_circuit();
    std::vector<Command> cmds = c->get_commands();
    REQUIRE(cmds.size() == 1);
    REQUIRE(check_multiplexor(op_map, *c));
  }
  GIVEN("MultiplexorBox with zero control") {
    ctrl_op_map_t op_map = {{{}, get_op_ptr(OpType::H)}};
    MultiplexorBox multiplexor(op_map);
    std::shared_ptr<Circuit> c = multiplexor.to_circuit();
    std::vector<Command> cmds = c->get_commands();
    REQUIRE(cmds.size() == 1);
    REQUIRE(check_multiplexor(op_map, *c));
  }
}

SCENARIO("MultiplexedRotationBox decomposition", "[boxes]") {
  GIVEN("Controlled Rz construction") {
    ctrl_op_map_t op_map = {
        {{1, 1}, get_op_ptr(OpType::Rz, 0.3)},
        {{0, 1}, get_op_ptr(OpType::Rz, 1.4)},
        {{1, 0}, get_op_ptr(OpType::Rz, 0.7)}};
    MultiplexedRotationBox multiplexor(op_map);
    std::shared_ptr<Circuit> c = multiplexor.to_circuit();
    std::vector<Command> cmds = c->get_commands();
    REQUIRE(cmds.size() == 8);
    for (auto cmd : cmds) {
      REQUIRE(
          (cmd.get_op_ptr()->get_type() == OpType::Rz ||
           cmd.get_op_ptr()->get_type() == OpType::CX));
    }
    REQUIRE(check_multiplexor(op_map, *c));
  }
  GIVEN("Controlled Ry construction") {
    ctrl_op_map_t op_map = {
        {{1, 1, 0, 1, 0, 0}, get_op_ptr(OpType::Ry, 0.3)},
        {{0, 1, 1, 1, 1, 0}, get_op_ptr(OpType::Ry, 1.4)},
        {{1, 0, 1, 1, 1, 0}, get_op_ptr(OpType::Ry, 0.7)}};
    MultiplexedRotationBox multiplexor(op_map);
    std::shared_ptr<Circuit> c = multiplexor.to_circuit();
    std::vector<Command> cmds = c->get_commands();
    REQUIRE(cmds.size() == 128);
    for (auto cmd : cmds) {
      REQUIRE(
          (cmd.get_op_ptr()->get_type() == OpType::Ry ||
           cmd.get_op_ptr()->get_type() == OpType::CX));
    }
    REQUIRE(check_multiplexor(op_map, *c));
  }
  GIVEN("Controlled Rx construction") {
    ctrl_op_map_t op_map = {
        {{1, 1}, get_op_ptr(OpType::Rx, 0.3)},
        {{0, 1}, get_op_ptr(OpType::Rx, 1.4)},
        {{1, 0}, get_op_ptr(OpType::Rx, 0.7)}};
    MultiplexedRotationBox multiplexor(op_map);
    std::shared_ptr<Circuit> c = multiplexor.to_circuit();
    std::vector<Command> cmds = c->get_commands();
    REQUIRE(cmds.size() == 10);
    for (auto cmd : cmds) {
      REQUIRE(
          (cmd.get_op_ptr()->get_type() == OpType::H ||
           cmd.get_op_ptr()->get_type() == OpType::Rz ||
           cmd.get_op_ptr()->get_type() == OpType::CX));
    }
    REQUIRE(check_multiplexor(op_map, *c));
  }
  GIVEN("MultiplexedRotationBox with one control") {
    ctrl_op_map_t op_map = {
        {{1}, get_op_ptr(OpType::Rz, 0.3)}, {{0}, get_op_ptr(OpType::Rz, 1.4)}};
    MultiplexedRotationBox multiplexor(op_map);
    std::shared_ptr<Circuit> c = multiplexor.to_circuit();
    std::vector<Command> cmds = c->get_commands();
    REQUIRE(cmds.size() == 4);
    REQUIRE(check_multiplexor(op_map, *c));
  }
  GIVEN("MultiplexedRotationBox with zero control") {
    ctrl_op_map_t op_map = {{{}, get_op_ptr(OpType::Rx, 0.3)}};
    MultiplexedRotationBox multiplexor(op_map);
    std::shared_ptr<Circuit> c = multiplexor.to_circuit();
    std::vector<Command> cmds = c->get_commands();
    REQUIRE(cmds.size() == 1);
    REQUIRE(check_multiplexor(op_map, *c));
  }
  GIVEN("MultiplexedRotationBox with symbols") {
    Sym a = SymTable::fresh_symbol("a");
    Expr expr_a(a);
    Sym b = SymTable::fresh_symbol("b");
    Expr expr_b(b);
    ctrl_op_map_t op_map = {
        {{1, 1, 0}, get_op_ptr(OpType::Ry, expr_a)},
        {{0, 1, 1}, get_op_ptr(OpType::Ry, 1.4)},
        {{1, 0, 1}, get_op_ptr(OpType::Ry, expr_b)}};
    ctrl_op_map_t numerical_map = {
        {{1, 1, 0}, get_op_ptr(OpType::Ry, 0.3)},
        {{0, 1, 1}, get_op_ptr(OpType::Ry, 1.4)},
        {{1, 0, 1}, get_op_ptr(OpType::Ry, 1.8)}};
    MultiplexedRotationBox multiplexor(op_map);
    std::shared_ptr<Circuit> c_sb = multiplexor.to_circuit();
    std::vector<Command> cmds = c_sb->get_commands();
    REQUIRE(cmds.size() == 16);
    SymEngine::map_basic_basic smap;
    smap[a] = Expr(0.3);
    smap[b] = Expr(1.8);
    c_sb->symbol_substitution(smap);
    REQUIRE(check_multiplexor(numerical_map, *c_sb));
  }
}

SCENARIO("MultiplexedU2Box decomposition", "[boxes]") {
  GIVEN("MultiplexedU2Box with zero control") {
    ctrl_op_map_t op_map = {{{}, get_op_ptr(OpType::H)}};
    MultiplexedU2Box multiplexor(op_map);
    std::shared_ptr<Circuit> c = multiplexor.to_circuit();
    std::vector<Command> cmds = c->get_commands();
    REQUIRE(cmds.size() == 1);
    REQUIRE(check_multiplexor(op_map, *c));
  }
  GIVEN("MultiplexedU2Box with one control") {
    ctrl_op_map_t op_map = {
        {{1}, get_op_ptr(OpType::H)}, {{0}, get_op_ptr(OpType::X)}};
    MultiplexedU2Box multiplexor(op_map);
    std::shared_ptr<Circuit> c = multiplexor.to_circuit();
    std::vector<Command> cmds = c->get_commands();
    REQUIRE(cmds.size() == 3 + 1);
    for (auto cmd : cmds) {
      REQUIRE(
          (cmd.get_op_ptr()->get_type() == OpType::Unitary1qBox ||
           cmd.get_op_ptr()->get_type() == OpType::CX ||
           cmd.get_op_ptr()->get_type() == OpType::MultiplexedRotationBox));
    }
    REQUIRE(check_multiplexor(op_map, *c));
  }
  GIVEN("MultiplexedU2Box with two controls") {
    Circuit c0(1);
    c0.add_op<unsigned>(OpType::TK1, {0.2374, 1.0353, 0.5372}, {0});
    Eigen::Matrix2cd m = tket_sim::get_unitary(c0);
    Unitary1qBox mbox(m);
    Op_ptr mbox_op = std::make_shared<Unitary1qBox>(mbox);
    ctrl_op_map_t op_map = {
        {{1, 1}, mbox_op},
        {{0, 1}, get_op_ptr(OpType::X)},
        {{1, 0}, get_op_ptr(OpType::TK1, std::vector<Expr>{0.3, 1.8, 3.4})}};
    MultiplexedU2Box multiplexor(op_map);
    std::shared_ptr<Circuit> c = multiplexor.to_circuit();
    std::vector<Command> cmds = c->get_commands();
    REQUIRE(cmds.size() == 7 + 2);
    for (auto cmd : cmds) {
      REQUIRE(
          (cmd.get_op_ptr()->get_type() == OpType::Unitary1qBox ||
           cmd.get_op_ptr()->get_type() == OpType::CX ||
           cmd.get_op_ptr()->get_type() == OpType::MultiplexedRotationBox));
    }
    REQUIRE(check_multiplexor(op_map, *c));
  }
  GIVEN("MultiplexedU2Box with 5 controls") {
    ctrl_op_map_t op_map = {
        {{1, 1, 0, 1, 0}, get_op_ptr(OpType::H)},
        {{0, 1, 1, 0, 0}, get_op_ptr(OpType::X)},
        {{1, 0, 0, 0, 1}, get_op_ptr(OpType::Z)},
        {{1, 1, 0, 1, 1}, get_op_ptr(OpType::X)},
        {{1, 0, 1, 0, 1},
         get_op_ptr(OpType::TK1, std::vector<Expr>{0.3, 1.8, 3.4})},
    };
    MultiplexedU2Box multiplexor(op_map);
    std::shared_ptr<Circuit> c = multiplexor.to_circuit();
    std::vector<Command> cmds = c->get_commands();
    REQUIRE(cmds.size() == 63 + 5);
    for (auto cmd : cmds) {
      REQUIRE(
          (cmd.get_op_ptr()->get_type() == OpType::Unitary1qBox ||
           cmd.get_op_ptr()->get_type() == OpType::CX ||
           cmd.get_op_ptr()->get_type() == OpType::MultiplexedRotationBox));
    }
    REQUIRE(check_multiplexor(op_map, *c));
  }
  GIVEN("MultiplexedU2Box without the final diagonal") {
    ctrl_op_map_t op_map = {
        {{1}, get_op_ptr(OpType::H)}, {{0}, get_op_ptr(OpType::X)}};
    MultiplexedU2Box multiplexor(op_map, false);
    std::shared_ptr<Circuit> c = multiplexor.to_circuit();
    std::vector<Command> cmds = c->get_commands();
    REQUIRE(cmds.size() == 3);
  }
}

SCENARIO("Exception handling", "[boxes]") {
  GIVEN("Empty op_map") {
    ctrl_op_map_t op_map;
    REQUIRE_THROWS_MATCHES(
        MultiplexorBox(op_map), std::invalid_argument,
        MessageContains("No Ops provided"));
  }
  GIVEN("Classical wire") {
    Circuit c0(2, 1);
    c0.add_op<unsigned>(OpType::Measure, {0, 0});
    CircBox cbox(c0);
    Op_ptr op0 = std::make_shared<CircBox>(cbox);
    ctrl_op_map_t op_map = {
        {{1, 1}, op0},
        {{0, 1}, get_op_ptr(OpType::CX)},
        {{1, 0}, get_op_ptr(OpType::TK2, std::vector<Expr>{0.2, 0.4, 0.4})}};
    REQUIRE_THROWS_MATCHES(
        MultiplexorBox(op_map), BadOpType,
        MessageContains("Quantum control of classical wires not supported"));
  }
  GIVEN("Bitstrings are too long") {
    std::vector<bool> bits(33);
    ctrl_op_map_t op_map = {{bits, get_op_ptr(OpType::Rx, 1.4)}};
    REQUIRE_THROWS_MATCHES(
        MultiplexedRotationBox(op_map), std::invalid_argument,
        MessageContains("Bitstrings longer than 32 are not supported"));
  }
  GIVEN("Unmatched bitstrings") {
    ctrl_op_map_t op_map = {
        {{0, 1}, get_op_ptr(OpType::H)}, {{1}, get_op_ptr(OpType::X)}};
    REQUIRE_THROWS_MATCHES(
        MultiplexorBox(op_map), std::invalid_argument,
        MessageContains("Bitstrings must have the same width"));
  }
  GIVEN("Unmatched op sizes") {
    ctrl_op_map_t op_map = {
        {{0, 1}, get_op_ptr(OpType::H)}, {{1, 0}, get_op_ptr(OpType::CX)}};
    REQUIRE_THROWS_MATCHES(
        MultiplexorBox(op_map), std::invalid_argument,
        MessageContains("Ops must have the same width"));
  }
  GIVEN("Mixed rotation axis") {
    ctrl_op_map_t op_map = {
        {{1}, get_op_ptr(OpType::Rz, 0.3)}, {{0}, get_op_ptr(OpType::Rx, 1.4)}};
    REQUIRE_THROWS_MATCHES(
        MultiplexedRotationBox(op_map), std::invalid_argument,
        MessageContains("Ops must have the same rotation type"));
  }
  GIVEN("Non-rotation type") {
    ctrl_op_map_t op_map = {{{1}, get_op_ptr(OpType::H)}};
    REQUIRE_THROWS_MATCHES(
        MultiplexedRotationBox(op_map), BadOpType,
        MessageContains("Ops must be either Rx, Ry, or Rz"));
  }
  GIVEN("MultiplexedU2Box unsupported gate") {
    ctrl_op_map_t op_map = {
        {{0, 1}, get_op_ptr(OpType::H)}, {{1, 0}, get_op_ptr(OpType::CX)}};
    REQUIRE_THROWS_MATCHES(
        MultiplexedU2Box(op_map), BadOpType,
        MessageContains(
            "Ops must be single-qubit unitary gate types or Unitary1qBox"));
  }
  GIVEN("Decompose symbolic MultiplexedU2Box") {
    Sym a = SymTable::fresh_symbol("a");
    Expr expr_a(a);
    ctrl_op_map_t op_map = {{{1}, get_op_ptr(OpType::Ry, expr_a)}};
    MultiplexedU2Box u2box(op_map);
    REQUIRE_THROWS_MATCHES(
        u2box.to_circuit(), Unsupported,
        MessageContains("Can't decompose symbolic MultiplexedU2Box"));
  }
}

TEMPLATE_TEST_CASE(
    "Auxiliary methods", "[boxes]", MultiplexorBox, MultiplexedRotationBox,
    MultiplexedU2Box) {
  GIVEN("symbol_substitution") {
    Sym a = SymTable::fresh_symbol("a");
    Expr expr_a(a);
    ctrl_op_map_t op_map = {{{0}, get_op_ptr(OpType::Rz, expr_a)}};
    ctrl_op_map_t num_op_map = {{{0}, get_op_ptr(OpType::Rz, 1.34)}};
    TestType multiplexor(op_map);
    SymEngine::map_basic_basic smap;
    smap[a] = Expr(1.34);
    const TestType new_box =
        static_cast<const TestType &>(*multiplexor.symbol_substitution(smap));

    std::shared_ptr<Circuit> c = new_box.to_circuit();
    REQUIRE(check_multiplexor(num_op_map, *c));
  }
  GIVEN("free_symbols") {
    Sym a = SymTable::fresh_symbol("a");
    Sym b = SymTable::fresh_symbol("b");
    Expr expr_a(a);
    Expr expr_b(b);
    ctrl_op_map_t op_map = {
        {{0, 1}, get_op_ptr(OpType::Rz, expr_a)},
        {{1, 1}, get_op_ptr(OpType::Rz, expr_b)},
        {{1, 0}, get_op_ptr(OpType::Rz, expr_a)}};
    TestType multiplexor(op_map);
    const SymSet symbols = multiplexor.free_symbols();
    REQUIRE(symbols.size() == 2);
    REQUIRE(symbols.find(a) != symbols.end());
    REQUIRE(symbols.find(b) != symbols.end());
  }
  GIVEN("Rotation Dagger & transpose") {
    ctrl_op_map_t op_map = {
        {{0, 1}, get_op_ptr(OpType::Rz, 3.7)},
        {{1, 1}, get_op_ptr(OpType::Rz, 1)},
        {{1, 0}, get_op_ptr(OpType::Rz, 2.5)}};
    TestType multiplexor(op_map);
    // Test dagger
    const TestType dag_box =
        static_cast<const TestType &>(*multiplexor.dagger());
    std::shared_ptr<Circuit> c = dag_box.to_circuit();
    REQUIRE(check_multiplexor(op_map, c->dagger()));
    // Test transpose
    const TestType transpose_box =
        static_cast<const TestType &>(*multiplexor.transpose());
    std::shared_ptr<Circuit> d = transpose_box.to_circuit();
    REQUIRE(check_multiplexor(op_map, d->transpose()));
  }
}

SCENARIO("MultiplexorBox Dagger & transpose", "[boxes]") {
  ctrl_op_map_t op_map = {
      {{1, 1}, get_op_ptr(OpType::TK2, std::vector<Expr>{0.3, 1.8, 3.4})},
      {{0, 1}, get_op_ptr(OpType::CX)},
      {{1, 0}, get_op_ptr(OpType::TK2, std::vector<Expr>{0.2, 0.4, 0.4})}};
  MultiplexorBox multiplexor(op_map);
  // Test dagger
  const MultiplexorBox dag_box =
      static_cast<const MultiplexorBox &>(*multiplexor.dagger());
  std::shared_ptr<Circuit> c = dag_box.to_circuit();
  REQUIRE(check_multiplexor(op_map, c->dagger()));
  // Test transpose
  const MultiplexorBox transpose_box =
      static_cast<const MultiplexorBox &>(*multiplexor.transpose());
  std::shared_ptr<Circuit> d = transpose_box.to_circuit();
  REQUIRE(check_multiplexor(op_map, d->transpose()));
}

}  // namespace test_Multiplexor
}  // namespace tket
