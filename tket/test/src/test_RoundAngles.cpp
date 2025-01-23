// Copyright Quantinuum
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

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <vector>

#include "tket/Circuit/Circuit.hpp"
#include "tket/Predicates/CompilationUnit.hpp"
#include "tket/Predicates/PassGenerators.hpp"

namespace tket {
namespace test_RoundAngles {

SCENARIO("Rounding angles") {
  GIVEN("Circuit with pi/8 angle") {
    Circuit c(2);
    c.add_op<unsigned>(OpType::H, {0});
    c.add_op<unsigned>(OpType::CRz, 0.125, {0, 1});
    CompilationUnit cu(c);
    REQUIRE_FALSE(RoundAngles(3)->apply(cu));
    REQUIRE(cu.get_circ_ref().n_gates() == 2);
    REQUIRE(RoundAngles(1)->apply(cu));
    REQUIRE(cu.get_circ_ref().n_gates() == 1);
  }
  GIVEN("Circuit with no parameters") {
    Circuit c(2);
    c.add_op<unsigned>(OpType::H, {0});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    CompilationUnit cu(c);
    REQUIRE_FALSE(RoundAngles(3)->apply(cu));
  }
  GIVEN("Request to only round to zero") {
    Circuit c(2);
    c.add_op<unsigned>(OpType::H, {0});
    c.add_op<unsigned>(OpType::CRz, 0.250001, {0, 1});
    CompilationUnit cu(c);
    REQUIRE_FALSE(RoundAngles(2, true)->apply(cu));
    REQUIRE(RoundAngles(2, false)->apply(cu));
  }
  GIVEN("Gate with a tiny angle") {
    Circuit c(2);
    c.add_op<unsigned>(OpType::H, {0});
    c.add_op<unsigned>(OpType::CRz, 0.000001, {0, 1});
    CompilationUnit cu(c);
    REQUIRE(RoundAngles(16)->apply(cu));
    const Circuit &c1 = cu.get_circ_ref();
    std::vector<Command> cmds = c1.get_commands();
    REQUIRE(cmds.size() == 1);
    REQUIRE(cmds[0].get_op_ptr()->get_type() == OpType::H);
  }
  GIVEN("Multi-parameter gate with one non-tiny angle") {
    Circuit c(2);
    c.add_op<unsigned>(OpType::TK2, {0.2, 0.001, -0.7}, {0, 1});
    CompilationUnit cu(c);
    REQUIRE(RoundAngles(4)->apply(cu));
    const Circuit &c1 = cu.get_circ_ref();
    std::vector<Command> cmds = c1.get_commands();
    REQUIRE(cmds.size() == 1);
    Op_ptr op = cmds[0].get_op_ptr();
    REQUIRE(op->get_type() == OpType::TK2);
    REQUIRE(op->get_params()[1] == 0.);
  }
  GIVEN("A nearly-Clifford circuit") {
    Circuit c(2, 2);
    c.add_op<unsigned>(OpType::H, {0});
    c.add_op<unsigned>(OpType::Ry, 0.5001, {0});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::TK1, {-0.4999, 1.5001, 0.999}, {0});
    c.add_measure(0, 0);
    c.add_measure(1, 1);
    CompilationUnit cu(c);
    REQUIRE(RoundAngles(1)->apply(cu));
    const Circuit &c1 = cu.get_circ_ref();
    std::vector<Command> cmds = c1.get_commands();
    REQUIRE(std::all_of(cmds.begin(), cmds.end(), [](const Command &cmd) {
      return cmd.get_op_ptr()->is_clifford() ||
             cmd.get_op_ptr()->get_type() == OpType::Measure;
    }));
  }
  GIVEN("An invalid precision parameter") {
    Circuit c(2);
    c.add_op<unsigned>(OpType::H, {0});
    c.add_op<unsigned>(OpType::CRz, 0.000001, {0, 1});
    CompilationUnit cu(c);
    REQUIRE_THROWS(RoundAngles(32)->apply(cu));
  }
  GIVEN("A conditional op") {
    Circuit c(3, 1);
    c.add_op<unsigned>(OpType::H, {0});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_measure(0, 0);
    c.add_conditional_gate<unsigned>(
        OpType::TK2, {0.499, 0.501, 0.5}, {1, 2}, {0}, 1);
    CompilationUnit cu(c);
    REQUIRE(RoundAngles(3)->apply(cu));
    const Circuit &c1 = cu.get_circ_ref();
    std::vector<Command> cmds = c1.get_commands();
    REQUIRE(cmds.size() == 4);
    Op_ptr op = cmds[3].get_op_ptr();
    REQUIRE(op->get_type() == OpType::Conditional);
    const Conditional &cond = static_cast<const Conditional &>(*op);
    REQUIRE(cond.get_op()->get_params() == std::vector<Expr>{0.5, 0.5, 0.5});
  }
}

}  // namespace test_RoundAngles
}  // namespace tket
