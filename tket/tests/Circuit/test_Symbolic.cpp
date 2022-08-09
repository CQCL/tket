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

#include <Circuit/Circuit.hpp>
#include <Transformations/BasicOptimisation.hpp>
#include <Transformations/CliffordOptimisation.hpp>
#include <Transformations/OptimisationPass.hpp>
#include <Transformations/Transform.hpp>
#include <catch2/catch_test_macros.hpp>
#include <vector>

#include "OpType/OpType.hpp"
#include "Simulation/CircuitSimulator.hpp"
#include "symengine/eval_double.h"

namespace tket {
namespace test_Symbolic {

// Check (by substituting a selection of values) equivalence of two single-qubit
// circuits containing (at most) a single symbol "a".
static void check_equiv(const Circuit &circ, const Circuit &circ1) {
  static const std::vector<double> as = {0.,  0.4, 0.8, 1.2, 1.6, 2.0,
                                         2.4, 2.8, 3.2, 3.6, 4.0};
  Sym asym = SymEngine::symbol("a");
  for (double a : as) {
    INFO("circ:\n" << circ << "circ1:\n" << circ1 << "a = " << a);
    symbol_map_t smap = {{asym, a}};
    Circuit c = circ;
    c.symbol_substitution(smap);
    Eigen::MatrixXcd u = tket_sim::get_unitary(c);
    Circuit c1 = circ1;
    c1.symbol_substitution(smap);
    Eigen::MatrixXcd u1 = tket_sim::get_unitary(c1);
    CHECK(u.isApprox(u1));
  }
}

SCENARIO("Symbolic squashing, correctness") {
  Sym asym = SymEngine::symbol("a");
  Expr alpha(asym);
  Sym bsym = SymEngine::symbol("b");
  Expr beta(bsym);

  GIVEN("squash_1qb_to_pqp") {
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::Ry, 0.5, {0});
    circ.add_op<unsigned>(OpType::Rz, alpha, {0});
    circ.add_op<unsigned>(OpType::Ry, 0.5, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.5, {0});
    circ.add_op<unsigned>(OpType::Ry, 0.5, {0});
    circ.add_op<unsigned>(OpType::Rz, 1, {0});
    circ.add_op<unsigned>(OpType::Ry, 0.5, {0});

    Circuit circ1 = circ;
    Transforms::squash_1qb_to_pqp(OpType::Ry, OpType::Rz, true).apply(circ1);
    check_equiv(circ, circ1);
  }

  GIVEN("singleq_clifford_sweep (1)") {
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::U3, {2 * alpha, 0, 1.5}, {0});
    circ.add_op<unsigned>(OpType::Z, {0});
    circ.add_op<unsigned>(OpType::X, {0});

    Circuit circ1 = circ;
    Transforms::singleq_clifford_sweep().apply(circ1);
    check_equiv(circ, circ1);
  }

  GIVEN("singleq_clifford_sweep (2)") {
    Circuit circ(3);

    circ.add_op<unsigned>(OpType::U3, {alpha, 0, 0.5}, {2});
    circ.add_op<unsigned>(OpType::Vdg, {0});
    circ.add_op<unsigned>(OpType::Sdg, {2});
    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::S, {2});
    circ.add_op<unsigned>(OpType::S, {0});
    circ.add_op<unsigned>(OpType::V, {2});
    circ.add_op<unsigned>(OpType::V, {0});
    circ.add_op<unsigned>(OpType::U3, {0.5, 0, 0}, {2});
    circ.add_op<unsigned>(OpType::Rz, 0.5, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 2});
    circ.add_op<unsigned>(OpType::U3, {0.5, 1.5, 1}, {2});
    circ.add_op<unsigned>(OpType::Sdg, {2});
    circ.add_op<unsigned>(OpType::X, {2});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::V, {1});
    circ.add_op<unsigned>(OpType::Z, {2});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::X, {2});
    circ.add_op<unsigned>(OpType::Sdg, {0});
    circ.add_op<unsigned>(OpType::S, {1});
    circ.add_op<unsigned>(OpType::V, {2});
    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::Sdg, {1});
    circ.add_op<unsigned>(OpType::Vdg, {2});
    circ.add_op<unsigned>(OpType::V, {0});
    circ.add_op<unsigned>(OpType::S, {1});
    circ.add_op<unsigned>(OpType::Vdg, {2});
    circ.add_op<unsigned>(OpType::S, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::S, {2});
    circ.add_op<unsigned>(OpType::V, {2});
    circ.add_op<unsigned>(OpType::U3, {0.5, 0, 0}, {1});
    circ.add_op<unsigned>(OpType::Z, {2});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::Z, {2});
    circ.add_op<unsigned>(OpType::Rz, 0.5, {2});
    circ.add_op<unsigned>(OpType::CX, {2, 1});
    circ.add_op<unsigned>(OpType::U3, {0.5, 1.5, 1}, {1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::Sdg, {2});
    circ.add_op<unsigned>(OpType::Z, {1});
    circ.add_op<unsigned>(OpType::Vdg, {2});
    circ.add_op<unsigned>(OpType::Vdg, {1});
    circ.add_op<unsigned>(OpType::S, {2});
    circ.add_op<unsigned>(OpType::Sdg, {1});
    circ.add_op<unsigned>(OpType::TK1, {1, 0.5, 3}, {2});
    circ.add_op<unsigned>(OpType::X, {1});
    circ.add_op<unsigned>(OpType::S, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    circ.add_op<unsigned>(OpType::Z, {0});

    Circuit circ1 = circ;
    Transforms::singleq_clifford_sweep().apply(circ1);
    check_equiv(circ, circ1);
  }

  GIVEN("Edge case where symengine atan2 returns nan (1)") {
    // https://github.com/CQCL/tket/issues/304
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::Rx, {alpha}, {0});
    circ.add_op<unsigned>(OpType::Ry, {beta}, {0});
    Transforms::synthesise_tket().apply(circ);
    symbol_map_t smap = {{asym, 0}, {bsym, 0}};
    circ.symbol_substitution(smap);
    std::vector<Command> cmds = circ.get_commands();
    CHECK(cmds.size() == 1);
    Op_ptr op = cmds[0].get_op_ptr();
    CHECK(op->get_type() == OpType::TK1);
    std::vector<Expr> params = op->get_params();
    CHECK_NOTHROW(SymEngine::eval_double(params[0]));
    CHECK_NOTHROW(SymEngine::eval_double(params[1]));
    CHECK_NOTHROW(SymEngine::eval_double(params[2]));
    CHECK(approx_0(params[1]));
    CHECK(approx_0(params[0] + params[2]));
  }

  GIVEN("Edge case where symengine atan2 returns nan (2)") {
    // https://github.com/CQCL/tket/issues/304
    // and
    // https://github.com/symengine/symengine/issues/1875
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::Rx, {alpha}, {0});
    circ.add_op<unsigned>(OpType::Ry, {alpha}, {0});
    Transforms::synthesise_tket().apply(circ);
    symbol_map_t smap = {{asym, 0}};
    circ.symbol_substitution(smap);
    std::vector<Command> cmds = circ.get_commands();
    CHECK(cmds.size() == 1);
    Op_ptr op = cmds[0].get_op_ptr();
    CHECK(op->get_type() == OpType::TK1);
    std::vector<Expr> params = op->get_params();
    CHECK_NOTHROW(SymEngine::eval_double(params[0]));
    CHECK_NOTHROW(SymEngine::eval_double(params[1]));
    CHECK_NOTHROW(SymEngine::eval_double(params[2]));
    CHECK(approx_0(params[1]));
    CHECK(approx_0(params[0] + params[2]));
  }
}

}  // namespace test_Symbolic
}  // namespace tket
