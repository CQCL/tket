// Copyright 2019-2024 Cambridge Quantum Computing
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
#include <cstdint>
#include <memory>
#include <nlohmann/json_fwd.hpp>
#include <sstream>

#include "tket/Circuit/Circuit.hpp"
#include "tket/Circuit/Command.hpp"
#include "tket/OpType/EdgeType.hpp"
#include "tket/Ops/ClExpr.hpp"
#include "tket/Ops/OpPtr.hpp"
#include "tket/Utils/UnitID.hpp"

namespace tket {

SCENARIO("Circuit containing a ClExprOp") {
  GIVEN("A simple classical expression") {
    // AND of two bits:
    ClExpr expr(ClOp::BitAnd, {ClBitVar{0}, ClBitVar{1}});
    // First two bits are inputs; last bit is output:
    WiredClExpr wexpr(expr, {{0, 0}, {1, 1}}, {}, {2});
    Op_ptr op = std::make_shared<const ClExprOp>(wexpr);
    REQUIRE(op->get_signature() == op_signature_t(3, EdgeType::Classical));
    Circuit circ(0, 3);
    circ.add_op<unsigned>(op, {0, 1, 2});
    std::vector<Command> cmds = circ.get_commands();
    REQUIRE(cmds.size() == 1);
  }
  GIVEN("A complicated classical expression") {
    // d[0,1,2] <-- (a[2,1,0] + b[2,3,4]) / (c[1,0,3] * d[0,1,2])
    ClExpr numer(ClOp::RegAdd, {ClRegVar{0}, ClRegVar{1}});
    ClExpr denom(ClOp::RegMul, {ClRegVar{2}, ClRegVar{3}});
    ClExpr expr(ClOp::RegDiv, {numer, denom});
    std::vector<unsigned> a_pos{0, 3, 4};
    std::vector<unsigned> b_pos{1, 11, 5};
    std::vector<unsigned> c_pos{10, 2, 7};
    std::vector<unsigned> d_pos{8, 9, 6};
    WiredClExpr wexpr(
        expr, {}, {{0, a_pos}, {1, b_pos}, {2, c_pos}, {3, d_pos}}, d_pos);
    std::vector<unsigned> e_pos{0, 1, 2};
    REQUIRE_THROWS_AS(
        WiredClExpr(
            expr, {}, {{0, a_pos}, {1, b_pos}, {2, e_pos}, {3, d_pos}}, d_pos),
        ClExprWiringError);
    Op_ptr op = std::make_shared<const ClExprOp>(wexpr);
    Circuit circ;
    register_t preg = circ.add_c_register("p", 6);
    register_t qreg = circ.add_c_register("q", 6);
    circ.add_op<Bit>(
        op, {Bit{"p", 2}, Bit{"q", 2}, Bit{"p", 1}, Bit{"q", 3}, Bit{"p", 0},
             Bit{"q", 4}, Bit{"p", 5}, Bit{"q", 5}, Bit{"p", 4}, Bit{"q", 0},
             Bit{"p", 3}, Bit{"q", 1}});
    std::vector<Command> cmds = circ.get_commands();
    REQUIRE(cmds.size() == 1);
  }
}

SCENARIO("Serialization and stringification") {
  GIVEN("ClOp") {
    ClOp op = ClOp::RegEq;
    std::stringstream ss;
    ss << op;
    REQUIRE(ss.str() == "eq");
    nlohmann::json j = op;
    ClOp op1 = j.get<ClOp>();
    REQUIRE(op1 == op);
  }
  GIVEN("All ClOps") {
    std::stringstream ss;
    ss << ClOp::INVALID << " " << ClOp::BitAnd << " " << ClOp::BitOr << " "
       << ClOp::BitXor << " " << ClOp::BitEq << " " << ClOp::BitNeq << " "
       << ClOp::BitNot << " " << ClOp::BitZero << " " << ClOp::BitOne << " "
       << ClOp::RegAnd << " " << ClOp::RegOr << " " << ClOp::RegXor << " "
       << ClOp::RegEq << " " << ClOp::RegNeq << " " << ClOp::RegNot << " "
       << ClOp::RegZero << " " << ClOp::RegOne << " " << ClOp::RegLt << " "
       << ClOp::RegGt << " " << ClOp::RegLeq << " " << ClOp::RegGeq << " "
       << ClOp::RegAdd << " " << ClOp::RegSub << " " << ClOp::RegMul << " "
       << ClOp::RegDiv << " " << ClOp::RegPow << " " << ClOp::RegLsh << " "
       << ClOp::RegRsh << " " << ClOp::RegNeg;
    REQUIRE(
        ss.str() ==
        "INVALID and or xor eq neq not zero one and or xor eq neq not zero one "
        "lt gt leq geq add sub mul div pow lsh rsh neg");
  }
  GIVEN("ClBitVar") {
    ClBitVar var{3};
    std::stringstream ss;
    ss << var;
    REQUIRE(ss.str() == "b3");
    nlohmann::json j = var;
    ClBitVar var1 = j.get<ClBitVar>();
    REQUIRE(var1 == var);
  }
  GIVEN("ClRegVar") {
    ClRegVar var{4};
    std::stringstream ss;
    ss << var;
    REQUIRE(ss.str() == "r4");
    nlohmann::json j = var;
    ClRegVar var1 = j.get<ClRegVar>();
    REQUIRE(var1 == var);
  }
  GIVEN("ClExprVar") {
    ClExprVar var_bit = ClBitVar{3};
    ClExprVar var_reg = ClRegVar{4};
    std::stringstream ss;
    ss << var_bit << ", " << var_reg;
    REQUIRE(ss.str() == "b3, r4");
    nlohmann::json j_bit = var_bit;
    nlohmann::json j_reg = var_reg;
    ClExprVar var_bit1 = j_bit.get<ClExprVar>();
    ClExprVar var_reg1 = j_reg.get<ClExprVar>();
    REQUIRE(var_bit1 == var_bit);
    REQUIRE(var_reg1 == var_reg);
  }
  GIVEN("ClExprTerm") {
    ClExprTerm term_int = uint64_t{7};
    ClExprTerm term_var = ClRegVar{5};
    std::stringstream ss;
    ss << term_int << ", " << term_var;
    REQUIRE(ss.str() == "7, r5");
    nlohmann::json j_int = term_int;
    nlohmann::json j_var = term_var;
    ClExprTerm term_int1 = j_int.get<ClExprTerm>();
    ClExprTerm term_var1 = j_var.get<ClExprTerm>();
    REQUIRE(term_int1 == term_int);
    REQUIRE(term_var1 == term_var);
  }
  GIVEN("Vector of ClExprArg (1)") {
    std::vector<ClExprArg> args{ClRegVar{2}, uint64_t{3}};
    nlohmann::json j = args;
    std::vector<ClExprArg> args1 = j.get<std::vector<ClExprArg>>();
    REQUIRE(args == args1);
  }
  GIVEN("ClExpr (1)") {
    // r0 + 7
    ClExpr expr(ClOp::RegAdd, {ClRegVar{0}, uint64_t{7}});
    std::stringstream ss;
    ss << expr;
    REQUIRE(ss.str() == "add(r0, 7)");
    nlohmann::json j = expr;
    ClExpr expr1 = j.get<ClExpr>();
    REQUIRE(expr1 == expr);
  }
  GIVEN("Vector of ClExprArg (2)") {
    ClExpr expr(ClOp::RegAdd, {ClRegVar{0}, uint64_t{8}});
    std::vector<ClExprArg> args{expr};
    nlohmann::json j = args;
    std::vector<ClExprArg> args1 = j.get<std::vector<ClExprArg>>();
    REQUIRE(args == args1);
  }
  GIVEN("ClExpr (2)") {
    // (r0 + r1) / (r2 * 3)
    ClExpr numer(ClOp::RegAdd, {ClRegVar{0}, ClRegVar{1}});
    ClExpr denom(ClOp::RegMul, {ClRegVar{2}, uint64_t{3}});
    ClExpr expr(ClOp::RegDiv, {numer, denom});
    std::stringstream ss;
    ss << expr;
    REQUIRE(ss.str() == "div(add(r0, r1), mul(r2, 3))");
    nlohmann::json j = expr;
    ClExpr expr1 = j.get<ClExpr>();
    REQUIRE(expr1 == expr);
  }
  GIVEN("WiredClExpr") {
    ClExpr numer(ClOp::RegAdd, {ClRegVar{0}, ClRegVar{1}});
    ClExpr denom(ClOp::RegMul, {ClRegVar{2}, ClRegVar{3}});
    ClExpr expr(ClOp::RegDiv, {numer, denom});
    std::vector<unsigned> a_pos{0, 3, 4};
    std::vector<unsigned> b_pos{1, 11, 5};
    std::vector<unsigned> c_pos{10, 2, 7};
    std::vector<unsigned> d_pos{8, 9, 6};
    WiredClExpr wexpr(
        expr, {}, {{0, a_pos}, {1, b_pos}, {2, c_pos}, {3, d_pos}}, d_pos);
    std::stringstream ss;
    ss << wexpr;
    REQUIRE(
        ss.str() ==
        "div(add(r0, r1), mul(r2, r3)) [r0:(0,3,4), r1:(1,11,5), r2:(10,2,7), "
        "r3:(8,9,6) --> (8,9,6)]");
    nlohmann::json j = wexpr;
    WiredClExpr wexpr1 = j.get<WiredClExpr>();
    REQUIRE(wexpr1 == wexpr);
  }
  GIVEN("ClExprOp") {
    ClExpr numer(ClOp::RegAdd, {ClRegVar{0}, ClRegVar{1}});
    ClExpr denom(ClOp::RegMul, {ClRegVar{2}, ClRegVar{3}});
    ClExpr expr(ClOp::RegDiv, {numer, denom});
    std::vector<unsigned> a_pos{0, 3, 4};
    std::vector<unsigned> b_pos{1, 11, 5};
    std::vector<unsigned> c_pos{10, 2, 7};
    std::vector<unsigned> d_pos{8, 9, 6};
    WiredClExpr wexpr(
        expr, {}, {{0, a_pos}, {1, b_pos}, {2, c_pos}, {3, d_pos}}, d_pos);
    Op_ptr op = std::make_shared<ClExprOp>(wexpr);
    nlohmann::json j = op;
    Op_ptr op1 = j.get<Op_ptr>();
    const ClExprOp& exprop = static_cast<const ClExprOp&>(*op1);
    REQUIRE(exprop.get_wired_expr() == wexpr);
    Op_ptr op2 = op->symbol_substitution({});
    REQUIRE(op2->free_symbols().empty());
  }
}

}  // namespace tket
