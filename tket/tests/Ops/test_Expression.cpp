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

#include <symengine/eval.h>

#include <catch2/catch_test_macros.hpp>

#include "../testutil.hpp"
#include "Utils/Expression.hpp"

namespace tket {
namespace test_Expression {

SCENARIO("Basic Expr evaluation", "[ops]") {
  GIVEN("A constant") {
    Expr e(2.5);
    REQUIRE(test_equiv_val(e, 0.5));
  }
  GIVEN("A symbol") {
    Sym s = SymEngine::symbol("a");
    Expr e = Expr(s);
    REQUIRE(!eval_expr_mod(e));
    SymEngine::map_basic_basic smap;
    smap[s] = Expr(3.4);
    Expr ee = e.subs(smap);
    REQUIRE(test_equiv_val(ee, 1.4));
  }
  GIVEN("A non-empty sum") {
    Sym s = SymEngine::symbol("b");
    Expr e = 0.2 + Expr(s) + 0.5 + Expr(s);
    REQUIRE(!eval_expr_mod(e));
    SymEngine::map_basic_basic smap;
    smap[s] = Expr(0.3);
    Expr ee = e.subs(smap);
    REQUIRE(test_equiv_val(ee, 1.3));
  }
  GIVEN("A non-empty product") {
    Sym s = SymEngine::symbol("b");
    Expr e = 0.2 * Expr(s) * 0.5 * Expr(s);
    REQUIRE(!eval_expr_mod(e));
    SymEngine::map_basic_basic smap;
    smap[s] = Expr(3.);
    Expr ee = e.subs(smap);
    REQUIRE(test_equiv_val(ee, 0.9));
  }
  GIVEN("A more complicated expression") {
    Sym s = SymEngine::symbol("d");
    Expr e = -0.3 + (3.4 * Expr(SymEngine::sin(Expr(s) - 2.3)));
    SymEngine::map_basic_basic smap;
    smap[s] = Expr(2.3);
    Expr ee = e.subs(smap);
    REQUIRE(test_equiv_val(ee, 1.7));
  }
}

SCENARIO("Expression uniqueness", "[ops]") {
  GIVEN("Two equivalent constants") {
    Expr a(0.5);
    Expr b(2 * 3. / 4 - 1);
    b = SymEngine::evalf(b, 53);
    REQUIRE(a == b);
  }
  GIVEN("Two different constants") {
    Expr a(2.);
    Expr b0(2);
    Expr b1(3.);
    REQUIRE(a != b0);
    REQUIRE(a != b1);
  }
  GIVEN("Two identical symbols") {
    Expr a("alpha");
    Expr b("alpha");
    REQUIRE(a == b);
  }
  GIVEN("Two different symbols") {
    Expr a("alpha");
    Expr b("beta");
    REQUIRE(a != b);
  }
  GIVEN("Parsed atan2") {
    Expr a("alpha");
    Expr b("beta");
    Expr at(SymEngine::atan2(a, b));
    Expr at2 = Expr("atan2(alpha, beta)");
    REQUIRE(at == at2);
  }
}

}  // namespace test_Expression
}  // namespace tket
