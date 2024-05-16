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

#include "tket/Circuit/CircPool.hpp"
#include "tket/Circuit/CircUtils.hpp"
#include "tket/Circuit/Simulation/CircuitSimulator.hpp"
#include "tket/Predicates/Predicates.hpp"

namespace tket {
namespace test_CircPool {

SCENARIO("Simple CircPool identities") {
  Circuit orig, res;

  GIVEN("tk1_to_tk1") {
    orig = Circuit(1);
    orig.add_op<unsigned>(OpType::TK1, {0.2, 0.3, 0.4}, {0});
    res = CircPool::tk1_to_tk1(0.2, 0.3, 0.4);
  }
  GIVEN("CCX") {
    orig = Circuit(3);
    orig.add_op<unsigned>(OpType::CCX, {0, 1, 2});
    res = CircPool::CCX();
  }
  GIVEN("BRIDGE") {
    orig = Circuit(3);
    orig.add_op<unsigned>(OpType::BRIDGE, {0, 1, 2});
    res = CircPool::BRIDGE();
  }
  GIVEN("H_CZ_H") {
    orig = Circuit(2);
    orig.add_op<unsigned>(OpType::CX, {0, 1});
    res = CircPool::H_CZ_H();
  }
  GIVEN("CX_using_AAMS") {
    orig = Circuit(2);
    orig.add_op<unsigned>(OpType::CX, {0, 1});
    res = CircPool::CX_using_AAMS();
  }

  auto u_orig = tket_sim::get_unitary(orig);
  auto u_res = tket_sim::get_unitary(res);
  REQUIRE(u_res.isApprox(u_orig));
}

SCENARIO("TK2_using_normalised_TK2") {
  Expr e1, e2, e3;
  Sym asym = SymEngine::symbol("a");
  Expr a(asym);
  Sym bsym = SymEngine::symbol("b");
  Expr b(bsym);
  Sym csym = SymEngine::symbol("c");
  Expr c(csym);

  GIVEN("Normalised concrete angles (1)") {
    e1 = .3;
    e2 = .1;
    e3 = .05;
  }
  GIVEN("Normalised concrete angles (2)") {
    e1 = 0.32;
    e2 = 0.31;
    e3 = -0.3;
  }
  GIVEN("Not normalised concrete angles (1)") {
    e1 = .3;
    e2 = .4;
    e3 = .45;
  }
  GIVEN("Not normalised concrete angles (2)") {
    e1 = .3;
    e2 = 1.4;
    e3 = .489;
  }
  GIVEN("Not normalised concrete angles (3)") {
    e1 = 2.3;
    e2 = 3.4;
    e3 = .489;
  }
  GIVEN("Not normalised concrete angles (4)") {
    e1 = .3;
    e2 = -.2;
    e3 = .1;
  }
  GIVEN("Not normalised concrete angles (5)") {
    e1 = -.3;
    e2 = -.2;
    e3 = .1;
  }
  GIVEN("Not normalised concrete angles (6)") {
    e1 = .3;
    e2 = .2;
    e3 = -.3;
  }
  GIVEN("Not normalised concrete angles (7)") {
    e1 = 0;
    e2 = 0;
    e3 = -1.2;
  }
  GIVEN("Not normalised concrete angles (8)") {
    e1 = 0.1;
    e2 = 0.3;
    e3 = 0.2;
  }
  GIVEN("Symbolic angles (1)") {
    e1 = a;
    e2 = 3.4;
    e3 = .489;
  }
  GIVEN("Symbolic angles (2)") {
    e1 = a;
    e2 = b;
    e3 = 2.42;
  }
  GIVEN("Symbolic angles (3)") {
    e1 = 2.3;
    e2 = b;
    e3 = 1.489;
  }
  GIVEN("Symbolic angles (4)") {
    e1 = 2.3;
    e2 = 123.08174;
    e3 = c;
  }
  GIVEN("Symbolic angles (5)") {
    e1 = a;
    e2 = 123.08174;
    e3 = c;
  }
  GIVEN("Symbolic angles (6)") {
    e1 = 0.10012;
    e2 = b;
    e3 = c;
  }

  Circuit orig(2);
  orig.add_op<unsigned>(OpType::TK2, {e1, e2, e3}, {0, 1});
  Circuit res = CircPool::TK2_using_normalised_TK2(e1, e2, e3);

  REQUIRE(NormalisedTK2Predicate().verify(res));

  // check unitary identity
  auto symset = orig.free_symbols();
  std::vector<Sym> symbols(symset.begin(), symset.end());
  if (symbols.empty()) {
    auto u_orig = tket_sim::get_unitary(orig);
    auto u_res = tket_sim::get_unitary(res);
    REQUIRE(u_res.isApprox(u_orig));
  } else {
    std::vector<double> rands{0.1231, 2.3124, 34.23, 2.23, 3.15, 1.2, 0.93};
    // substitute random values for symbolics and check equality
    unsigned i = 0;
    while (i + symbols.size() <= rands.size()) {
      symbol_map_t symmap;
      for (unsigned j = 0; j < symbols.size(); ++j) {
        symmap[symbols[j]] = rands[i + j];
      }
      Circuit res_sub = res;
      Circuit orig_sub = orig;
      res_sub.symbol_substitution(symmap);
      orig_sub.symbol_substitution(symmap);
      auto u_orig = tket_sim::get_unitary(orig_sub);
      auto u_res = tket_sim::get_unitary(res_sub);
      REQUIRE(u_res.isApprox(u_orig));
      ++i;
    }
  }
}

SCENARIO("TK2_using_ZZMax") {
  double a, b, c;
  GIVEN("General angles") {
    a = 1.2;
    b = 2.3;
    c = 3.4;
  }
  GIVEN("c = 0") {
    a = 0.3;
    b = 1.2;
    c = 0.;
  }
  GIVEN("b = c = 0") {
    a = -1.9;
    b = c = 0.;
  }
  GIVEN("a = b = c = 0") { a = b = c = 0.; }
  Circuit circ(2);
  Circuit orig(2);
  orig.add_op<unsigned>(OpType::TK2, {a, b, c}, {0, 1});
  Circuit res = CircPool::TK2_using_ZZMax(a, b, c);
  Eigen::Matrix4cd u_orig = tket_sim::get_unitary(orig);
  Eigen::Matrix4cd u_res = tket_sim::get_unitary(res);
  CHECK(u_res.isApprox(u_orig));
}

SCENARIO("Test remove_noops") {
  GIVEN("A circuit with noops") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::U1, 0., {0});
    circ.add_op<unsigned>(OpType::Rx, 0., {0});
    circ.add_op<unsigned>(OpType::U3, {0., 0., 0.}, {0});
    circ.add_op<unsigned>(OpType::TK1, {0., 0., 0.}, {0});
    circ.add_op<unsigned>(OpType::TK2, {0., 0., 0.}, {0, 1});

    circ.remove_noops();

    REQUIRE(circ == Circuit(2));
  }
  GIVEN("A circuit with mix of ops and noops") {
    Circuit circ(2), circ2(2);
    circ.add_op<unsigned>(OpType::U1, 0., {0});
    circ.add_op<unsigned>(OpType::Rx, 0., {0});
    circ.add_op<unsigned>(OpType::U2, {0., 1.2}, {0});
    circ2.add_op<unsigned>(OpType::U2, {0., 1.2}, {0});
    circ.add_op<unsigned>(OpType::U3, {0., 0., 0.}, {0});
    circ.add_op<unsigned>(OpType::U3, {0.1, 0.2, 1.2}, {0});
    circ2.add_op<unsigned>(OpType::U3, {0.1, 0.2, 1.2}, {0});
    circ.add_op<unsigned>(OpType::TK1, {0., 0., 0.}, {0});
    circ.add_op<unsigned>(OpType::TK2, {0., 0., 0.}, {0, 1});
    circ.add_op<unsigned>(OpType::TK2, {0.1, 0.3, 2.1}, {0, 1});
    circ2.add_op<unsigned>(OpType::TK2, {0.1, 0.3, 2.1}, {0, 1});

    circ.remove_noops();

    REQUIRE(circ == circ2);
  }
}

SCENARIO(
    "Rx_using_GPI, Ry_using_GPI, Rz_using_GPI, XXPhase_using_AAMS,"
    "YYPhase_using_AAMS, ZZPhase_using_AAMS") {
  Expr e1;
  Sym asym = SymEngine::symbol("a");
  Expr a(asym);

  GIVEN("Normalised concrete angles (1)") { e1 = .3; }
  GIVEN("Normalised concrete angles (2)") { e1 = -0.32; }
  GIVEN("Not normalised concrete angles (1)") { e1 = 1.4; }
  GIVEN("Not normalised concrete angles (2)") { e1 = -5.7; }
  GIVEN("Symbolic angles (1)") { e1 = a; }
 
  Circuit rx_orig(1);
  rx_orig.add_op<unsigned>(OpType::Rx, {e1}, {0});
  Circuit ry_orig(1);
  ry_orig.add_op<unsigned>(OpType::Ry, {e1}, {0});
  Circuit rz_orig(1);
  rz_orig.add_op<unsigned>(OpType::Rz, {e1}, {0});
  Circuit xxphase_orig(2);
  xxphase_orig.add_op<unsigned>(OpType::XXPhase, {e1}, {0, 1});
  Circuit yyphase_orig(2);
  yyphase_orig.add_op<unsigned>(OpType::YYPhase, {e1}, {0, 1});
  Circuit zzphase_orig(2);
  zzphase_orig.add_op<unsigned>(OpType::ZZPhase, {e1}, {0, 1});
  std::vector<Circuit> circuits_orig = {
    rx_orig, ry_orig, rz_orig, xxphase_orig, yyphase_orig, zzphase_orig};
 
  Circuit rx_res = CircPool::Rx_using_GPI(e1);
  Circuit ry_res = CircPool::Ry_using_GPI(e1);
  Circuit rz_res = CircPool::Rz_using_GPI(e1);
  Circuit xxphase_res = CircPool::XXPhase_using_AAMS(e1);
  Circuit yyphase_res = CircPool::YYPhase_using_AAMS(e1);
  Circuit zzphase_res = CircPool::ZZPhase_using_AAMS(e1);
  std::vector<Circuit> circuits_res = {rx_res,      ry_res,      rz_res,
	                                     xxphase_res, yyphase_res, zzphase_res};

  // check unitary identity
  auto symset = circuits_orig[0].free_symbols();
  std::vector<Sym> symbols(symset.begin(), symset.end());
  if (symbols.empty()) {
    for (unsigned k = 0; k < 6; k++) {
      auto u_orig = tket_sim::get_unitary(circuits_orig[k]);
      auto u_res = tket_sim::get_unitary(circuits_res[k]);
      REQUIRE(u_res.isApprox(u_orig));
    }
  } else {
    std::vector<double> rands{0.1231, 2.3124, 34.23, 2.23, 3.15, 1.2, 0.93};
    // substitute random values for symbolics and check equality
    unsigned i = 0;
    while (i + symbols.size() <= rands.size()) {
      symbol_map_t symmap;
      for (unsigned j = 0; j < symbols.size(); ++j) {
        symmap[symbols[j]] = rands[i + j];
      }
      for (unsigned k = 0; k < 6; ++k) {
        Circuit orig_sub = circuits_orig[k];
        orig_sub.symbol_substitution(symmap);
        auto u_orig = tket_sim::get_unitary(orig_sub);
        Circuit res_sub = circuits_res[k];
        res_sub.symbol_substitution(symmap);
        auto u_res = tket_sim::get_unitary(res_sub);
        REQUIRE(u_res.isApprox(u_orig));
      }
      ++i;
    }
  }
}

SCENARIO("TK1_using_GPI, TK2_using_AAMS") {
  Expr e1, e2, e3;
  Sym asym = SymEngine::symbol("a");
  Expr a(asym);
  Sym bsym = SymEngine::symbol("b");
  Expr b(bsym);
  Sym csym = SymEngine::symbol("c");
  Expr c(csym);

  GIVEN("Normalised concrete angles (1)") {
    e1 = .3;
    e2 = .1;
    e3 = .05;
  }
  GIVEN("Normalised concrete angles (2)") {
    e1 = 0.32;
    e2 = 0.31;
    e3 = -0.3;
  }
  GIVEN("Not normalised concrete angles (1)") {
    e1 = .3;
    e2 = .4;
    e3 = .45;
  }
  GIVEN("Not normalised concrete angles (2)") {
    e1 = .3;
    e2 = 1.4;
    e3 = .489;
  }
  GIVEN("Not normalised concrete angles (3)") {
    e1 = 2.3;
    e2 = 3.4;
    e3 = .489;
  }
  GIVEN("Not normalised concrete angles (4)") {
    e1 = .3;
    e2 = -.2;
    e3 = .1;
  }
  GIVEN("Not normalised concrete angles (5)") {
    e1 = -.3;
    e2 = -.2;
    e3 = .1;
  }
  GIVEN("Not normalised concrete angles (6)") {
    e1 = .3;
    e2 = .2;
    e3 = -.3;
  }
  GIVEN("Not normalised concrete angles (7)") {
    e1 = 0;
    e2 = 0;
    e3 = -1.2;
  }
  GIVEN("Not normalised concrete angles (8)") {
    e1 = 0.1;
    e2 = 0.3;
    e3 = 0.2;
  }
  GIVEN("Symbolic angles (1)") {
    e1 = a;
    e2 = 3.4;
    e3 = .489;
  }
  GIVEN("Symbolic angles (2)") {
    e1 = a;
    e2 = b;
    e3 = 2.42;
  }
  GIVEN("Symbolic angles (3)") {
    e1 = 2.3;
    e2 = b;
    e3 = 1.489;
  }
  GIVEN("Symbolic angles (4)") {
    e1 = 2.3;
    e2 = 123.08174;
    e3 = c;
  }
  GIVEN("Symbolic angles (5)") {
    e1 = a;
    e2 = 123.08174;
    e3 = c;
  }
  GIVEN("Symbolic angles (6)") {
    e1 = 0.10012;
    e2 = b;
    e3 = c;
  }

  Circuit tk1_orig(1);
  tk1_orig.add_op<unsigned>(OpType::TK1, {e1, e2, e3}, {0});
  Circuit tk1_res = CircPool::TK1_using_GPI(e1, e2, e3);

  Circuit tk2_orig(2);
  tk2_orig.add_op<unsigned>(OpType::TK2, {e1, e2, e3}, {0, 1});
  Circuit tk2_res = CircPool::TK2_using_AAMS(e1, e2, e3);

  // check unitary identity
  auto symset = tk1_orig.free_symbols();
  std::vector<Sym> symbols(symset.begin(), symset.end());
  if (symbols.empty()) {
    auto u_tk1_orig = tket_sim::get_unitary(tk1_orig);
    auto u_tk2_orig = tket_sim::get_unitary(tk2_orig);
    auto u_tk1_res = tket_sim::get_unitary(tk1_res);
    auto u_tk2_res = tket_sim::get_unitary(tk2_res);
    REQUIRE(u_tk1_res.isApprox(u_tk1_orig));
    REQUIRE(u_tk2_res.isApprox(u_tk2_orig));
  } else {
    std::vector<double> rands{0.1231, 2.3124, 34.23, 2.23, 3.15, 1.2, 0.93};
    // substitute random values for symbolics and check equality
    unsigned i = 0;
    while (i + symbols.size() <= rands.size()) {
      symbol_map_t symmap;
      for (unsigned j = 0; j < symbols.size(); ++j) {
        symmap[symbols[j]] = rands[i + j];
      }
      Circuit tk1_orig_sub = tk1_orig;
      Circuit tk2_orig_sub = tk2_orig;
      tk1_orig_sub.symbol_substitution(symmap);
      tk2_orig_sub.symbol_substitution(symmap);
      auto u_tk1_orig = tket_sim::get_unitary(tk1_orig_sub);
      auto u_tk2_orig = tket_sim::get_unitary(tk2_orig_sub);
      Circuit tk1_res_sub = tk1_res;
      Circuit tk2_res_sub = tk2_res;
      tk1_res_sub.symbol_substitution(symmap);
      tk2_res_sub.symbol_substitution(symmap);
      auto u_tk1_res = tket_sim::get_unitary(tk1_res_sub);
      auto u_tk2_res = tket_sim::get_unitary(tk2_res_sub);
      REQUIRE(u_tk1_res.isApprox(u_tk1_orig));
      REQUIRE(u_tk2_res.isApprox(u_tk2_orig));
      ++i;
    }
  }
}

}  // namespace test_CircPool
}  // namespace tket
