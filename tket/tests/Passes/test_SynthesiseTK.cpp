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

#include "../testutil.hpp"
#include "Circuit/Circuit.hpp"
#include "OpType/OpType.hpp"
#include "Predicates/CompilationUnit.hpp"
#include "Predicates/PassLibrary.hpp"
#include "Utils/Expression.hpp"

namespace tket {
namespace test_SynthesiseTK {

void check_synthesise_tk(const Circuit& c) {
  Circuit c0 = c;
  CompilationUnit cu(c0);
  SynthesiseTK()->apply(cu);
  Circuit c1 = cu.get_circ_ref();

  // Substitute "arbitrary" values for symbols in both circuits
  const SymSet symbols = c.free_symbols();
  symbol_map_t smap;
  unsigned i = 0;
  for (const Sym& s : symbols) {
    smap[s] = PI * (i + 1) / ((i + 2) * (i + 3));
    i++;
  }
  c0.symbol_substitution(smap);
  c1.symbol_substitution(smap);

  REQUIRE(test_unitary_comparison(c0, c1));

  OpTypeSet tk_types = {OpType::TK1, OpType::TK2};
  GateSetPredicate pred(tk_types);
  REQUIRE(pred.verify(c1));
}

SCENARIO("SynthesiseTK correctness") {
  GIVEN("A simple circuit") {
    Circuit c(2);
    c.add_op<unsigned>(OpType::H, {0});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    check_synthesise_tk(c);
  }
  GIVEN("A bigger circuit") {
    Circuit c(2);
    c.add_op<unsigned>(OpType::H, {0});
    c.add_op<unsigned>(OpType::CX, {0, 1});
    c.add_op<unsigned>(OpType::X, {0});
    c.add_op<unsigned>(OpType::CY, {1, 0});
    c.add_op<unsigned>(OpType::Y, {1});
    c.add_op<unsigned>(OpType::CZ, {0, 1});
    c.add_op<unsigned>(OpType::T, {0});
    c.add_op<unsigned>(OpType::SWAP, {0, 1});
    c.add_op<unsigned>(OpType::S, {1});
    c.add_op<unsigned>(OpType::CH, {1, 0});
    c.add_op<unsigned>(OpType::Sdg, {0});
    c.add_op<unsigned>(OpType::CV, {1, 0});
    c.add_op<unsigned>(OpType::Tdg, {1});
    c.add_op<unsigned>(OpType::CVdg, {0, 1});
    c.add_op<unsigned>(OpType::Vdg, {0});
    c.add_op<unsigned>(OpType::CSX, {1, 0});
    c.add_op<unsigned>(OpType::SX, {1});
    c.add_op<unsigned>(OpType::CSXdg, {0, 1});
    c.add_op<unsigned>(OpType::SXdg, {0});
    c.add_op<unsigned>(OpType::ECR, {1, 0});
    c.add_op<unsigned>(OpType::noop, {1});
    c.add_op<unsigned>(OpType::ZZMax, {0, 1});
    c.add_op<unsigned>(OpType::Sycamore, {1, 0});
    c.add_op<unsigned>(OpType::ISWAPMax, {0, 1});
    check_synthesise_tk(c);
  }
  GIVEN("A circuit with (>2)-qubit gates") {
    Circuit c(4);
    c.add_op<unsigned>(OpType::CCX, {0, 1, 2});
    c.add_op<unsigned>(OpType::CnX, {0, 1, 2, 3});
    c.add_op<unsigned>(OpType::CnRy, 0.25, {0, 1, 2, 3});
    c.add_op<unsigned>(OpType::CSWAP, {1, 2, 3});
    c.add_op<unsigned>(OpType::BRIDGE, {1, 3, 0});
    check_synthesise_tk(c);
  }
  GIVEN("A circuit with symbolic gates") {
    Sym asym = SymEngine::symbol("a");
    Expr a(asym);
    Circuit c(3);
    c.add_op<unsigned>(OpType::Rx, a, {0});
    c.add_op<unsigned>(OpType::Ry, a + 0.1, {0});
    c.add_op<unsigned>(OpType::Rz, a + 0.2, {1});
    c.add_op<unsigned>(OpType::XXPhase, a + 0.3, {0, 1});
    c.add_op<unsigned>(OpType::YYPhase, a + 0.1, {1, 0});
    c.add_op<unsigned>(OpType::ZZPhase, a + 0.2, {0, 1});
    c.add_op<unsigned>(OpType::U1, a, {1});
    c.add_op<unsigned>(OpType::U2, {a + 0.1, a + 0.2}, {0});
    c.add_op<unsigned>(OpType::U3, {a + 0.3, a + 0.4, a + 0.5}, {1});
    c.add_op<unsigned>(OpType::TK1, {a + 0.2, a + 0.3, a + 0.4}, {0});
    c.add_op<unsigned>(OpType::CRz, a, {0, 1});
    c.add_op<unsigned>(OpType::CRx, a + 0.1, {1, 0});
    c.add_op<unsigned>(OpType::CRy, a + 0.2, {0, 1});
    c.add_op<unsigned>(OpType::CU1, a + 0.3, {1, 0});
    c.add_op<unsigned>(OpType::CU3, {a + 0.1, a + 0.2, a}, {0, 1});
    c.add_op<unsigned>(OpType::ISWAP, a, {0, 1});
    c.add_op<unsigned>(OpType::PhasedISWAP, {a, a - 0.1}, {1, 0});
    c.add_op<unsigned>(OpType::XXPhase3, a, {0, 1, 2});
    c.add_op<unsigned>(OpType::PhasedX, {a, a + 0.2}, {0});
    c.add_op<unsigned>(OpType::NPhasedX, {a + 0.3, a}, {0, 1, 2});
    c.add_op<unsigned>(OpType::CnRy, a - 0.3, {1, 2, 0});
    c.add_op<unsigned>(OpType::ESWAP, a, {1, 2});
    c.add_op<unsigned>(OpType::FSim, {a + 0.1, a + 0.2}, {2, 0});
    check_synthesise_tk(c);
  }
}

}  // namespace test_SynthesiseTK
}  // namespace tket
