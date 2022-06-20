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

#include "Circuit/Circuit.hpp"
#include "OpType/OpType.hpp"
#include "Predicates/CompilationUnit.hpp"
#include "Predicates/PassLibrary.hpp"

namespace tket {
namespace test_SynthesiseTket {

SCENARIO("SynthesiseTket with conditionals") {
  Circuit c(2);
  c.add_c_register("c", 1);
  c.add_conditional_gate<unsigned>(OpType::CnRy, {0.25}, {0, 1}, {0}, 0);
  c.add_conditional_gate<unsigned>(OpType::Ry, {0.125}, {1}, {0}, 0);
  CompilationUnit cu(c);
  SynthesiseTket()->apply(cu);
  Circuit c1 = cu.get_circ_ref();
  REQUIRE(c1.n_gates() == 5);
}

}  // namespace test_SynthesiseTket
}  // namespace tket
