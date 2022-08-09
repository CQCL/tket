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
#include <fstream>

#include "Converters2/Converters.hpp"
#include "PauliGraph2/PauliGraph.hpp"
#include "testutil.hpp"

namespace tket {
namespace test_PauliGraph2 {

using namespace pg;

SCENARIO("Correct creation of refactored PauliGraphs") {
  GIVEN("A Clifford circuit") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::S, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::Vdg, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    PauliGraph pg = circuit_to_pauli_graph2(circ);
    REQUIRE(pg.n_vertices() == 0);
  }
  GIVEN("A 1qb circuit") {
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    circ.add_op<unsigned>(OpType::Rx, 0.6, {0});
    circ.add_op<unsigned>(OpType::Ry, 1.2, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {0});
    PauliGraph pg = circuit_to_pauli_graph2(circ);
    REQUIRE(pg.n_vertices() == 4);
    THEN("Print diagram to file") {
      std::ofstream dot_file("pauligraph.dot");
      pg.to_graphviz(dot_file);
      dot_file.close();
      // remove("pauligraph.dot");
    }
  }
}

}  // namespace test_PauliGraph2
}  // namespace tket
