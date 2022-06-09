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

#include "Mapping/LexiRoute.hpp"
#include "Mapping/MappingManager.hpp"
#include "Mapping/Verification.hpp"
#include "Placement/Placement.hpp"
#include "testutil.hpp"

namespace tket {
SCENARIO(
    "Test validity of circuit against architecture using "
    "respects_connectivity_constraints method.",
    "[routing]") {
  Architecture arc({{1, 0}, {1, 2}});

  GIVEN("A simple CX circuit and a line_placement map.") {
    Circuit circ(5);
    add_2qb_gates(circ, OpType::CX, {{0, 1}, {0, 3}, {2, 4}, {1, 4}, {0, 4}});
    Architecture test_arc({{0, 1}, {1, 2}, {2, 3}, {3, 4}});
    LinePlacement lp_obj(test_arc);
    lp_obj.place(circ);
    MappingManager mm(std::make_shared<Architecture>(test_arc));

    REQUIRE(
        mm.route_circuit(circ, {std::make_shared<LexiRouteRoutingMethod>()}));
    CHECK(respects_connectivity_constraints(circ, test_arc, false));
  }
  GIVEN("A failing case, undirected") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CX, {0, 2});
    reassign_boundary(circ);
    REQUIRE_FALSE(respects_connectivity_constraints(circ, arc, false));
  }
  GIVEN("A working case, undirected") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    reassign_boundary(circ);
    REQUIRE(respects_connectivity_constraints(circ, arc, false));
  }
  GIVEN("A failing case, directed") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    reassign_boundary(circ);
    REQUIRE_FALSE(respects_connectivity_constraints(circ, arc, true));
  }
  GIVEN("A working case, directed") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CX, {1, 0});
    reassign_boundary(circ);
    REQUIRE(respects_connectivity_constraints(circ, arc, true));
  }
  GIVEN("A failing case, undirected, with SWAP") {
    Circuit circ(3);
    Vertex swap_v = circ.add_op<unsigned>(OpType::SWAP, {1, 2});

    EdgeVec swap_outs = circ.get_all_out_edges(swap_v);
    circ.dag[swap_outs[0]].ports.first = 1;
    circ.dag[swap_outs[1]].ports.first = 0;

    circ.add_op<unsigned>(OpType::CX, {0, 1});
    reassign_boundary(circ);
    REQUIRE_FALSE(respects_connectivity_constraints(circ, arc, false));
  }
  GIVEN("A working case, undirected, with SWAP") {
    Circuit circ(3);
    Vertex swap_v = circ.add_op<unsigned>(OpType::SWAP, {1, 2});

    EdgeVec swap_outs = circ.get_all_out_edges(swap_v);
    circ.dag[swap_outs[0]].ports.first = 1;
    circ.dag[swap_outs[1]].ports.first = 0;

    circ.add_op<unsigned>(OpType::CX, {0, 2});
    reassign_boundary(circ);
    REQUIRE(respects_connectivity_constraints(circ, arc, false));
  }
  GIVEN("A failing case, directed, with SWAP") {
    Circuit circ(3);
    Vertex swap_v = circ.add_op<unsigned>(OpType::SWAP, {1, 0});

    EdgeVec swap_outs = circ.get_all_out_edges(swap_v);
    circ.dag[swap_outs[0]].ports.first = 1;
    circ.dag[swap_outs[1]].ports.first = 0;

    circ.add_op<unsigned>(OpType::CX, {1, 0});
    reassign_boundary(circ);
    REQUIRE_FALSE(respects_connectivity_constraints(circ, arc, true));
  }
  GIVEN("A working case, directed, with SWAP") {
    Circuit circ(3);
    Vertex swap_v = circ.add_op<unsigned>(OpType::SWAP, {1, 0});

    EdgeVec swap_outs = circ.get_all_out_edges(swap_v);
    circ.dag[swap_outs[0]].ports.first = 1;
    circ.dag[swap_outs[1]].ports.first = 0;

    circ.add_op<unsigned>(OpType::CX, {0, 1});
    reassign_boundary(circ);
    REQUIRE(respects_connectivity_constraints(circ, arc, false));
  }
}
}  // namespace tket
