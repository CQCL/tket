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
#include <memory>

#include "Circuit/Boxes.hpp"
#include "Gate/SymTable.hpp"
#include "Placement/Placement.hpp"
#include "Predicates/CompilationUnit.hpp"
#include "Predicates/Predicates.hpp"
#include "testutil.hpp"

namespace tket {
namespace test_Predicates {

SCENARIO("Test out basic Predicate useage") {
  GIVEN("GateSetPredicate") {
    OpTypeSet ots = {OpType::CX};
    PredicatePtr gsp = std::make_shared<GateSetPredicate>(ots);

    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    REQUIRE(gsp->verify(circ));
    circ.add_op<unsigned>(OpType::Collapse, {0});
    REQUIRE(!gsp->verify(circ));

    OpTypeSet ots2 = {OpType::CX, OpType::Z};
    PredicatePtr gsp2 = std::make_shared<GateSetPredicate>(ots2);
    REQUIRE(gsp->implies(*gsp2));

    OpTypeSet ots3 = {OpType::CX, OpType::Ry};
    PredicatePtr gsp3 = std::make_shared<GateSetPredicate>(ots3);
    REQUIRE(!gsp2->implies(*gsp3));
  }
  GIVEN("NoClassicalControlPredicate") {
    PredicatePtr pp = std::make_shared<NoClassicalControlPredicate>();
    Circuit circ(1, 1);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_measure(0, 0);
    REQUIRE(pp->verify(circ));
    WHEN("Add a conditional gate") {
      circ.add_conditional_gate<unsigned>(OpType::X, {}, uvec{0}, {0}, 1);
      REQUIRE_FALSE(pp->verify(circ));
    }
    WHEN("Add a CircBox without any conditionals") {
      CircBox cbox(circ);
      Circuit larger(2, 2);
      larger.add_op<unsigned>(OpType::CX, {0, 1});
      larger.add_box(cbox, {0, 0});
      REQUIRE(pp->verify(larger));
    }
    WHEN("Add a CircBox with conditionals") {
      circ.add_conditional_gate<unsigned>(OpType::X, {}, uvec{0}, {0}, 1);
      CircBox cbox(circ);
      Circuit larger(2, 2);
      larger.add_op<unsigned>(OpType::CX, {0, 1});
      larger.add_box(cbox, {0, 0});
      REQUIRE_FALSE(pp->verify(larger));
    }
    PredicatePtr pp2 = std::make_shared<NoClassicalControlPredicate>();
    REQUIRE(pp->implies(*pp2));
  }
  GIVEN("NoClassicalBitsPredicate") {
    PredicatePtr pp = std::make_shared<NoClassicalBitsPredicate>();
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::X, {0});
    REQUIRE(pp->verify(circ));
    Vertex in = circ.add_vertex(OpType::ClInput);
    Vertex out = circ.add_vertex(OpType::ClOutput);
    circ.add_edge({in, 0}, {out, 0}, EdgeType::Classical);
    circ.boundary.insert({Bit(0), in, out});
    REQUIRE(!pp->verify(circ));
  }
  GIVEN("MaxTwoQubitGatesPredicate") {
    PredicatePtr pp = std::make_shared<MaxTwoQubitGatesPredicate>();
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CZ, {0, 1});
    REQUIRE(pp->verify(circ));
    circ.add_op<unsigned>(OpType::CCX, {0, 1, 2});
    REQUIRE(!pp->verify(circ));
  }
  GIVEN("NoFastFeedforwardPredicate") {
    PredicatePtr pp = std::make_shared<NoFastFeedforwardPredicate>();
    Circuit circ(2, 2);
    circ.add_conditional_gate<unsigned>(OpType::H, {}, uvec{0}, {0}, 0);
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {0}, 0);
    WHEN("Normal gate sequencing") {
      circ.add_measure(1, 0);
      REQUIRE(pp->verify(circ));
      circ.add_conditional_gate<unsigned>(OpType::X, {}, uvec{0}, {0}, 0);
      REQUIRE_FALSE(pp->verify(circ));
    }
    WHEN("Adding a CircBox that measures") {
      Circuit inner(1, 1);
      inner.add_measure(0, 0);
      CircBox cbox(inner);
      circ.add_box(cbox, {1, 1});
      circ.add_conditional_gate<unsigned>(OpType::X, {}, uvec{1}, {0}, 0);
      REQUIRE(pp->verify(circ));
      circ.add_conditional_gate<unsigned>(OpType::Y, {}, uvec{0}, {1}, 0);
      REQUIRE_FALSE(pp->verify(circ));
    }
    WHEN("Adding a CircBox that needs feed-forward") {
      Circuit inner(1, 1);
      inner.add_conditional_gate<unsigned>(OpType::X, {}, uvec{0}, {0}, 0);
      CircBox cbox(inner);
      circ.add_measure(1, 0);
      circ.add_box(cbox, {1, 1});
      REQUIRE(pp->verify(circ));
      circ.add_box(cbox, {0, 0});
      REQUIRE_FALSE(pp->verify(circ));
    }
  }
  GIVEN("DefaultRegisterPredicate") {
    PredicatePtr pp = std::make_shared<DefaultRegisterPredicate>();
    Circuit circ;
    REQUIRE(pp->verify(circ));
    circ.add_q_register(q_default_reg(), 3);
    circ.add_c_register(c_default_reg(), 2);
    REQUIRE(pp->verify(circ));
    Qubit unusual("unusual", 4);
    circ.add_qubit(unusual);
    REQUIRE(!pp->verify(circ));
    circ.rename_units<Qubit, Qubit>({{unusual, Qubit(7)}});
    REQUIRE(pp->verify(circ));
  }
  GIVEN("GlobalPhasedXPredicate") {
    PredicatePtr pp = std::make_shared<GlobalPhasedXPredicate>();
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    REQUIRE(pp->verify(circ));
    WHEN("Add a non-global NPhasedX gate") {
      circ.add_op<unsigned>(OpType::NPhasedX, {0.2, 0.3}, {0, 1});
      REQUIRE_FALSE(pp->verify(circ));
    }
    WHEN("Add two global NPhasedX gates") {
      circ.add_op<unsigned>(OpType::NPhasedX, {0.2, 0.3}, {0, 1, 2});
      circ.add_op<unsigned>(OpType::NPhasedX, {0.5, 0.2}, {1, 0, 2});
      REQUIRE(pp->verify(circ));
    }
    WHEN("Add some global and some non-global NPhasedX gates") {
      circ.add_op<unsigned>(OpType::NPhasedX, {0.2, 0.3}, {0, 1, 2});
      circ.add_op<unsigned>(OpType::NPhasedX, {0.5, 0.2}, {1, 0, 2});
      REQUIRE(pp->verify(circ));
      circ.add_op<unsigned>(OpType::NPhasedX, {0.5, 0.2}, {1, 0});
      circ.add_op<unsigned>(OpType::NPhasedX, {0.2, 0.3}, {0, 1, 2});
      REQUIRE_FALSE(pp->verify(circ));
    }
  }
}

SCENARIO("Make sure combining predicates for `implies` throws as expected") {
  PredicatePtr pp1 = std::make_shared<MaxTwoQubitGatesPredicate>();
  PredicatePtr pp2 = std::make_shared<NoClassicalBitsPredicate>();
  REQUIRE_THROWS_AS(pp1->implies(*pp2), IncorrectPredicate);
}

SCENARIO("Test CliffordCircuitPredicate") {
  Circuit circ(8);
  circ.add_op<unsigned>(OpType::S, {1});
  circ.add_op<unsigned>(OpType::Rx, 1.5, {2});
  circ.add_op<unsigned>(OpType::CX, {1, 7});
  circ.add_op<unsigned>(OpType::CX, {2, 4});
  circ.add_op<unsigned>(OpType::Rz, 0.5, {1});
  circ.add_op<unsigned>(OpType::Rx, 0.5, {2});
  circ.add_op<unsigned>(OpType::CX, {1, 3});
  circ.add_op<unsigned>(OpType::CX, {5, 6});
  circ.add_op<unsigned>(OpType::CX, {6, 7});
  circ.add_op<unsigned>(OpType::H, {2});
  circ.add_barrier({3, 4, 5});
  circ.add_op<unsigned>(OpType::Rx, -0.5, {0});
  circ.add_op<unsigned>(OpType::Ry, 1.5, {1});
  circ.add_op<unsigned>(OpType::Rz, 0.5, {2});
  circ.add_op<unsigned>(OpType::U1, 1.0, {3});
  circ.add_op<unsigned>(OpType::U2, {-0.5, 1.5}, {4});
  circ.add_op<unsigned>(OpType::U3, {0., 1.5, 4.5}, {5});
  circ.add_op<unsigned>(OpType::TK1, {-0.5, 1.5, 4.}, {6});
  circ.add_op<unsigned>(OpType::TK2, {1.5, 2.5, -1.}, {7, 0});
  circ.add_op<unsigned>(OpType::XXPhase, -0.5, {1, 2});
  circ.add_op<unsigned>(OpType::YYPhase, 0.5, {2, 3});
  circ.add_op<unsigned>(OpType::ZZPhase, 0., {3, 4});
  circ.add_op<unsigned>(OpType::XXPhase3, 1.0, {4, 5, 6});
  circ.add_op<unsigned>(OpType::PhasedX, {-0.5, 0.5}, {5});
  circ.add_op<unsigned>(OpType::NPhasedX, {1.5, 1.5}, {6, 7});
  circ.add_op<unsigned>(OpType::ISWAP, 1.0, {0, 1});
  circ.add_op<unsigned>(OpType::ESWAP, 2.0, {2, 3});
  circ.add_op<unsigned>(OpType::PhasedISWAP, {1.5, 0.}, {4, 5});
  circ.add_op<unsigned>(OpType::FSim, {0.5, 1.}, {6, 7});
  CircBox cbox(circ);
  Circuit circ1(8);
  circ1.add_box(cbox, {0, 1, 2, 3, 4, 5, 6, 7});
  PauliExpBox pebox({Pauli::Y, Pauli::Z}, 0.5);
  circ1.add_box(pebox, {0, 1});
  Circuit setup(2);
  Sym a = SymTable::fresh_symbol("a");
  setup.add_op<unsigned>(OpType::Rx, {a}, {0});
  setup.add_op<unsigned>(OpType::CX, {0, 1});
  setup.add_op<unsigned>(OpType::Ry, 0.5, {0});
  composite_def_ptr_t def = CompositeGateDef::define_gate("g", setup, {a});
  CustomGate cgbox(def, {1.5});
  circ1.add_box(cgbox, {2, 3});
  Eigen::Matrix2cd U;
  U << 0.5 - 0.5 * i_, 0.5 + 0.5 * i_, 0.5 + 0.5 * i_, 0.5 - 0.5 * i_;
  Unitary1qBox u1box(U);
  circ1.add_box(u1box, {4});
  PredicatePtr ccp = std::make_shared<CliffordCircuitPredicate>();
  REQUIRE(ccp->verify(circ1));
  Circuit circ2(2);
  circ2.add_op<unsigned>(OpType::TK2, {1.5, 2.5, -1.01}, {0, 1});
  REQUIRE(!ccp->verify(circ2));
}

SCENARIO("Test routing-related predicates' meet and implication") {
  Node n0("test", 0);
  Node n1("test", 1);
  Node n2("test", 2);
  Node n3("test", 3);
  Architecture arc1({{n0, n1}, {n1, n2}});
  Architecture arc2({{n0, n1}, {n1, n2}, {n0, n2}});
  Architecture arc3({{n0, n2}, {n0, n1}});
  Architecture arc4({{n2, n0}, {n0, n1}});
  Architecture arc5({n0, n1, n2, n3});
  arc5.add_connection(n0, n1);
  arc5.add_connection(n1, n2);

  Circuit circ(3);
  circ.add_op<unsigned>(OpType::CX, {0, 1});
  circ.add_op<unsigned>(OpType::CX, {1, 2});
  circ.add_op<unsigned>(OpType::BRIDGE, {2, 1, 0});
  reassign_boundary(circ, node_vector_t{n0, n1, n2});

  Circuit circ2(3);
  add_2qb_gates(circ2, OpType::CX, {{0, 1}, {0, 2}, {1, 2}});
  reassign_boundary(circ2, node_vector_t{n0, n1, n2});
  GIVEN("Connectivity Predicates") {
    PredicatePtr con1 = std::make_shared<ConnectivityPredicate>(arc1);
    PredicatePtr con2 = std::make_shared<ConnectivityPredicate>(arc2);
    PredicatePtr con3 = std::make_shared<ConnectivityPredicate>(arc3);
    PredicatePtr con4 = std::make_shared<ConnectivityPredicate>(arc4);
    PredicatePtr con5 = std::make_shared<ConnectivityPredicate>(arc5);
    WHEN("Test implies") {
      REQUIRE(con1->implies(*con2));
      REQUIRE(con4->implies(*con3));  // directedness doesn't matter
      REQUIRE_FALSE(con1->implies(*con3));
    }
    WHEN("Test implies (isolated nodes)") {
      // https://github.com/CQCL/tket/issues/88
      REQUIRE(con1->implies(*con5));
      REQUIRE_FALSE(con5->implies(*con1));
    }
    WHEN("Test meet") {
      PredicatePtr meet_a = con1->meet(*con2);
      REQUIRE(meet_a->verify(circ));
      REQUIRE_FALSE(meet_a->verify(circ2));
    }
  }
  GIVEN("Directedness Predicates") {
    PredicatePtr con1 = std::make_shared<DirectednessPredicate>(arc1);
    PredicatePtr con2 = std::make_shared<DirectednessPredicate>(arc2);
    PredicatePtr con3 = std::make_shared<DirectednessPredicate>(arc3);
    PredicatePtr con4 = std::make_shared<DirectednessPredicate>(arc4);
    WHEN("Test verify") { REQUIRE_FALSE(con1->verify(circ)); }
    WHEN("Test implies") {
      REQUIRE(con1->implies(*con2));
      REQUIRE_FALSE(con4->implies(*con3));  // directedness *does* matter
      REQUIRE_FALSE(con1->implies(*con3));
      REQUIRE_FALSE(con4->implies(*con1));
    }
    WHEN("Test meet") {
      PredicatePtr meet_a = con1->meet(*con2);
      REQUIRE_FALSE(meet_a->verify(circ));
      REQUIRE_FALSE(meet_a->verify(circ2));
    }
  }
}

SCENARIO("Test basic functionality of CompilationUnit") {
  GIVEN("A satisfied/unsatisfied predicate in a CompilationUnit") {
    OpTypeSet ots = {OpType::CX};
    PredicatePtr gsp = std::make_shared<GateSetPredicate>(ots);
    PredicatePtrMap ppm{CompilationUnit::make_type_pair(gsp)};

    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});

    CompilationUnit cu(circ, ppm);
    REQUIRE(cu.check_all_predicates());
    PredicateCache pc = cu.get_cache_ref();
    REQUIRE(pc.size() == 1);
    REQUIRE(pc.begin()->second.second);  // cache has predicate satisfied

    OpTypeSet ots2 = {OpType::CZ};
    PredicatePtr gsp2 = std::make_shared<GateSetPredicate>(ots2);
    PredicatePtrMap ppm2{CompilationUnit::make_type_pair(gsp2)};

    CompilationUnit cu2(circ, ppm2);
    REQUIRE(!cu2.check_all_predicates());
    PredicateCache pc2 = cu2.get_cache_ref();
    REQUIRE(pc2.size() == 1);
    REQUIRE(!pc2.begin()->second.second);  // cache has predicate unsatisfied
  }
}

SCENARIO("Test PlacementPredicate") {
  GIVEN(
      "Does the GraphPlacement class correctly modify Circuits and return "
      "maps?") {
    Architecture test_arc({{0, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 5}});

    PredicatePtr placement_pred =
        std::make_shared<PlacementPredicate>(test_arc);

    Circuit test_circ(6);
    add_2qb_gates(
        test_circ, OpType::CX,
        {{0, 1}, {2, 1}, {3, 1}, {2, 5}, {3, 4}, {0, 5}});
    WHEN("Base Placement") {
      Placement base_p(test_arc);
      REQUIRE(!placement_pred->verify(test_circ));
      base_p.place(test_circ);
      REQUIRE(placement_pred->verify(test_circ));
    }
    WHEN("Line Placement") {
      LinePlacement line_p(test_arc);
      REQUIRE(!placement_pred->verify(test_circ));
      line_p.place(test_circ);
      REQUIRE(placement_pred->verify(test_circ));
    }
    WHEN("Graph Placement") {
      GraphPlacement graph_p(test_arc);
      REQUIRE(!placement_pred->verify(test_circ));
      graph_p.place(test_circ);
      REQUIRE(placement_pred->verify(test_circ));
    }
    WHEN("Noise Placement") {
      NoiseAwarePlacement noise_p(test_arc);
      REQUIRE(!placement_pred->verify(test_circ));
      noise_p.place(test_circ);
      REQUIRE(placement_pred->verify(test_circ));
    }
  }
}

SCENARIO("Verifying whether or not circuits have mid-circuit measurements") {
  PredicatePtr mid_meas_pred = std::make_shared<NoMidMeasurePredicate>();
  GIVEN("No measurements") {
    Circuit c(2, 2);
    c.add_op<unsigned>(OpType::Z, {0});
    c.add_op<unsigned>(OpType::CX, {1, 0});
    c.add_op<unsigned>(OpType::Rx, 0.3, {1});
    REQUIRE(mid_meas_pred->verify(c));
  }
  GIVEN("Some mid-measurements") {
    Circuit c(2, 2);
    c.add_op<unsigned>(OpType::Z, {0});
    c.add_op<unsigned>(OpType::Measure, {0, 0});
    c.add_op<unsigned>(OpType::CX, {1, 0});
    c.add_op<unsigned>(OpType::Rx, 0.3, {1});
    c.add_op<unsigned>(OpType::Measure, {1, 1});
    REQUIRE_FALSE(mid_meas_pred->verify(c));
  }
  GIVEN("All measurements at end") {
    Circuit c(2, 2);
    c.add_op<unsigned>(OpType::Z, {0});
    c.add_op<unsigned>(OpType::CX, {1, 0});
    c.add_op<unsigned>(OpType::Rx, 0.3, {1});
    c.add_op<unsigned>(OpType::Measure, {0, 0});
    c.add_op<unsigned>(OpType::Measure, {1, 1});
    REQUIRE(mid_meas_pred->verify(c));
  }
  GIVEN("Measures by output with feedforward") {
    Circuit c(2, 2);
    c.add_op<unsigned>(OpType::CZ, {0, 1});
    c.add_op<unsigned>(OpType::Measure, {0, 0});
    c.add_conditional_gate<unsigned>(OpType::Z, {}, {1}, {0}, 1);
    REQUIRE_FALSE(mid_meas_pred->verify(c));
  }
  GIVEN("Measures at end on same bit") {
    Circuit c(2, 2);
    c.add_op<unsigned>(OpType::Z, {0});
    c.add_op<unsigned>(OpType::CX, {1, 0});
    c.add_op<unsigned>(OpType::Rx, 0.3, {1});
    c.add_op<unsigned>(OpType::Measure, {0, 0});
    c.add_op<unsigned>(OpType::Measure, {1, 0});
    REQUIRE_FALSE(mid_meas_pred->verify(c));
  }
  GIVEN("Measure in CircBox") {
    Circuit inner(1, 1);
    inner.add_op<unsigned>(OpType::X, {0});
    inner.add_op<unsigned>(OpType::Measure, {0, 0});
    CircBox cbox(inner);
    Circuit c(2, 2);
    c.add_box(cbox, {0, 0});
    c.add_op<unsigned>(OpType::Measure, {1, 1});
    REQUIRE(mid_meas_pred->verify(c));
    c.add_op<unsigned>(OpType::Z, {0});
    REQUIRE_FALSE(mid_meas_pred->verify(c));
  }
  GIVEN("Subsequent gates in CircBox") {
    Circuit inner(1);
    inner.add_op<unsigned>(OpType::X, {0});
    CircBox cbox(inner);
    Circuit c(2, 2);
    c.add_op<unsigned>(OpType::Z, {0});
    c.add_op<unsigned>(OpType::CX, {1, 0});
    c.add_op<unsigned>(OpType::Rx, 0.3, {1});
    c.add_op<unsigned>(OpType::Measure, {0, 0});
    c.add_box(cbox, {0});
    REQUIRE_FALSE(mid_meas_pred->verify(c));
  }
  GIVEN("Identity CircBox after measures") {
    Circuit inner(2, 1);
    inner.add_op<unsigned>(OpType::X, {0});
    CircBox cbox(inner);
    Circuit c(2, 2);
    c.add_op<unsigned>(OpType::Z, {0});
    c.add_op<unsigned>(OpType::CX, {1, 0});
    c.add_op<unsigned>(OpType::Rx, 0.3, {1});
    c.add_op<unsigned>(OpType::Measure, {0, 0});
    c.add_box(cbox, {1, 0, 0});
    REQUIRE(mid_meas_pred->verify(c));
  }
}

}  // namespace test_Predicates
}  // namespace tket
