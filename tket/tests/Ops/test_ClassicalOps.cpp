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
#include "Ops/ClassicalOps.hpp"
#include "Transformations/BasicOptimisation.hpp"
#include "Transformations/CliffordOptimisation.hpp"
#include "Transformations/PauliOptimisation.hpp"
#include "Transformations/PhaseOptimisation.hpp"
#include "Transformations/Transform.hpp"

namespace tket {
namespace test_ClassicalOps {

SCENARIO("Check that classical bundles work as expected") {
  GIVEN("Out bundles on a trivial circuit") {
    Circuit circ(3, 1);
    Vertex cx =
        circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {0}, 0);
    Vertex cz =
        circ.add_conditional_gate<unsigned>(OpType::CZ, {}, {1, 2}, {0}, 1);
    Vertex h =
        circ.add_conditional_gate<unsigned>(OpType::H, {}, uvec{0}, {0}, 0);
    circ.assert_valid();
    std::vector<EdgeVec> out_c_edges =
        circ.get_b_out_bundles(circ.c_inputs()[0]);
    REQUIRE(out_c_edges.size() == 1);
    REQUIRE(out_c_edges.begin()->size() == 3);
    Vertex cl = circ.c_outputs()[0];
    REQUIRE(circ.n_out_edges_of_type(cl, EdgeType::Boolean) == 0);
    REQUIRE(
        circ.get_nth_b_out_bundle(circ.c_inputs()[0], 0) ==
        out_c_edges.front());
  }
  GIVEN("In bundles on a trivial circuit") {
    Circuit circ(2, 3);
    Vertex cx =
        circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {0, 1}, 2);
    Vertex h =
        circ.add_conditional_gate<unsigned>(OpType::H, {}, {0}, {0, 1}, 3);
    circ.assert_valid();
    EdgeVec in_c_edges = circ.get_in_edges_of_type(cx, EdgeType::Boolean);
    REQUIRE(in_c_edges.size() == 2);
  }
  GIVEN("Test with a measure delaying some slicing") {
    Circuit circ(2, 1);
    circ.add_measure(0, 0);
    circ.add_conditional_gate<unsigned>(OpType::X, {}, uvec{1}, {0}, 1);
    circ.assert_valid();
    SliceVec slices = circ.get_slices();
    REQUIRE(slices.size() == 2);
  }
  GIVEN("Quantum teleportation circuit!") {
    Circuit circ(3, 2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_measure(0, 0);
    circ.add_measure(1, 1);
    circ.add_conditional_gate<unsigned>(OpType::X, {}, uvec{2}, {1}, 1);
    circ.add_conditional_gate<unsigned>(OpType::Z, {}, uvec{2}, {0}, 1);
    circ.assert_valid();
    std::vector<std::vector<OpType>> types;
    for (Circuit::SliceIterator si = circ.slice_begin(); si != circ.slice_end();
         ++si) {
      Slice sl = *si;
      std::vector<OpType> slice_types;
      for (const Vertex& v : sl) {
        slice_types.push_back(circ.get_OpType_from_Vertex(v));
      }
      types.push_back(slice_types);
    }
    std::vector<std::vector<OpType>> correct_types{
        {OpType::CX},
        {OpType::Measure, OpType::H},
        {OpType::Measure, OpType::Conditional},
        {OpType::Conditional}};
    REQUIRE(types == correct_types);
  }
  GIVEN("QASM-style entanglement swapping circuit") {
    Circuit circ(4, 2);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::CX, {2, 3});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_measure(1, 0);
    circ.add_measure(2, 1);
    circ.add_conditional_gate<unsigned>(OpType::Z, {}, {3}, {0, 1}, 1);
    circ.add_conditional_gate<unsigned>(OpType::Z, {}, {3}, {0, 1}, 3);
    circ.add_conditional_gate<unsigned>(OpType::X, {}, {3}, {0, 1}, 2);
    circ.add_conditional_gate<unsigned>(OpType::X, {}, {3}, {0, 1}, 3);
    circ.add_measure(0, 0);
    circ.add_measure(3, 1);
    circ.assert_valid();
    std::vector<std::vector<OpType>> types;
    for (Circuit::SliceIterator si = circ.slice_begin(); si != circ.slice_end();
         ++si) {
      Slice sl = *si;
      std::vector<OpType> slice_types;
      for (const Vertex& v : sl) {
        slice_types.push_back(circ.get_OpType_from_Vertex(v));
      }
      types.push_back(slice_types);
    }
    std::vector<std::vector<OpType>> correct_types{
        {OpType::H, OpType::H}, {OpType::CX, OpType::CX},
        {OpType::CX},           {OpType::Measure, OpType::H},
        {OpType::Measure},      {OpType::Conditional},
        {OpType::Conditional},  {OpType::Conditional},
        {OpType::Conditional},  {OpType::Measure, OpType::Measure}};
    REQUIRE(types == correct_types);
  }
}

SCENARIO("Test successor and predecessor methods with mixed circuits") {
  GIVEN("A purely quantum circuit") {
    Circuit circ(2);
    Vertex cx = circ.add_op<unsigned>(OpType::CX, {1, 0});
    VertexVec correct = {};
    REQUIRE(circ.get_successors_of_type(cx, EdgeType::Boolean) == correct);
    REQUIRE(circ.get_predecessors_of_type(cx, EdgeType::Boolean) == correct);
    REQUIRE(
        circ.get_successors_of_type(cx, EdgeType::Quantum) ==
        circ.get_successors(cx));
    REQUIRE(
        circ.get_predecessors_of_type(cx, EdgeType::Quantum) ==
        circ.get_predecessors(cx));
  }
  GIVEN("A purely classical circuit") {
    Circuit circ(0, 3);
    Circuit cl(0, 3);
    CircBox cbox(cl);
    Vertex cv = circ.add_box(cbox, {0, 1, 2});
    circ.assert_valid();
    VertexVec ins = circ.c_inputs();
    VertexVec outs = circ.c_outputs();
    VertexVec correct = {ins[0], ins[1], ins[2]};
    REQUIRE(circ.get_predecessors_of_type(cv, EdgeType::Classical) == correct);
    correct = {outs[0], outs[1], outs[2]};
    REQUIRE(circ.get_successors_of_type(cv, EdgeType::Classical) == correct);
    correct = {};
    REQUIRE(circ.get_predecessors_of_type(cv, EdgeType::Boolean) == correct);
    REQUIRE(circ.get_successors_of_type(cv, EdgeType::Boolean) == correct);
  }
  GIVEN("All together now") {
    Circuit circ(2, 1);
    VertexVec ins = circ.q_inputs();
    Vertex cin = circ.c_inputs()[0];
    VertexVec outs = circ.q_outputs();
    Vertex cout = circ.c_outputs()[0];
    Vertex x =
        circ.add_conditional_gate<unsigned>(OpType::X, {}, uvec{0}, {0}, 0);
    Vertex cx =
        circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {}, 0);
    Vertex m = circ.add_measure(0, 0);
    Vertex z =
        circ.add_conditional_gate<unsigned>(OpType::Z, {}, uvec{1}, {0}, 1);
    Vertex cz =
        circ.add_conditional_gate<unsigned>(OpType::CZ, {}, {1, 0}, {}, 0);
    circ.assert_valid();
    THEN("Check controlled quantum gate") {
      VertexVec correct = {cin, ins[0]};
      CHECK(circ.get_predecessors(x) == correct);
      correct = {ins[0]};
      CHECK(circ.get_predecessors_of_type(x, EdgeType::Quantum) == correct);
      correct = {cin};
      CHECK(circ.get_predecessors_of_type(x, EdgeType::Boolean) == correct);
      correct = {cx};
      CHECK(circ.get_successors(x) == correct);
      CHECK(circ.get_successors_of_type(x, EdgeType::Quantum) == correct);
      correct = {};
      CHECK(circ.get_successors_of_type(x, EdgeType::Boolean) == correct);
    }
    THEN("Check purely quantum gate") {
      VertexVec correct = {x, ins[1]};
      CHECK(circ.get_predecessors(cx) == correct);
      CHECK(circ.get_predecessors_of_type(cx, EdgeType::Quantum) == correct);
      correct = {};
      CHECK(circ.get_predecessors_of_type(cx, EdgeType::Boolean) == correct);
      correct = {m, z};
      CHECK(circ.get_successors(cx) == correct);
      CHECK(circ.get_successors_of_type(cx, EdgeType::Quantum) == correct);
      correct = {};
      CHECK(circ.get_successors_of_type(cx, EdgeType::Boolean) == correct);
    }
    THEN("Check measurement") {
      VertexVec correct = {cx, cin};
      CHECK(circ.get_predecessors(m) == correct);
      correct = {cx};
      CHECK(circ.get_predecessors_of_type(m, EdgeType::Quantum) == correct);
      correct = {};
      CHECK(circ.get_predecessors_of_type(m, EdgeType::Boolean) == correct);
      correct = {cz, cout, z};
      CHECK(circ.get_successors(m) == correct);
      correct = {cz};
      CHECK(circ.get_successors_of_type(m, EdgeType::Quantum) == correct);
      correct = {z};
      CHECK(circ.get_successors_of_type(m, EdgeType::Boolean) == correct);
    }
    THEN("Check controlled quantum gate after measurement") {
      VertexVec correct = {m, cx};
      CHECK(circ.get_predecessors(z) == correct);
      correct = {cx};
      CHECK(circ.get_predecessors_of_type(z, EdgeType::Quantum) == correct);
      correct = {m};
      CHECK(circ.get_predecessors_of_type(z, EdgeType::Boolean) == correct);
      correct = {cz};
      CHECK(circ.get_successors(z) == correct);
      CHECK(circ.get_successors_of_type(z, EdgeType::Quantum) == correct);
      correct = {};
      CHECK(circ.get_successors_of_type(z, EdgeType::Boolean) == correct);
    }
    THEN("Check purely quantum gate after measurement") {
      VertexVec correct = {z, m};
      CHECK(circ.get_predecessors(cz) == correct);
      CHECK(circ.get_predecessors_of_type(cz, EdgeType::Quantum) == correct);
      correct = {};
      CHECK(circ.get_predecessors_of_type(cz, EdgeType::Boolean) == correct);
      correct = {outs[1], outs[0]};
      CHECK(circ.get_successors(cz) == correct);
      CHECK(circ.get_successors_of_type(cz, EdgeType::Quantum) == correct);
      correct = {};
      CHECK(circ.get_successors_of_type(cz, EdgeType::Boolean) == correct);
    }
    THEN("Check boundaries") {
      VertexVec correct = {x};
      CHECK(circ.get_successors(ins[0]) == correct);
      correct = {cx};
      CHECK(circ.get_successors(ins[1]) == correct);
      correct = {m, x};
      CHECK(circ.get_successors(cin) == correct);
      correct = {cz};
      CHECK(circ.get_predecessors(outs[0]) == correct);
      CHECK(circ.get_predecessors(outs[1]) == correct);
      correct = {m};
      CHECK(circ.get_predecessors(cout) == correct);
    }
  }
}

SCENARIO("Try slicing on a Circuit with classical data") {
  Circuit circ(2, 1);
  Circuit cl(0, 1);
  CircBox cbox(cl);
  circ.add_box(cbox, uvec{0});
  circ.add_conditional_gate<unsigned>(OpType::H, {}, uvec{0}, {0}, 1);
  circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {0}, 1);
  circ.assert_valid();

  SliceVec sv = circ.get_slices();
  REQUIRE(sv.size() == 3);
  for (const Slice& s : sv) REQUIRE(s.size() == 1);
  check_command_types(
      circ, {OpType::CircBox, OpType::Conditional, OpType::Conditional});

  GIVEN(
      "A circuit with ClOutput does not include the ClOutput in the "
      "commands") {
    Circuit circ(1, 1);
    circ.add_measure(0, 0);
    circ.assert_valid();
    check_command_types(circ, {OpType::Measure});
  }
}

SCENARIO("Remove redundancies on a circuit with classical controls") {
  Circuit circ(2, 2);
  WHEN(
      "Conditional gates are considered boxes and so are not hit by "
      "remove_redundancies") {
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {0}, 0);
    circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {0}, 1);
    REQUIRE(!Transforms::remove_redundancies().apply(circ));
    circ.assert_valid();
    REQUIRE(circ.n_gates() == 2);
  }
}

SCENARIO(
    "Test pauli gadget opt. throws properly on circuit with classical "
    "controls") {
  GIVEN("Circuit to run pairwise Transforms on") {
    Circuit circ(4, 1);
    Vertex cx1 = circ.add_op<unsigned>(OpType::CX, {0, 3});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::CZ, {0, 3});
    circ.add_op<unsigned>(OpType::CX, {3, 2});
    WHEN(
        "Pauli gadget optimisation does not throw on simple unitary "
        "circuit") {
      REQUIRE_NOTHROW(Transforms::pairwise_pauli_gadgets().apply(circ));
    }
    WHEN(
        "Pauli gadget optimisation does throw on classically controlled "
        "circuit") {
      circ.add_conditional_gate<unsigned>(OpType::CX, {}, {3, 2}, {0}, 0);
      REQUIRE_THROWS_AS(
          Transforms::pairwise_pauli_gadgets().apply(circ), CircuitInvalidity);
    }
    WHEN(
        "Pauli gadget optimisation does not throw on circuit with "
        "measures at the end") {
      circ.add_measure(0, 0);
      REQUIRE_NOTHROW(Transforms::pairwise_pauli_gadgets().apply(circ));
      circ.assert_valid();
    }
  }

  GIVEN("Circuit with classical controls on CXs next to phase-gadgets") {
    Circuit circ(5, 1);
    WHEN("Normal phase gadget optimisation") {
      circ.add_op<unsigned>(OpType::CX, {0, 1});
      circ.add_op<unsigned>(OpType::PhaseGadget, 0.3, {1, 2, 3, 4});
      circ.add_op<unsigned>(OpType::CX, {0, 1});
      REQUIRE(Transforms::smash_CX_PhaseGadgets().apply(circ));
      REQUIRE(circ.n_gates() == 1);
      REQUIRE(circ.count_gates(OpType::PhaseGadget) == 1);
    }
    WHEN("Add classical wire to first cx") {
      circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {0}, 1);
      circ.add_op<unsigned>(OpType::PhaseGadget, 0.3, {1, 2, 3, 4});
      circ.add_op<unsigned>(OpType::CX, {0, 1});
      REQUIRE(!Transforms::smash_CX_PhaseGadgets().apply(circ));
    }
    WHEN("Add classical wire to second cx") {
      circ.add_op<unsigned>(OpType::CX, {0, 1});
      circ.add_op<unsigned>(OpType::PhaseGadget, 0.3, {1, 2, 3, 4});
      circ.add_conditional_gate<unsigned>(OpType::CX, {}, {0, 1}, {0}, 0);
      REQUIRE(!Transforms::smash_CX_PhaseGadgets().apply(circ));
      circ.assert_valid();
    }
  }
}

SCENARIO("PI-copy rule") {
  GIVEN("A circuit with PI-copy rules to be done") {
    Circuit circ(2, 1);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    WHEN("Do pi copy rule") {
      circ.add_op<unsigned>(OpType::X, {0});
      REQUIRE(Transforms::singleq_clifford_sweep().apply(circ));
      REQUIRE(circ.n_gates() == 3);
    }
    WHEN("Add classical wires") {
      circ.add_conditional_gate<unsigned>(OpType::X, {}, uvec{0}, {0}, 1);
      REQUIRE(!Transforms::singleq_clifford_sweep().apply(circ));
      circ.assert_valid();
    }
  }
}

SCENARIO("Classical wires for appending circuits") {
  GIVEN("Append easy circuit") {
    Circuit circ(2, 1);
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_measure(0, 0);
    Circuit circ2(1, 1);
    circ2.add_conditional_gate<unsigned>(OpType::X, {}, uvec{0}, {0}, 1);
    REQUIRE_NOTHROW(circ.append_qubits(circ2, {0}, {0}));
    circ.assert_valid();
    REQUIRE(circ.n_gates() == 3);
    REQUIRE(circ.get_slices().size() == 3);
    REQUIRE(circ.n_bits() == 1);
  }
  GIVEN("Another circuit") {
    Circuit circ(6, 5);
    for (unsigned i = 0; i < 5; ++i) {
      circ.add_measure(Qubit(q_default_reg(), i), Bit(c_default_reg(), i));
    }
    Circuit circ2(1, 2);
    circ2.add_op<unsigned>(OpType::X, {0});
    circ2.add_op<unsigned>(OpType::Z, {0});

    REQUIRE_NOTHROW(circ.append_qubits(circ2, {5}, {1, 0}));
    circ.assert_valid();
    REQUIRE(circ.n_gates() == 7);
    REQUIRE(circ.depth() == 2);
  }
  GIVEN("Invalid circuit appending") {
    Circuit circ(1, 0);
    circ.add_op<unsigned>(OpType::Y, {0});

    Circuit circ2(1, 2);
    Vertex z =
        circ2.add_conditional_gate<unsigned>(OpType::Z, {}, uvec{0}, {0}, 1);
    WHEN("Append with an invalid cmap") {
      REQUIRE_THROWS_AS(
          circ.append_qubits(circ2, {0}, {0, 0}), CircuitInvalidity);
    }
    circ2.assert_valid();
  }
}

SCENARIO("Reverse slicing of mixed circuits") {
  GIVEN("A circuit with a measure and no classically-controlled gates") {
    Circuit circ(4, 1);
    Vertex x = circ.add_op<unsigned>(OpType::X, {0});
    Vertex y = circ.add_op<unsigned>(OpType::Y, {1});
    Vertex m = circ.add_measure(1, 0);
    Vertex z = circ.add_op<unsigned>(OpType::Z, {3});
    SliceVec backwards = circ.get_reverse_slices();
    REQUIRE(backwards.size() == 2);
    Slice correct = {m, x, z};
    REQUIRE(backwards[0] == correct);
    correct = {y};
    REQUIRE(backwards[1] == correct);
  }
  GIVEN("A circuit with classically-controlled gates") {
    Circuit circ(2, 1);
    Vertex x = circ.add_op<unsigned>(OpType::X, {0});
    Vertex y =
        circ.add_conditional_gate<unsigned>(OpType::Y, {}, uvec{1}, {0}, 1);
    circ.assert_valid();
    REQUIRE_NOTHROW(circ.get_reverse_slices());
  }
}

SCENARIO("Pure classical operations") {
  GIVEN("A pure classical circuit") {
    std::vector<uint32_t> and_table = {0, 1, 2, 7, 0, 1, 2, 7};
    std::shared_ptr<ClassicalTransformOp> and_ttop =
        std::make_shared<ClassicalTransformOp>(3, and_table);
    for (unsigned i = 0; i < 2; i++) {
      for (unsigned j = 0; j < 2; j++) {
        for (unsigned k = 0; k < 2; k++) {
          std::vector<bool> y = and_ttop->eval({(bool)i, (bool)j, (bool)k});
          REQUIRE(y[0] == i);
          REQUIRE(y[1] == j);
          REQUIRE(y[2] == (i & j));
        }
      }
    }

    uint32_t a = 2, b = 6;
    std::shared_ptr<RangePredicateOp> rpop =
        std::make_shared<RangePredicateOp>(3, a, b);
    for (uint32_t x = 0; x < 8; x++) {
      REQUIRE(
          rpop->eval(
              {(bool)(x & 1), (bool)((x >> 1) & 1), (bool)((x >> 2) & 1)})[0] ==
          (x >= a && x <= b));
    }

    Circuit circ(1, 4);
    circ.add_op<unsigned>(OpType::H, {0});
    Vertex v_and_ttop_0 = circ.add_op<unsigned>(and_ttop, {0, 1, 2});
    Vertex v_and_ttop_1 = circ.add_op<unsigned>(and_ttop, {1, 2, 3});
    Vertex v_rpop = circ.add_op<unsigned>(rpop, {0, 1, 2, 3});
    Vertex v_andop = circ.add_op<unsigned>(AndOp(), {2, 3, 0});
    Vertex v_orop = circ.add_op<unsigned>(OrOp(), {0, 1, 2});
    Vertex v_notop = circ.add_op<unsigned>(NotOp(), {2, 3});
    Vertex v_clx = circ.add_op<unsigned>(ClassicalX(), {1});
    Vertex v_clcx = circ.add_op<unsigned>(ClassicalCX(), {0, 1});
    Vertex v_andwop = circ.add_op<unsigned>(AndWithOp(), {2, 3});
    Vertex v_orwop = circ.add_op<unsigned>(OrWithOp(), {1, 0});
    circ.assert_valid();
    REQUIRE(circ.get_commands().size() == 11);
    std::vector<EdgeVec> and_ttop_1_b_out =
        circ.get_b_out_bundles(v_and_ttop_1);
    REQUIRE(and_ttop_1_b_out.size() == 3);
    REQUIRE(and_ttop_1_b_out[0].size() == 2);
    REQUIRE(and_ttop_1_b_out[1].size() == 2);
    REQUIRE(and_ttop_1_b_out[2].size() == 0);
    REQUIRE(circ.n_out_edges_of_type(v_and_ttop_0, EdgeType::Boolean) == 1);
    REQUIRE(circ.n_out_edges_of_type(v_and_ttop_0, EdgeType::Classical) == 3);
    REQUIRE(circ.n_in_edges_of_type(v_orwop, EdgeType::Boolean) == 1);
    REQUIRE(circ.n_in_edges_of_type(v_orwop, EdgeType::Classical) == 1);
    REQUIRE(circ.n_in_edges_of_type(v_orop, EdgeType::Boolean) == 2);
    REQUIRE(circ.n_in_edges_of_type(v_orop, EdgeType::Classical) == 1);
    REQUIRE(circ.n_out_edges_of_type(v_orop, EdgeType::Boolean) == 2);
    REQUIRE(circ.n_out_edges_of_type(v_orop, EdgeType::Classical) == 1);
  }
  GIVEN("A multi-bit operation") {
    Circuit circ(0, 6);
    std::shared_ptr<MultiBitOp> mbop = std::make_shared<MultiBitOp>(AndOp(), 2);
    circ.add_op<unsigned>(mbop, {0, 1, 2, 3, 4, 5});
    circ.assert_valid();
    REQUIRE(circ.count_gates(OpType::MultiBit) == 1);
    std::vector<bool> x = {0, 1, 1, 1};
    std::vector<bool> y = mbop->eval(x);
    REQUIRE(y.size() == 2);
    REQUIRE(!y[0]);
    REQUIRE(y[1]);
  }
}

}  // namespace test_ClassicalOps
}  // namespace tket
