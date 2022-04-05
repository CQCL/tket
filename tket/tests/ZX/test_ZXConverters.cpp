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

#include <catch2/catch.hpp>

#include "../testutil.hpp"
#include "Converters/Converters.hpp"
#include "Gate/SymTable.hpp"

namespace tket {
namespace zx {
namespace test_ZXConverters {
SCENARIO("Check converting gates to spiders") {
  GIVEN("X") {
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::X, {0});
    ZXDiagram zx = circuit_to_zx(circ);
    ZXVertVec boundary = zx.get_boundary(ZXType::Input, QuantumType::Quantum);
    ZXVert input = boundary[0];
    ZXVert x = zx.neighbours(input)[0];
    PhasedGen x_gen = zx.get_vertex_ZXGen<PhasedGen>(x);
    REQUIRE(x_gen.get_type() == ZXType::XSpider);
    REQUIRE(x_gen.get_qtype() == QuantumType::Quantum);
    REQUIRE(x_gen.get_param() == 1);
    REQUIRE(test_equiv_expr_c(zx.get_scalar(), 1));
    REQUIRE(zx.n_vertices() == 3);
    REQUIRE(zx.n_wires() == 2);
    REQUIRE_NOTHROW(zx.check_validity());
  }
  GIVEN("Rx") {
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::Rx, 0.3, {0});
    ZXDiagram zx = circuit_to_zx(circ);
    ZXVertVec boundary = zx.get_boundary(ZXType::Input, QuantumType::Quantum);
    ZXVert input = boundary[0];
    ZXVert x = zx.neighbours(input)[0];
    PhasedGen x_gen = zx.get_vertex_ZXGen<PhasedGen>(x);
    REQUIRE(x_gen.get_type() == ZXType::XSpider);
    REQUIRE(x_gen.get_qtype() == QuantumType::Quantum);
    REQUIRE(x_gen.get_param() == 0.3);
    REQUIRE(test_equiv_expr_c(zx.get_scalar(), 1));
    REQUIRE(zx.n_vertices() == 3);
    REQUIRE(zx.n_wires() == 2);
    REQUIRE_NOTHROW(zx.check_validity());
  }
  GIVEN("Z") {
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::Z, {0});
    ZXDiagram zx = circuit_to_zx(circ);
    ZXVertVec boundary = zx.get_boundary(ZXType::Input, QuantumType::Quantum);
    ZXVert input = boundary[0];
    ZXVert z = zx.neighbours(input)[0];
    PhasedGen z_gen = zx.get_vertex_ZXGen<PhasedGen>(z);
    REQUIRE(z_gen.get_type() == ZXType::ZSpider);
    REQUIRE(z_gen.get_qtype() == QuantumType::Quantum);
    REQUIRE(z_gen.get_param() == 1);
    REQUIRE(test_equiv_expr_c(zx.get_scalar(), 1));
    REQUIRE(zx.n_vertices() == 3);
    REQUIRE(zx.n_wires() == 2);
    REQUIRE_NOTHROW(zx.check_validity());
  }
  GIVEN("Rz") {
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::Rz, 0.4, {0});
    ZXDiagram zx = circuit_to_zx(circ);
    ZXVertVec boundary = zx.get_boundary(ZXType::Input, QuantumType::Quantum);
    ZXVert input = boundary[0];
    ZXVert z = zx.neighbours(input)[0];
    PhasedGen z_gen = zx.get_vertex_ZXGen<PhasedGen>(z);
    REQUIRE(z_gen.get_type() == ZXType::ZSpider);
    REQUIRE(z_gen.get_qtype() == QuantumType::Quantum);
    REQUIRE(z_gen.get_param() == 0.4);
    REQUIRE(test_equiv_expr_c(zx.get_scalar(), 1));
    REQUIRE(zx.n_vertices() == 3);
    REQUIRE(zx.n_wires() == 2);
    REQUIRE_NOTHROW(zx.check_validity());
  }
  GIVEN("H") {
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::H, {0});
    ZXDiagram zx = circuit_to_zx(circ);
    ZXVertVec boundary = zx.get_boundary(ZXType::Input, QuantumType::Quantum);
    ZXVert input = boundary[0];
    ZXVert h = zx.neighbours(input)[0];
    ZXGen_ptr h_ptr = zx.get_vertex_ZXGen_ptr(h);
    REQUIRE(h_ptr->get_type() == ZXType::Hbox);
    REQUIRE(h_ptr->get_qtype() == QuantumType::Quantum);
    REQUIRE(test_equiv_expr_c(zx.get_scalar(), 0.5));
    REQUIRE(zx.n_vertices() == 3);
    REQUIRE(zx.n_wires() == 2);
    REQUIRE_NOTHROW(zx.check_validity());
  }
  GIVEN("CX") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    ZXDiagram zx = circuit_to_zx(circ);
    ZXVertVec boundary = zx.get_boundary(ZXType::Input, QuantumType::Quantum);
    ZXVert input0 = boundary[0];
    ZXVert input1 = boundary[1];
    ZXVert ctr = zx.neighbours(input0)[0];
    ZXVert targ = zx.neighbours(input1)[0];
    std::optional<Wire> w = zx.wire_between(ctr, targ);
    PhasedGen ctr_gen = zx.get_vertex_ZXGen<PhasedGen>(ctr);
    PhasedGen targ_gen = zx.get_vertex_ZXGen<PhasedGen>(targ);
    REQUIRE(ctr_gen.get_param() == 0);
    REQUIRE(targ_gen.get_param() == 0);
    REQUIRE(ctr_gen.get_type() == ZXType::ZSpider);
    REQUIRE(ctr_gen.get_qtype() == QuantumType::Quantum);
    REQUIRE(targ_gen.get_type() == ZXType::XSpider);
    REQUIRE(targ_gen.get_qtype() == QuantumType::Quantum);
    REQUIRE(test_equiv_expr_c(zx.get_scalar(), 2));
    REQUIRE(zx.n_vertices() == 6);
    REQUIRE(zx.n_wires() == 5);
    REQUIRE(zx.get_qtype(w.value()) == QuantumType::Quantum);
    REQUIRE(zx.get_wire_type(w.value()) == ZXWireType::Basic);
    REQUIRE_NOTHROW(zx.check_validity());
  }
  GIVEN("CZ") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CZ, {0, 1});
    ZXDiagram zx = circuit_to_zx(circ);
    ZXVertVec boundary = zx.get_boundary(ZXType::Input, QuantumType::Quantum);
    ZXVert input0 = boundary[0];
    ZXVert input1 = boundary[1];
    ZXVert ctr = zx.neighbours(input0)[0];
    ZXVert targ = zx.neighbours(input1)[0];
    std::optional<Wire> w = zx.wire_between(ctr, targ);
    PhasedGen ctr_gen = zx.get_vertex_ZXGen<PhasedGen>(ctr);
    PhasedGen targ_gen = zx.get_vertex_ZXGen<PhasedGen>(targ);
    REQUIRE(ctr_gen.get_param() == 0);
    REQUIRE(targ_gen.get_param() == 0);
    REQUIRE(ctr_gen.get_type() == ZXType::ZSpider);
    REQUIRE(ctr_gen.get_qtype() == QuantumType::Quantum);
    REQUIRE(targ_gen.get_type() == ZXType::ZSpider);
    REQUIRE(targ_gen.get_qtype() == QuantumType::Quantum);
    REQUIRE(test_equiv_expr_c(zx.get_scalar(), 1));
    REQUIRE(zx.n_vertices() == 6);
    REQUIRE(zx.n_wires() == 5);
    REQUIRE(zx.get_qtype(w.value()) == QuantumType::Quantum);
    REQUIRE(zx.get_wire_type(w.value()) == ZXWireType::H);
    REQUIRE_NOTHROW(zx.check_validity());
  }
  GIVEN("Measure") {
    Circuit circ(1, 1);
    circ.add_op<unsigned>(OpType::Measure, {0, 0});
    ZXDiagram zx = circuit_to_zx(circ);
    ZXVertVec q_boundary = zx.get_boundary(ZXType::Input, QuantumType::Quantum);
    ZXVertVec c_boundary =
        zx.get_boundary(ZXType::Input, QuantumType::Classical);
    ZXVert input0 = q_boundary[0];
    ZXVert input1 = c_boundary[0];
    ZXVert vert0 = zx.neighbours(input0)[0];
    ZXVert vert1 = zx.neighbours(input1)[0];
    PhasedGen gen0 = zx.get_vertex_ZXGen<PhasedGen>(vert0);
    PhasedGen gen1 = zx.get_vertex_ZXGen<PhasedGen>(vert1);
    REQUIRE(zx.degree(vert0) == 3);
    REQUIRE(zx.degree(vert1) == 1);
    REQUIRE(gen0.get_param() == 0);
    REQUIRE(gen1.get_param() == 0);
    REQUIRE(gen0.get_type() == ZXType::ZSpider);
    REQUIRE(gen0.get_qtype() == QuantumType::Classical);
    REQUIRE(gen1.get_type() == ZXType::ZSpider);
    REQUIRE(gen1.get_qtype() == QuantumType::Classical);
    REQUIRE(test_equiv_expr_c(zx.get_scalar(), 1));
    REQUIRE(zx.n_vertices() == 6);
    REQUIRE(zx.n_wires() == 4);
    REQUIRE_NOTHROW(zx.check_validity());
  }
  GIVEN("Reset") {
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::Reset, {0});
    ZXDiagram zx = circuit_to_zx(circ);
    ZXVertVec in_boundary =
        zx.get_boundary(ZXType::Input, QuantumType::Quantum);
    ZXVertVec out_boundary =
        zx.get_boundary(ZXType::Output, QuantumType::Quantum);
    ZXVert input = in_boundary[0];
    ZXVert output = out_boundary[0];
    ZXVert discard = zx.neighbours(input)[0];
    ZXVert init = zx.neighbours(output)[0];
    PhasedGen gen0 = zx.get_vertex_ZXGen<PhasedGen>(discard);
    PhasedGen gen1 = zx.get_vertex_ZXGen<PhasedGen>(init);
    REQUIRE(zx.degree(discard) == 1);
    REQUIRE(zx.degree(init) == 1);
    REQUIRE(gen0.get_param() == 0);
    REQUIRE(gen1.get_param() == 0);
    REQUIRE(gen0.get_type() == ZXType::ZSpider);
    REQUIRE(gen0.get_qtype() == QuantumType::Classical);
    REQUIRE(gen1.get_type() == ZXType::XSpider);
    REQUIRE(gen1.get_qtype() == QuantumType::Quantum);
    REQUIRE(test_equiv_expr_c(zx.get_scalar(), 0.5));
    REQUIRE(zx.n_vertices() == 4);
    REQUIRE(zx.n_wires() == 2);
    REQUIRE_NOTHROW(zx.check_validity());
  }
  GIVEN("Collapse") {
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::Collapse, {0});
    ZXDiagram zx = circuit_to_zx(circ);
    ZXVertVec boundary = zx.get_boundary(ZXType::Input, QuantumType::Quantum);
    ZXVert input = boundary[0];
    ZXVert z = zx.neighbours(input)[0];
    PhasedGen z_gen = zx.get_vertex_ZXGen<PhasedGen>(z);
    REQUIRE(z_gen.get_type() == ZXType::ZSpider);
    REQUIRE(z_gen.get_qtype() == QuantumType::Classical);
    REQUIRE(z_gen.get_param() == 0.);
    REQUIRE(test_equiv_expr_c(zx.get_scalar(), 1));
    REQUIRE(zx.n_vertices() == 3);
    REQUIRE(zx.n_wires() == 2);
    REQUIRE_NOTHROW(zx.check_validity());
  }
  GIVEN("Barrier") {
    Circuit circ(1, 1);
    circ.add_barrier({0}, {0});
    circ.add_barrier({0});
    circ.add_barrier({0}, {});
    ZXDiagram zx = circuit_to_zx(circ);
    ZXVertVec q_in_boundary =
        zx.get_boundary(ZXType::Input, QuantumType::Quantum);
    ZXVertVec c_in_boundary =
        zx.get_boundary(ZXType::Input, QuantumType::Classical);
    ZXVert q_input = q_in_boundary[0];
    ZXVert q_next = zx.neighbours(q_input)[0];
    ZXVert c_input = c_in_boundary[0];
    ZXVert c_next = zx.neighbours(c_input)[0];
    ZXVertVec q_out_boundary =
        zx.get_boundary(ZXType::Output, QuantumType::Quantum);
    ZXVertVec c_out_boundary =
        zx.get_boundary(ZXType::Output, QuantumType::Classical);
    REQUIRE(q_out_boundary[0] == q_next);
    REQUIRE(c_out_boundary[0] == c_next);
    REQUIRE(zx.n_vertices() == 4);
    REQUIRE(zx.n_wires() == 2);
    REQUIRE_NOTHROW(zx.check_validity());
  }
}

SCENARIO("Check converting circuits to diagrams") {
  GIVEN("A empty circuit") {
    Circuit circ;
    ZXDiagram zx = circuit_to_zx(circ);
    REQUIRE(zx.n_vertices() == 0);
    REQUIRE(zx.n_wires() == 0);
    REQUIRE(test_equiv_expr_c(zx.get_scalar(), 1));
    REQUIRE_NOTHROW(zx.check_validity());
  }

  GIVEN("A circuit with no gates") {
    Circuit circ(3, 1);
    Sym a = SymTable::fresh_symbol("a");
    Expr ea(a);
    circ.add_phase(ea);
    ZXDiagram zx = circuit_to_zx(circ);
    REQUIRE(zx.n_vertices() == 8);
    REQUIRE(zx.n_wires() == 4);
    REQUIRE(zx.count_vertices(ZXType::Input, QuantumType::Quantum) == 3);
    REQUIRE(zx.count_vertices(ZXType::Input, QuantumType::Classical) == 1);
    REQUIRE(zx.count_vertices(ZXType::Output, QuantumType::Quantum) == 3);
    REQUIRE(zx.count_vertices(ZXType::Output, QuantumType::Classical) == 1);
    REQUIRE(test_equiv_expr_c(zx.get_scalar(), 1));
    REQUIRE_NOTHROW(zx.check_validity());
  }

  GIVEN("A circuit with barriers") {
    Circuit circ(3, 2);
    circ.add_op<unsigned>(OpType::X, {2});
    circ.add_barrier({0, 2}, {0, 1});
    circ.add_barrier({0, 1});
    circ.add_barrier({}, {0});
    circ.add_op<unsigned>(OpType::X, {1});
    circ.add_op<unsigned>(OpType::Measure, {0, 0});
    circ.add_barrier({1}, {0, 1});
    ZXDiagram zx = circuit_to_zx(circ);
    REQUIRE(zx.n_vertices() == 14);
    REQUIRE(zx.n_wires() == 9);
    REQUIRE(zx.count_vertices(ZXType::Input, QuantumType::Quantum) == 3);
    REQUIRE(zx.count_vertices(ZXType::Input, QuantumType::Classical) == 2);
    REQUIRE(zx.count_vertices(ZXType::Output, QuantumType::Quantum) == 3);
    REQUIRE(zx.count_vertices(ZXType::Output, QuantumType::Classical) == 2);
    REQUIRE(zx.count_vertices(ZXType::XSpider, QuantumType::Quantum) == 2);
    REQUIRE(zx.count_vertices(ZXType::ZSpider, QuantumType::Classical) == 2);
    REQUIRE_NOTHROW(zx.check_validity());
  }

  GIVEN("A simple circuit") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.5, {1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::CZ, {1, 2});
    circ.add_op<unsigned>(OpType::CZ, {1, 0});
    ZXDiagram zx = circuit_to_zx(circ);
    REQUIRE(zx.n_vertices() == 15);
    REQUIRE(zx.count_vertices(ZXType::Input, QuantumType::Quantum) == 3);
    REQUIRE(zx.count_vertices(ZXType::Input, QuantumType::Classical) == 0);
    REQUIRE(zx.count_vertices(ZXType::Output, QuantumType::Quantum) == 3);
    REQUIRE(zx.count_vertices(ZXType::Output, QuantumType::Classical) == 0);
    REQUIRE(zx.count_vertices(ZXType::XSpider, QuantumType::Quantum) == 2);
    REQUIRE(zx.count_vertices(ZXType::ZSpider, QuantumType::Quantum) == 6);
    REQUIRE(zx.count_vertices(ZXType::Hbox, QuantumType::Quantum) == 1);
    REQUIRE(zx.count_wires(ZXWireType::H) == 2);
    REQUIRE(test_equiv_expr_c(zx.get_scalar(), 1));
    REQUIRE_NOTHROW(zx.check_validity());
  }

  GIVEN("A simple symbolic circuit") {
    Sym a = SymTable::fresh_symbol("a");
    Expr ea(a);
    Sym b = SymTable::fresh_symbol("b");
    Expr eb(b);
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::Rz, ea, {0});
    circ.add_phase(eb);
    ZXDiagram zx = circuit_to_zx(circ);
    REQUIRE(zx.n_vertices() == 3);
    REQUIRE(zx.count_vertices(ZXType::ZSpider, QuantumType::Quantum) == 1);
    REQUIRE(test_equiv_expr_c(zx.get_scalar(), 1));
    REQUIRE_NOTHROW(zx.check_validity());
  }

  GIVEN("A simple circuit with projective operations") {
    Circuit circ(3, 1);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Measure, {0, 0});
    circ.add_op<unsigned>(OpType::X, {2});
    circ.add_op<unsigned>(OpType::Reset, {2});
    circ.add_op<unsigned>(OpType::Collapse, {1});
    ZXDiagram zx = circuit_to_zx(circ);
    REQUIRE(zx.n_vertices() == 16);
    REQUIRE(zx.count_vertices(ZXType::Input, QuantumType::Quantum) == 3);
    REQUIRE(zx.count_vertices(ZXType::Input, QuantumType::Classical) == 1);
    REQUIRE(zx.count_vertices(ZXType::Output, QuantumType::Quantum) == 3);
    REQUIRE(zx.count_vertices(ZXType::Output, QuantumType::Classical) == 1);
    REQUIRE(zx.count_vertices(ZXType::XSpider, QuantumType::Quantum) == 3);
    REQUIRE(zx.count_vertices(ZXType::ZSpider, QuantumType::Quantum) == 1);
    REQUIRE(zx.count_vertices(ZXType::ZSpider, QuantumType::Classical) == 4);
    REQUIRE(test_equiv_expr_c(zx.get_scalar(), 1));
    REQUIRE_NOTHROW(zx.check_validity());
  }
}
}  // namespace test_ZXConverters
}  // namespace zx
}  // namespace tket