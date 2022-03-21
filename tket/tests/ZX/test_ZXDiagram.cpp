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
#include <fstream>

#include "../testutil.hpp"
#include "Gate/SymTable.hpp"
#include "ZX/ZXDiagram.hpp"
#include "ZX/ZXGenerator.hpp"

namespace tket {
namespace zx {
namespace test_ZXDiagram {

SCENARIO("Testing generator creation") {
  BoundaryGen input(ZXType::Input, QuantumType::Quantum);
  CHECK(input.get_name() == "Q-Input");
  CHECK(input.get_type() == ZXType::Input);
  CHECK(input.get_qtype() == QuantumType::Quantum);
  CHECK(input.valid_edge(std::nullopt, QuantumType::Quantum));
  CHECK_FALSE(input.valid_edge(0, QuantumType::Quantum));
  CHECK_FALSE(input.valid_edge(std::nullopt, QuantumType::Classical));

  PhasedGen zSpider(ZXType::ZSpider, 0.3, QuantumType::Classical);
  CHECK(zSpider.get_name() == "C-Z(0.3)");
  CHECK(zSpider.get_type() == ZXType::ZSpider);
  CHECK(zSpider.get_qtype() == QuantumType::Classical);
  CHECK(zSpider.free_symbols().empty());
  CHECK(zSpider.valid_edge(std::nullopt, QuantumType::Quantum));
  CHECK(zSpider.valid_edge(std::nullopt, QuantumType::Classical));
  CHECK_FALSE(zSpider.valid_edge(0, QuantumType::Quantum));

  PhasedGen xSpider(ZXType::XSpider, Expr("2*a"), QuantumType::Quantum);
  CHECK(xSpider.get_name() == "Q-X(2*a)");
  CHECK(xSpider.get_type() == ZXType::XSpider);
  CHECK(xSpider.get_qtype() == QuantumType::Quantum);
  CHECK(xSpider.free_symbols().size() == 1);
  CHECK(xSpider.valid_edge(std::nullopt, QuantumType::Quantum));
  CHECK_FALSE(xSpider.valid_edge(std::nullopt, QuantumType::Classical));
  SymEngine::map_basic_basic sub_map;
  Sym a = SymEngine::symbol("a");
  sub_map[a] = Expr(0.8);
  CHECK(xSpider.symbol_substitution(sub_map)->get_name() == "Q-X(1.6)");
  CHECK(
      *xSpider.symbol_substitution(sub_map) ==
      PhasedGen(ZXType::XSpider, 1.6, QuantumType::Quantum));

  CliffordGen px(ZXType::PX, true, QuantumType::Classical);
  CHECK(px.get_name() == "C-X(1)");
  CHECK(px.get_type() == ZXType::PX);
  CHECK(px.get_param() == true);
  CHECK(px.free_symbols().empty());
  CHECK(!(px == CliffordGen(ZXType::PX, false, QuantumType::Quantum)));
  CHECK(px == CliffordGen(ZXType::PX, true, QuantumType::Classical));

  // Should throw an error: type Triangle is not a BasicGen type
  REQUIRE_THROWS_AS(PhasedGen(ZXType::Triangle, 0.3), ZXError);

  DirectedGen tri(ZXType::Triangle, QuantumType::Classical);
  CHECK(tri.get_name() == "C-Tri");
  CHECK_FALSE(tri.valid_edge(std::nullopt, QuantumType::Classical));
  CHECK_FALSE(tri.valid_edge(2, QuantumType::Classical));
  CHECK(tri.valid_edge(0, QuantumType::Classical));
  CHECK_FALSE(tri.valid_edge(1, QuantumType::Quantum));
}

SCENARIO("Testing diagram creation & vertex/edge additions") {
  ZXDiagram diag(1, 1, 0, 0);
  CHECK(diag.get_scalar() == Expr(1.));
  CHECK_FALSE(diag.is_symbolic());
  diag.multiply_scalar(0.4);
  diag.multiply_scalar(Expr("2*a"));
  CHECK(diag.get_scalar() == Expr("0.8*a"));
  CHECK(diag.free_symbols().size() == 1);

  ZXVert zSpid_v = diag.add_vertex(ZXType::ZSpider, 0.1);
  ZXVert xSpid_v = diag.add_vertex(ZXType::XSpider, 3.4);
  ZXVert hbox_v = diag.add_vertex(
      ZXType::Hbox, 6.7 * Expr("b") + 3 * Expr(i_), QuantumType::Classical);

  REQUIRE_THROWS_AS(diag.add_vertex(ZXType::ZXBox, 3.), ZXError);

  REQUIRE_THROWS_WITH(
      diag.check_validity(), "Boundary vertex does not have degree 1");

  diag.add_wire(diag.get_boundary().at(0), zSpid_v);
  diag.add_wire(diag.get_boundary().at(1), xSpid_v);
  diag.add_wire(zSpid_v, xSpid_v);
  diag.add_wire(xSpid_v, zSpid_v, ZXWireType::H);
  Wire extra = diag.add_wire(diag.get_boundary().at(1), zSpid_v);
  REQUIRE_THROWS_WITH(
      diag.check_validity(), "Boundary vertex does not have degree 1");

  diag.remove_wire(extra);
  REQUIRE_NOTHROW(diag.check_validity());

  Wire wrong_port = diag.add_wire(
      hbox_v, xSpid_v, ZXWireType::Basic, QuantumType::Quantum, 0);
  REQUIRE_THROWS_WITH(
      diag.check_validity(), "Wire at a named port of an undirected vertex");
  diag.remove_wire(wrong_port);

  ZXVert tri_v = diag.add_vertex(ZXType::Triangle);
  diag.add_wire(tri_v, zSpid_v, ZXWireType::Basic, QuantumType::Quantum, 0);
  REQUIRE_THROWS_WITH(
      diag.check_validity(),
      "Not all ports of a directed vertex have wires connected");

  diag.add_wire(zSpid_v, tri_v);
  REQUIRE_THROWS_WITH(
      diag.check_validity(), "Wire at an unnamed port of a directed vertex");

  CHECK(diag.remove_wire(
      tri_v, zSpid_v,
      {ZXWireType::Basic, QuantumType::Quantum, std::nullopt, std::nullopt}));
  diag.add_wire(
      zSpid_v, tri_v, ZXWireType::Basic, QuantumType::Quantum, std::nullopt, 1);
  REQUIRE_NOTHROW(diag.check_validity());

  Wire extra_port =
      diag.add_wire(tri_v, zSpid_v, ZXWireType::Basic, QuantumType::Quantum, 1);
  REQUIRE_THROWS_WITH(
      diag.check_validity(), "Multiple wires on the same port of a vertex");
  diag.remove_wire(extra_port);

  ZXDiagram inner(1, 2, 1, 0);
  ZXVert inner_spid =
      inner.add_vertex(ZXType::ZSpider, 0.6, QuantumType::Classical);
  inner.add_wire(inner_spid, inner.get_boundary().at(0));
  inner.add_wire(inner_spid, inner.get_boundary().at(1));
  inner.add_wire(inner_spid, inner.get_boundary().at(2), ZXWireType::H);
  inner.add_wire(
      inner_spid, inner.get_boundary().at(3), ZXWireType::Basic,
      QuantumType::Classical);
  ZXGen_ptr box = std::make_shared<const ZXBox>(inner);

  ZXVert box_v = diag.add_vertex(box);
  diag.add_wire(box_v, hbox_v, ZXWireType::Basic, QuantumType::Quantum, 0);
  diag.add_wire(box_v, hbox_v, ZXWireType::Basic, QuantumType::Quantum, 1);
  diag.add_wire(box_v, xSpid_v, ZXWireType::Basic, QuantumType::Quantum, 2);
  Wire wrong_qtype =
      diag.add_wire(box_v, hbox_v, ZXWireType::Basic, QuantumType::Quantum, 3);
  REQUIRE_THROWS_WITH(
      diag.check_validity(),
      "QuantumType of wire is incompatible with the given port");

  diag.set_wire_qtype(wrong_qtype, QuantumType::Classical);
  REQUIRE_NOTHROW(diag.check_validity());

  CHECK(diag.free_symbols().size() == 2);
  SymEngine::map_basic_basic sub_map;
  Sym a = SymEngine::symbol("a");
  Sym b = SymEngine::symbol("b");
  sub_map[a] = Expr(0.8);
  diag.symbol_substitution(sub_map);
  CHECK(diag.free_symbols().size() == 1);
  sub_map[b] = Expr(0.4);
  diag.symbol_substitution(sub_map);
  CHECK(diag.free_symbols().size() == 0);
  CHECK(diag.get_name(hbox_v) == "C-H(2.68 + 3.0*I)");
  REQUIRE_NOTHROW(diag.check_validity());

  THEN("Print diagram to file") {
    std::ofstream dot_file("zxdiag.dot");
    dot_file << diag.to_graphviz_str();
    dot_file.close();
    remove("zxdiag.dot");
  }
}

SCENARIO("Test move constructors") {
  GIVEN("A single spider") {
    ZXDiagram diag(1, 2, 0, 0);
    ZXVertVec ins = diag.get_boundary(ZXType::Input);
    ZXVertVec outs = diag.get_boundary(ZXType::Output);
    ZXVert z = diag.add_vertex(ZXType::ZSpider, 0.3);
    diag.add_wire(ins[0], z);
    diag.add_wire(outs[0], z);
    diag.add_wire(outs[1], z);

    WHEN("Move constructor") {
      ZXDiagram d2(std::move(diag));
      REQUIRE_NOTHROW(d2.check_validity());
      CHECK(d2.n_vertices() == 4);
      CHECK(d2.n_wires() == 3);
    }
    WHEN("Move assignment") {
      ZXDiagram d2(4, 4, 1, 3);
      d2 = std::move(diag);
      REQUIRE_NOTHROW(d2.check_validity());
      CHECK(d2.n_vertices() == 4);
      CHECK(d2.n_wires() == 3);
    }
  }
}

SCENARIO("Check that diagram conversions achieve the correct form") {
  GIVEN("A mixed circuit") {
    ZXDiagram diag(2, 2, 1, 1);
    ZXVertVec ins = diag.get_boundary(ZXType::Input);
    ZXVertVec outs = diag.get_boundary(ZXType::Output);
    ZXVert qz = diag.add_vertex(ZXType::ZSpider, 0.3);
    ZXVert qx = diag.add_vertex(ZXType::XSpider);
    ZXVert cz = diag.add_vertex(ZXType::ZSpider, QuantumType::Classical);
    ZXDiagram inner(1, 1, 0, 0);
    ZXVert h = inner.add_vertex(ZXType::Hbox, Expr(i_));
    ZXVert tri = inner.add_vertex(ZXType::Triangle);
    ZXGen_ptr box = std::make_shared<const ZXBox>(inner);
    ZXVert b = diag.add_vertex(box);
    diag.add_wire(ins[0], qz);
    diag.add_wire(qz, outs[0]);
    diag.add_wire(qx, outs[1], ZXWireType::H);
    diag.add_wire(ins[1], cz, ZXWireType::Basic, QuantumType::Quantum);
    diag.add_wire(ins[2], cz, ZXWireType::H, QuantumType::Classical);
    diag.add_wire(outs[2], cz, ZXWireType::Basic, QuantumType::Classical);
    diag.add_wire(b, qz, ZXWireType::Basic, QuantumType::Quantum, 0);
    diag.add_wire(
        qx, b, ZXWireType::Basic, QuantumType::Quantum, std::nullopt, 1);
    THEN("Expand quantum vertices/edges into pairs of classical ones") {
      ZXDiagram doubled = diag.to_doubled_diagram();
      REQUIRE_NOTHROW(doubled.check_validity());
      CHECK(doubled.n_vertices() == 16);
      CHECK(doubled.n_wires() == 14);
      for (const ZXVert &b : doubled.get_boundary()) {
        CHECK(doubled.get_qtype(b) == QuantumType::Classical);
        CHECK(
            doubled.get_qtype(doubled.adj_wires(b)[0]) ==
            QuantumType::Classical);
        CHECK(
            doubled.get_qtype(doubled.neighbours(b)[0]) ==
            QuantumType::Classical);
      }
      ZXVertVec ins = doubled.get_boundary(ZXType::Input);
      CHECK(doubled.get_name(doubled.neighbours(ins[0])[0]) == "C-Z(0.3)");
      CHECK(doubled.get_name(doubled.neighbours(ins[1])[0]) == "C-Z(-0.3)");
    }
    GIVEN("Embedding classical boundaries into quantum states") {
      ZXDiagram embedded = diag.to_quantum_embedding();
      REQUIRE_NOTHROW(embedded.check_validity());
      CHECK(embedded.n_vertices() == 12);
      CHECK(embedded.n_wires() == 10);
      for (const ZXVert &b : embedded.get_boundary()) {
        CHECK(embedded.get_qtype(b) == QuantumType::Quantum);
      }
    }
  }
}

SCENARIO("Check converting gates to spiders") {
  GIVEN("X") {
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::X, {0});
    ZXDiagram zx(circ);
    ZXVertVec boundary = zx.get_boundary(ZXType::Input, QuantumType::Quantum);
    ZXVert input = boundary[0];
    ZXVert x = zx.neighbours(input)[0];
    ZXGen_ptr x_ptr = zx.get_vertex_ZXGen_ptr(x);
    REQUIRE(x_ptr->get_type() == ZXType::XSpider);
    REQUIRE(x_ptr->get_qtype() == QuantumType::Quantum);
    PhasedGen x_gen = static_cast<const PhasedGen &>(*x_ptr);
    REQUIRE(x_gen.get_param() == 1);
    REQUIRE(test_equiv_expr_c(zx.get_scalar(), 1));
    REQUIRE(zx.n_vertices() == 3);
    REQUIRE(zx.n_wires() == 2);
  }
  GIVEN("Rx") {
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::Rx, 0.3, {0});
    ZXDiagram zx(circ);
    ZXVertVec boundary = zx.get_boundary(ZXType::Input, QuantumType::Quantum);
    ZXVert input = boundary[0];
    ZXVert x = zx.neighbours(input)[0];
    ZXGen_ptr x_ptr = zx.get_vertex_ZXGen_ptr(x);
    REQUIRE(x_ptr->get_type() == ZXType::XSpider);
    REQUIRE(x_ptr->get_qtype() == QuantumType::Quantum);
    PhasedGen x_gen = static_cast<const PhasedGen &>(*x_ptr);
    REQUIRE(x_gen.get_param() == 0.3);
    REQUIRE(test_equiv_expr_c(zx.get_scalar(), exp(-i_ * PI * 0.3)));
    // False
    REQUIRE(zx.n_vertices() == 3);
    REQUIRE(zx.n_wires() == 2);
  }
  GIVEN("Z") {
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::Z, {0});
    ZXDiagram zx(circ);
    ZXVertVec boundary = zx.get_boundary(ZXType::Input, QuantumType::Quantum);
    ZXVert input = boundary[0];
    ZXVert z = zx.neighbours(input)[0];
    ZXGen_ptr z_ptr = zx.get_vertex_ZXGen_ptr(z);
    REQUIRE(z_ptr->get_type() == ZXType::ZSpider);
    REQUIRE(z_ptr->get_qtype() == QuantumType::Quantum);
    PhasedGen z_gen = static_cast<const PhasedGen &>(*z_ptr);
    REQUIRE(z_gen.get_param() == 1);
    REQUIRE(test_equiv_expr_c(zx.get_scalar(), 1));
    REQUIRE(zx.n_vertices() == 3);
    REQUIRE(zx.n_wires() == 2);
  }
  GIVEN("Rz") {
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::Rz, 0.4, {0});
    ZXDiagram zx(circ);
    ZXVertVec boundary = zx.get_boundary(ZXType::Input, QuantumType::Quantum);
    ZXVert input = boundary[0];
    ZXVert z = zx.neighbours(input)[0];
    ZXGen_ptr z_ptr = zx.get_vertex_ZXGen_ptr(z);
    REQUIRE(z_ptr->get_type() == ZXType::ZSpider);
    REQUIRE(z_ptr->get_qtype() == QuantumType::Quantum);
    PhasedGen z_gen = static_cast<const PhasedGen &>(*z_ptr);
    REQUIRE(z_gen.get_param() == 0.4);
    REQUIRE(test_equiv_expr_c(zx.get_scalar(), exp(-i_ * PI * 0.4)));
    REQUIRE(zx.n_vertices() == 3);
    REQUIRE(zx.n_wires() == 2);
  }
  GIVEN("H") {
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::H, {0});
    ZXDiagram zx(circ);
    ZXVertVec boundary = zx.get_boundary(ZXType::Input, QuantumType::Quantum);
    ZXVert input = boundary[0];
    ZXVert h = zx.neighbours(input)[0];
    ZXGen_ptr h_ptr = zx.get_vertex_ZXGen_ptr(h);
    REQUIRE(h_ptr->get_type() == ZXType::Hbox);
    REQUIRE(h_ptr->get_qtype() == QuantumType::Quantum);
    REQUIRE(test_equiv_expr_c(zx.get_scalar(), 0.5));
    REQUIRE(zx.n_vertices() == 3);
    REQUIRE(zx.n_wires() == 2);
  }
  GIVEN("CX") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    ZXDiagram zx(circ);
    ZXVertVec boundary = zx.get_boundary(ZXType::Input, QuantumType::Quantum);
    ZXVert input0 = boundary[0];
    ZXVert input1 = boundary[1];
    ZXVert ctr = zx.neighbours(input0)[0];
    ZXVert targ = zx.neighbours(input1)[0];
    ZXGen_ptr ctr_ptr = zx.get_vertex_ZXGen_ptr(ctr);
    ZXGen_ptr targ_ptr = zx.get_vertex_ZXGen_ptr(targ);
    std::optional<Wire> w = zx.wire_between(ctr, targ);
    PhasedGen ctr_gen = static_cast<const PhasedGen &>(*ctr_ptr);
    REQUIRE(ctr_gen.get_param() == 0);
    PhasedGen targ_gen = static_cast<const PhasedGen &>(*targ_ptr);
    REQUIRE(targ_gen.get_param() == 0);
    REQUIRE(ctr_ptr->get_type() == ZXType::ZSpider);
    REQUIRE(ctr_ptr->get_qtype() == QuantumType::Quantum);
    REQUIRE(targ_ptr->get_type() == ZXType::XSpider);
    REQUIRE(targ_ptr->get_qtype() == QuantumType::Quantum);
    REQUIRE(test_equiv_expr_c(zx.get_scalar(), 2));
    REQUIRE(zx.n_vertices() == 6);
    REQUIRE(zx.n_wires() == 5);
    REQUIRE(zx.get_qtype(w.value()) == QuantumType::Quantum);
    REQUIRE(zx.get_wire_type(w.value()) == ZXWireType::Basic);
  }
  GIVEN("CZ") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CZ, {0, 1});
    ZXDiagram zx(circ);
    ZXVertVec boundary = zx.get_boundary(ZXType::Input, QuantumType::Quantum);
    ZXVert input0 = boundary[0];
    ZXVert input1 = boundary[1];
    ZXVert ctr = zx.neighbours(input0)[0];
    ZXVert targ = zx.neighbours(input1)[0];
    ZXGen_ptr ctr_ptr = zx.get_vertex_ZXGen_ptr(ctr);
    ZXGen_ptr targ_ptr = zx.get_vertex_ZXGen_ptr(targ);
    std::optional<Wire> w = zx.wire_between(ctr, targ);
    PhasedGen ctr_gen = static_cast<const PhasedGen &>(*ctr_ptr);
    REQUIRE(ctr_gen.get_param() == 0);
    PhasedGen targ_gen = static_cast<const PhasedGen &>(*targ_ptr);
    REQUIRE(targ_gen.get_param() == 0);
    REQUIRE(ctr_ptr->get_type() == ZXType::ZSpider);
    REQUIRE(ctr_ptr->get_qtype() == QuantumType::Quantum);
    REQUIRE(targ_ptr->get_type() == ZXType::ZSpider);
    REQUIRE(targ_ptr->get_qtype() == QuantumType::Quantum);
    REQUIRE(test_equiv_expr_c(zx.get_scalar(), 1));
    REQUIRE(zx.n_vertices() == 6);
    REQUIRE(zx.n_wires() == 5);
    REQUIRE(zx.get_qtype(w.value()) == QuantumType::Quantum);
    REQUIRE(zx.get_wire_type(w.value()) == ZXWireType::H);
  }
  GIVEN("Measure") {
    Circuit circ(1, 1);
    circ.add_op<unsigned>(OpType::Measure, {0, 0});
    ZXDiagram zx(circ);
    ZXVertVec q_boundary = zx.get_boundary(ZXType::Input, QuantumType::Quantum);
    ZXVertVec c_boundary =
        zx.get_boundary(ZXType::Input, QuantumType::Classical);
    ZXVert input0 = q_boundary[0];
    ZXVert input1 = c_boundary[0];
    ZXVert vert0 = zx.neighbours(input0)[0];
    ZXVert vert1 = zx.neighbours(input1)[0];
    ZXGen_ptr gen_ptr0 = zx.get_vertex_ZXGen_ptr(vert0);
    ZXGen_ptr gen_ptr1 = zx.get_vertex_ZXGen_ptr(vert1);
    PhasedGen gen0 = static_cast<const PhasedGen &>(*gen_ptr0);
    REQUIRE(gen0.get_param() == 0);
    PhasedGen gen1 = static_cast<const PhasedGen &>(*gen_ptr1);
    REQUIRE(gen1.get_param() == 0);
    REQUIRE(gen_ptr0->get_type() == ZXType::ZSpider);
    REQUIRE(gen_ptr0->get_qtype() == QuantumType::Classical);
    REQUIRE(gen_ptr1->get_type() == ZXType::ZSpider);
    REQUIRE(gen_ptr1->get_qtype() == QuantumType::Classical);
    REQUIRE(test_equiv_expr_c(zx.get_scalar(), 1));
    REQUIRE(zx.n_vertices() == 6);
    REQUIRE(zx.n_wires() == 4);
  }
  GIVEN("Reset") {
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::Reset, {0});
    ZXDiagram zx(circ);
    ZXVertVec in_boundary =
        zx.get_boundary(ZXType::Input, QuantumType::Quantum);
    ZXVertVec out_boundary =
        zx.get_boundary(ZXType::Output, QuantumType::Quantum);
    ZXVert input = in_boundary[0];
    ZXVert output = out_boundary[0];
    ZXVert discard = zx.neighbours(input)[0];
    ZXVert init = zx.neighbours(output)[0];
    ZXGen_ptr gen_ptr0 = zx.get_vertex_ZXGen_ptr(discard);
    ZXGen_ptr gen_ptr1 = zx.get_vertex_ZXGen_ptr(init);
    PhasedGen gen0 = static_cast<const PhasedGen &>(*gen_ptr0);
    REQUIRE(gen0.get_param() == 0);
    PhasedGen gen1 = static_cast<const PhasedGen &>(*gen_ptr1);
    REQUIRE(gen1.get_param() == 0);
    REQUIRE(gen_ptr0->get_type() == ZXType::ZSpider);
    REQUIRE(gen_ptr0->get_qtype() == QuantumType::Classical);
    REQUIRE(gen_ptr1->get_type() == ZXType::XSpider);
    REQUIRE(gen_ptr1->get_qtype() == QuantumType::Quantum);
    REQUIRE(test_equiv_expr_c(zx.get_scalar(), 0.5));
    REQUIRE(zx.n_vertices() == 4);
    REQUIRE(zx.n_wires() == 2);
  }
  GIVEN("Collapse") {
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::Collapse, {0});
    ZXDiagram zx(circ);
    ZXVertVec boundary = zx.get_boundary(ZXType::Input, QuantumType::Quantum);
    ZXVert input = boundary[0];
    ZXVert z = zx.neighbours(input)[0];
    ZXGen_ptr z_ptr = zx.get_vertex_ZXGen_ptr(z);
    REQUIRE(z_ptr->get_type() == ZXType::ZSpider);
    REQUIRE(z_ptr->get_qtype() == QuantumType::Classical);
    PhasedGen z_gen = static_cast<const PhasedGen &>(*z_ptr);
    REQUIRE(z_gen.get_param() == 0.);
    REQUIRE(test_equiv_expr_c(zx.get_scalar(), 1));
    REQUIRE(zx.n_vertices() == 3);
    REQUIRE(zx.n_wires() == 2);
  }
}

SCENARIO("Check converting circuits to diagrams") {
  GIVEN("A empty circuit") {
    Circuit circ;
    ZXDiagram zx(circ);
    REQUIRE(zx.n_vertices() == 0);
    REQUIRE(zx.n_wires() == 0);
    REQUIRE(test_equiv_expr_c(zx.get_scalar(), 1));
  }

  GIVEN("A circuit with no gates") {
    Circuit circ(3, 1);
    Sym a = SymTable::fresh_symbol("a");
    Expr ea(a);
    circ.add_phase(ea);
    ZXDiagram zx(circ);
    REQUIRE(zx.n_vertices() == 8);
    REQUIRE(zx.n_wires() == 4);
    REQUIRE(zx.count_vertices(ZXType::Input, QuantumType::Quantum) == 3);
    REQUIRE(zx.count_vertices(ZXType::Input, QuantumType::Classical) == 1);
    REQUIRE(zx.count_vertices(ZXType::Output, QuantumType::Quantum) == 3);
    REQUIRE(zx.count_vertices(ZXType::Output, QuantumType::Classical) == 1);
    Expr e = exp(2. * i_ * PI * ea);
    REQUIRE(zx.get_scalar() == 1.0 * e);
  }

  GIVEN("A simple circuit") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.5, {1});
    circ.add_op<unsigned>(OpType::H, {2});
    circ.add_op<unsigned>(OpType::CZ, {1, 2});
    circ.add_op<unsigned>(OpType::CZ, {1, 0});
    ZXDiagram zx(circ);
    REQUIRE(zx.n_vertices() == 15);
    REQUIRE(zx.count_vertices(ZXType::Input, QuantumType::Quantum) == 3);
    REQUIRE(zx.count_vertices(ZXType::Input, QuantumType::Classical) == 0);
    REQUIRE(zx.count_vertices(ZXType::Output, QuantumType::Quantum) == 3);
    REQUIRE(zx.count_vertices(ZXType::Output, QuantumType::Classical) == 0);
    REQUIRE(zx.count_vertices(ZXType::XSpider, QuantumType::Quantum) == 2);
    REQUIRE(zx.count_vertices(ZXType::ZSpider, QuantumType::Quantum) == 6);
    REQUIRE(zx.count_vertices(ZXType::Hbox, QuantumType::Quantum) == 1);
    REQUIRE(zx.count_wires(ZXWireType::H) == 2);
    REQUIRE(test_equiv_expr_c(zx.get_scalar(), exp(-i_ * PI * 0.5)));
  }

  GIVEN("A simple symbolic circuit") {
    Sym a = SymTable::fresh_symbol("a");
    Expr ea(a);
    Sym b = SymTable::fresh_symbol("b");
    Expr eb(b);
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::Rz, ea, {0});
    circ.add_phase(eb);
    ZXDiagram zx(circ);
    REQUIRE(zx.n_vertices() == 3);
    REQUIRE(zx.count_vertices(ZXType::ZSpider, QuantumType::Quantum) == 1);
    Expr e = exp(-i_ * PI * ea + 2. * i_ * PI * eb);
    REQUIRE(zx.get_scalar() == 1.0 * e);
  }

  GIVEN("A circuit with a simple box") {
    Circuit c0(2);
    c0.add_op<unsigned>(OpType::Rx, 0.5, {0});
    c0.add_op<unsigned>(OpType::Rz, 0.3, {1});
    c0.add_op<unsigned>(OpType::CX, {0, 1});
    CircBox c0box(c0);
    Circuit circ(2);
    circ.add_box(c0box, {0, 1});
    ZXDiagram zx(circ);
    REQUIRE(zx.n_vertices() == 8);
    REQUIRE(test_equiv_expr_c(zx.get_scalar(), 2. * exp(-i_ * PI * 0.8)));
  }

  GIVEN("A circuit requires rebase") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CRz, 0.3, {0, 1});
    ZXDiagram zx(circ);
    REQUIRE(zx.n_vertices() == 10);
    REQUIRE(zx.count_vertices(ZXType::Input, QuantumType::Quantum) == 2);
    REQUIRE(zx.count_vertices(ZXType::Output, QuantumType::Quantum) == 2);
    REQUIRE(zx.count_vertices(ZXType::XSpider, QuantumType::Quantum) == 2);
    REQUIRE(zx.count_vertices(ZXType::ZSpider, QuantumType::Quantum) == 4);
    REQUIRE(test_equiv_expr_c(zx.get_scalar(), 4));
  }

  GIVEN("A simple circuit with projective operations") {
    Circuit circ(3, 1);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Measure, {0, 0});
    circ.add_op<unsigned>(OpType::X, {2});
    circ.add_op<unsigned>(OpType::Reset, {2});
    circ.add_op<unsigned>(OpType::Collapse, {1});
    ZXDiagram zx(circ);
    REQUIRE(zx.n_vertices() == 16);
    REQUIRE(zx.count_vertices(ZXType::Input, QuantumType::Quantum) == 3);
    REQUIRE(zx.count_vertices(ZXType::Input, QuantumType::Classical) == 1);
    REQUIRE(zx.count_vertices(ZXType::Output, QuantumType::Quantum) == 3);
    REQUIRE(zx.count_vertices(ZXType::Output, QuantumType::Classical) == 1);
    REQUIRE(zx.count_vertices(ZXType::XSpider, QuantumType::Quantum) == 3);
    REQUIRE(zx.count_vertices(ZXType::ZSpider, QuantumType::Quantum) == 1);
    REQUIRE(zx.count_vertices(ZXType::ZSpider, QuantumType::Classical) == 4);
    REQUIRE(test_equiv_expr_c(zx.get_scalar(), 1));
  }
}

}  // namespace test_ZXDiagram
}  // namespace zx
}  // namespace tket
