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

#include <boost/graph/adjacency_list.hpp>
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
    circ.add_barrier({}, {0});
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
  GIVEN("noop") {
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::noop, {0});
    ZXDiagram zx = circuit_to_zx(circ);
    ZXVertVec q_in_boundary =
        zx.get_boundary(ZXType::Input, QuantumType::Quantum);
    ZXVert q_input = q_in_boundary[0];
    ZXVert q_next = zx.neighbours(q_input)[0];
    ZXVertVec q_out_boundary =
        zx.get_boundary(ZXType::Output, QuantumType::Quantum);
    REQUIRE(q_out_boundary[0] == q_next);
    REQUIRE(zx.n_vertices() == 2);
    REQUIRE(zx.n_wires() == 1);
    REQUIRE_NOTHROW(zx.check_validity());
  }
  GIVEN("SWAP") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::SWAP, {0, 1});
    circ.add_op<unsigned>(OpType::SWAP, {1, 2});
    circ.add_op<unsigned>(OpType::Rz, 0.1, {0});
    circ.add_op<unsigned>(OpType::Rz, 0.2, {1});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {2});

    ZXDiagram zx = circuit_to_zx(circ);
    ZXVertVec q_in_boundary =
        zx.get_boundary(ZXType::Input, QuantumType::Quantum);
    ZXVert q_input_0 = q_in_boundary[0];
    ZXVert q_next_0 = zx.neighbours(q_input_0)[0];
    ZXVert q_input_1 = q_in_boundary[1];
    ZXVert q_next_1 = zx.neighbours(q_input_1)[0];
    ZXVert q_input_2 = q_in_boundary[2];
    ZXVert q_next_2 = zx.neighbours(q_input_2)[0];

    PhasedGen z0_gen = zx.get_vertex_ZXGen<PhasedGen>(q_next_0);
    PhasedGen z1_gen = zx.get_vertex_ZXGen<PhasedGen>(q_next_1);
    PhasedGen z2_gen = zx.get_vertex_ZXGen<PhasedGen>(q_next_2);

    REQUIRE(z0_gen.get_type() == ZXType::ZSpider);
    REQUIRE(z0_gen.get_qtype() == QuantumType::Quantum);
    REQUIRE(z0_gen.get_param() == 0.3);

    REQUIRE(z1_gen.get_type() == ZXType::ZSpider);
    REQUIRE(z1_gen.get_qtype() == QuantumType::Quantum);
    REQUIRE(z1_gen.get_param() == 0.1);

    REQUIRE(z2_gen.get_type() == ZXType::ZSpider);
    REQUIRE(z2_gen.get_qtype() == QuantumType::Quantum);
    REQUIRE(z2_gen.get_param() == 0.2);

    REQUIRE(zx.n_vertices() == 9);
    REQUIRE(zx.n_wires() == 6);
    REQUIRE_NOTHROW(zx.check_validity());
  }
  GIVEN("Create") {
    Circuit circ(1);
    circ.qubit_create(Qubit(0));
    ZXDiagram zx = circuit_to_zx(circ);
    ZXVertVec input_boundary =
        zx.get_boundary(ZXType::Input, QuantumType::Quantum);
    REQUIRE(input_boundary.size() == 0);
    ZXVertVec output_boundary =
        zx.get_boundary(ZXType::Output, QuantumType::Quantum);
    ZXVert output = output_boundary[0];
    ZXVert x = zx.neighbours(output)[0];
    PhasedGen x_gen = zx.get_vertex_ZXGen<PhasedGen>(x);
    REQUIRE(x_gen.get_type() == ZXType::XSpider);
    REQUIRE(x_gen.get_qtype() == QuantumType::Quantum);
    REQUIRE(x_gen.get_param() == 0);
    REQUIRE(test_equiv_expr_c(zx.get_scalar(), 0.5));
    REQUIRE(zx.n_vertices() == 2);
    REQUIRE(zx.n_wires() == 1);
    REQUIRE_NOTHROW(zx.check_validity());
  }
  GIVEN("Discard") {
    Circuit circ(1);
    circ.qubit_discard(Qubit(0));
    ZXDiagram zx = circuit_to_zx(circ);
    ZXVertVec input_boundary =
        zx.get_boundary(ZXType::Input, QuantumType::Quantum);
    ZXVertVec output_boundary =
        zx.get_boundary(ZXType::Output, QuantumType::Quantum);
    REQUIRE(output_boundary.size() == 0);
    ZXVert input = input_boundary[0];
    ZXVert z = zx.neighbours(input)[0];
    PhasedGen z_gen = zx.get_vertex_ZXGen<PhasedGen>(z);
    REQUIRE(z_gen.get_type() == ZXType::ZSpider);
    REQUIRE(z_gen.get_qtype() == QuantumType::Classical);
    REQUIRE(z_gen.get_param() == 0);
    REQUIRE(test_equiv_expr_c(zx.get_scalar(), 1));
    REQUIRE(zx.n_vertices() == 2);
    REQUIRE(zx.n_wires() == 1);
    REQUIRE_NOTHROW(zx.check_validity());
  }
  GIVEN("CircBox") {
    Circuit circ(4);
    Circuit inner(3);
    inner.add_op<unsigned>(OpType::Rx, 0.5, {2});
    inner.add_op<unsigned>(OpType::CX, {0, 1});
    CircBox cbox(inner);
    circ.add_box(cbox, {0, 1, 2});
    ZXDiagram zx = circuit_to_zx(circ);
  }
  // GIVEN("Conditional") {
  //   Circuit circ(1,1);
  //   circ.add_conditional_gate<unsigned>(OpType::X, {}, {0}, {0}, 1);
  //   ZXDiagram zx = circuit_to_zx(circ);
  // }
  GIVEN("Conditional test") {
    Circuit circ(1,1);
    // Vertex m = circ.add_op<unsigned>(OpType::Measure, {0, 0});
    Vertex c = circ.add_conditional_gate<unsigned>(OpType::X, {}, {0}, {0}, 1);
    // EdgeVec ins =
    //     circ.get_in_edges(m);
    // std::cout<<"\nIn total\n"<<ins.size()<<"\n";
    // EdgeVec outs =
    //     circ.get_all_out_edges(m);
    // std::cout<<"\nOut total\n"<<outs.size()<<"\n";
    // EdgeVec q_in_holes =
    //     circ.get_in_edges_of_type(m, EdgeType::Quantum);
    // std::cout<<"\nq_in_holes\n"<<q_in_holes.size()<<"\n";
    // EdgeVec c_in_holes =
    //     circ.get_in_edges_of_type(m, EdgeType::Classical);
    // std::cout<<"\nc_in_holes\n"<<c_in_holes.size()<<"\n";
    // EdgeVec b_in_holes =
    //     circ.get_in_edges_of_type(m, EdgeType::Boolean);
    // std::cout<<"\nb_in_holes\n"<<b_in_holes.size()<<"\n";
    // EdgeVec b_out_holes =
    //     circ.get_out_edges_of_type(m, EdgeType::Boolean);
    // std::cout<<"\nb_out_holes\n"<<b_out_holes.size()<<"\n";
    // EdgeVec c_out_holes =
    //     circ.get_out_edges_of_type(m, EdgeType::Classical);
    // std::cout<<"\nc_out_holes\n"<<c_out_holes.size()<<"\n";
    EdgeVec ins2 =
        circ.get_in_edges(c);
    std::cout<<"\nIn total\n"<<ins2.size()<<"\n";
    EdgeVec outs2 =
        circ.get_all_out_edges(c);
    std::cout<<"\nOut total\n"<<outs2.size()<<"\n";
    EdgeVec q_in_holes2 =
        circ.get_in_edges_of_type(c, EdgeType::Quantum);
    std::cout<<"\nq_in_holes2\n"<<q_in_holes2.size()<<"\n";
    EdgeVec c_in_holes2 =
        circ.get_in_edges_of_type(c, EdgeType::Classical);
    std::cout<<"\nc_in_holes2\n"<<c_in_holes2.size()<<"\n";
    EdgeVec b_in_holes2 =
        circ.get_in_edges_of_type(c, EdgeType::Boolean);
    std::cout<<"\nb_in_holes2\n"<<b_in_holes2.size()<<"\n";
    EdgeVec b_out_holes2 =
        circ.get_out_edges_of_type(c, EdgeType::Boolean);
    std::cout<<"\nb_out_holes2\n"<<b_out_holes2.size()<<"\n";
    EdgeVec c_out_holes2 =
        circ.get_out_edges_of_type(c, EdgeType::Classical);
    std::cout<<"\nc_out_holes2\n"<<c_out_holes2.size()<<"\n";

    // 
    Vertex in = circ.get_in(Bit(0));
    EdgeVec outs3 =
    circ.get_all_out_edges(in);
    std::cout<<"\nOut total3\n"<<outs3.size()<<"\n";
    Op_ptr op = circ.get_Op_ptr_from_Vertex(in);
    std::cout<<"\nSig\n"<<op->get_signature().size()<<"\n";
    Op_ptr op2 = circ.get_Op_ptr_from_Vertex(c);
    std::cout<<"\nSig2\n"<<op2->get_signature().size()<<"\n";
    for (auto e : outs3) {
      std::cout<<"\nSource port" << circ.get_source_port(e);
    }
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

  GIVEN("A circuit with spiderless ops") {
    Circuit circ(3, 2);
    circ.add_op<unsigned>(OpType::X, {2});
    circ.add_barrier({0, 2}, {0, 1});
    circ.add_op<unsigned>(OpType::noop, {0});
    circ.add_op<unsigned>(OpType::noop, {0});
    circ.add_barrier({0, 1});
    circ.add_op<unsigned>(OpType::noop, {1});
    circ.add_barrier({}, {0});
    circ.add_op<unsigned>(OpType::SWAP, {0, 1});
    circ.add_op<unsigned>(OpType::SWAP, {0, 1});
    circ.add_op<unsigned>(OpType::SWAP, {0, 1});
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