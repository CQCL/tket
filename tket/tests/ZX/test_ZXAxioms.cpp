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

#include "ZX/Rewrite.hpp"

namespace tket {
namespace zx {
namespace test_ZXAxioms {

SCENARIO("Simplify to a graph-like diagram") {
  GIVEN("A manually-constructed diagram") {
    /**
     * Diagram on https://arxiv.org/pdf/1902.03178.pdf, Figure 2
     * We have added an extra input / output pair for testing purposes.
     */
    ZXDiagram diag1(5, 5, 0, 0);
    ZXVertVec diag1_inputs = diag1.get_boundary(ZXType::Input);
    ZXVertVec diag1_outputs = diag1.get_boundary(ZXType::Output);

    ZXVert zSpid1 = diag1.add_vertex(ZXType::ZSpider);
    ZXVert zSpid2 = diag1.add_vertex(ZXType::ZSpider);
    ZXVert zSpid3 = diag1.add_vertex(ZXType::ZSpider);
    ZXVert phZSpid1 = diag1.add_vertex(ZXType::ZSpider, 0.5);
    ZXVert phZSpid2 = diag1.add_vertex(ZXType::ZSpider, 1.);
    ZXVert xSpid1 = diag1.add_vertex(ZXType::XSpider);
    ZXVert xSpid2 = diag1.add_vertex(ZXType::XSpider);
    ZXVert xSpid3 = diag1.add_vertex(ZXType::XSpider);

    Wire w11 = diag1.add_wire(diag1_inputs[0], zSpid1);
    Wire w12 = diag1.add_wire(zSpid1, phZSpid1);
    Wire w13 = diag1.add_wire(phZSpid1, zSpid2);
    Wire w14 = diag1.add_wire(zSpid2, diag1_outputs[0], ZXWireType::H);
    Wire w21 = diag1.add_wire(zSpid1, xSpid1);
    Wire w22 = diag1.add_wire(zSpid2, xSpid2);
    Wire w31 = diag1.add_wire(diag1_inputs[1], xSpid1, ZXWireType::H);
    Wire w32 = diag1.add_wire(xSpid1, zSpid3);
    Wire w33 = diag1.add_wire(zSpid3, xSpid2);
    Wire w34 = diag1.add_wire(xSpid2, phZSpid2);
    Wire w35 = diag1.add_wire(phZSpid2, diag1_outputs[1]);
    Wire w41 = diag1.add_wire(zSpid3, xSpid3);
    Wire w51 = diag1.add_wire(diag1_inputs[2], xSpid3, ZXWireType::H);
    Wire w52 = diag1.add_wire(xSpid3, diag1_outputs[2]);
    Wire w61 = diag1.add_wire(diag1_inputs[3], diag1_outputs[3], ZXWireType::H);
    Wire w81 =
        diag1.add_wire(diag1_inputs[4], diag1_outputs[4], ZXWireType::Basic);

    REQUIRE_NOTHROW(diag1.check_validity());

    unsigned n_matches;

    // Replace X with Z spiders
    CHECK(Rewrite::red_to_green().apply(diag1));
    REQUIRE(diag1.count_vertices(ZXType::XSpider) == 0);
    REQUIRE(diag1.count_vertices(ZXType::ZSpider) == 8);

    // Spider fusion
    CHECK(Rewrite::spider_fusion().apply(diag1));
    REQUIRE(diag1.count_vertices(ZXType::ZSpider) == 6);

    // Parallel edge pair removal
    CHECK_FALSE(Rewrite::parallel_h_removal().apply(diag1));

    // Remove hadamard edges connected directly to the boundaries
    CHECK(Rewrite::io_extension().apply(diag1));
    REQUIRE(diag1.count_vertices(ZXType::ZSpider) == 10);

    // Boundary vertices sharing spiders
    // Deal with directly connected in/outputs
    CHECK(Rewrite::separate_boundaries().apply(diag1));
    REQUIRE(diag1.count_vertices(ZXType::ZSpider) == 13);

    REQUIRE_NOTHROW(diag1.check_validity());
  }
}

SCENARIO("Testing Spider fusion") {
  GIVEN("A manually-constructed diagram") {
    ZXDiagram diag2(2, 1, 0, 0);
    ZXVertVec diag2_inputs = diag2.get_boundary(ZXType::Input);
    ZXVertVec diag2_outputs = diag2.get_boundary(ZXType::Output);

    ZXVert spid1 = diag2.add_vertex(ZXType::ZSpider, 0.1);
    ZXVert spid2 = diag2.add_vertex(ZXType::ZSpider, 0.3);
    ZXVert spid3 = diag2.add_vertex(ZXType::ZSpider);
    ZXVert spid4 = diag2.add_vertex(ZXType::ZSpider, 0.5);
    ZXVert spid5 = diag2.add_vertex(ZXType::ZSpider);

    Wire i1s1 = diag2.add_wire(
        diag2_inputs[0], spid1, ZXWireType::Basic, QuantumType::Quantum);
    Wire i2s5 = diag2.add_wire(
        diag2_inputs[1], spid5, ZXWireType::H, QuantumType::Quantum);
    Wire s1s2 =
        diag2.add_wire(spid1, spid2, ZXWireType::H, QuantumType::Quantum);
    Wire s2s3 =
        diag2.add_wire(spid2, spid3, ZXWireType::Basic, QuantumType::Quantum);
    Wire s3s2 =
        diag2.add_wire(spid3, spid2, ZXWireType::H, QuantumType::Quantum);
    Wire s3s4 =
        diag2.add_wire(spid3, spid4, ZXWireType::H, QuantumType::Quantum);
    Wire s4s5 =
        diag2.add_wire(spid4, spid5, ZXWireType::Basic, QuantumType::Quantum);
    Wire s5s1 =
        diag2.add_wire(spid5, spid1, ZXWireType::Basic, QuantumType::Quantum);
    Wire s3o1 = diag2.add_wire(
        spid3, diag2_outputs[0], ZXWireType::Basic, QuantumType::Quantum);
    Wire loops3s3 =
        diag2.add_wire(spid3, spid3, ZXWireType::Basic, QuantumType::Quantum);
    Wire hloops3s3 =
        diag2.add_wire(spid3, spid3, ZXWireType::H, QuantumType::Quantum);

    REQUIRE_NOTHROW(diag2.check_validity());

    unsigned n_matches;

    // Remove self-loops
    CHECK(Rewrite::self_loop_removal().apply(diag2));

    // Spider fusion
    CHECK(Rewrite::spider_fusion().apply(diag2));
    REQUIRE(diag2.count_vertices(ZXType::ZSpider) == 2);

    // Remove self-loops after fusion
    CHECK(Rewrite::self_loop_removal().apply(diag2));

    // Parallel edge pair removal
    CHECK(Rewrite::parallel_h_removal().apply(diag2));

    // Remove hadamard edges connected directly to the boundaries
    CHECK(Rewrite::io_extension().apply(diag2));

    REQUIRE_NOTHROW(diag2.check_validity());
  }
  GIVEN("A scalar diagram") {
    ZXDiagram diag3(0, 0, 0, 0);

    ZXVert v1 = diag3.add_vertex(ZXType::ZSpider);
    ZXVert v2 = diag3.add_vertex(ZXType::ZSpider);
    ZXVert v3 = diag3.add_vertex(ZXType::ZSpider, 3.22);
    ZXVert v4 = diag3.add_vertex(ZXType::ZSpider);
    ZXVert v5 = diag3.add_vertex(ZXType::ZSpider);
    ZXVert v6 = diag3.add_vertex(ZXType::ZSpider);

    Wire v1v4 = diag3.add_wire(v1, v4, ZXWireType::H);
    Wire v4v5 = diag3.add_wire(v4, v5, ZXWireType::Basic);
    Wire v5v4 = diag3.add_wire(v5, v4, ZXWireType::H);
    Wire v5v6 = diag3.add_wire(v5, v6, ZXWireType::Basic);
    Wire v6v3 = diag3.add_wire(v6, v3, ZXWireType::H);
    Wire v3v2 = diag3.add_wire(v3, v2, ZXWireType::Basic);
    Wire v2v3 = diag3.add_wire(v2, v3, ZXWireType::H);
    Wire v2v1 = diag3.add_wire(v2, v1, ZXWireType::Basic);

    REQUIRE_NOTHROW(diag3.check_validity());

    unsigned n_matches;

    // Self-loop finding
    CHECK_FALSE(Rewrite::self_loop_removal().apply(diag3));

    // Spider fusion
    CHECK(Rewrite::spider_fusion().apply(diag3));
    REQUIRE(diag3.count_vertices(ZXType::ZSpider) == 2);

    // Self loop removal after fusion
    CHECK(Rewrite::self_loop_removal().apply(diag3));

    // Parallel edge pair removal
    CHECK(Rewrite::parallel_h_removal().apply(diag3));

    // Remove hadamard edges connected directly to the boundaries
    CHECK_FALSE(Rewrite::io_extension().apply(diag3));

    // Deal with directly connected in/outputs
    CHECK_FALSE(Rewrite::separate_boundaries().apply(diag3));

    REQUIRE(diag3.count_vertices(ZXType::ZSpider) == 2);
    REQUIRE_NOTHROW(diag3.check_validity());
  }
}

SCENARIO("ZXBox decomposition") {
  GIVEN("Nested ZXBoxes") {
    ZXDiagram innermost(1, 0, 0, 2);
    ZXVertVec innermost_ins = innermost.get_boundary(ZXType::Input);
    ZXVertVec innermost_outs = innermost.get_boundary(ZXType::Output);
    ZXVert innermost_spid =
        innermost.add_vertex(ZXType::XSpider, QuantumType::Classical);
    innermost.add_wire(
        innermost_ins[0], innermost_spid, ZXWireType::Basic,
        QuantumType::Quantum);
    innermost.add_wire(
        innermost_outs[0], innermost_spid, ZXWireType::Basic,
        QuantumType::Classical);
    innermost.add_wire(
        innermost_outs[1], innermost_spid, ZXWireType::Basic,
        QuantumType::Classical);
    ZXGen_ptr inner_box_gen = std::make_shared<const ZXBox>(innermost);

    ZXDiagram inner(0, 2, 0, 0);
    ZXVertVec inner_outs = inner.get_boundary();
    ZXVert inner_box = inner.add_vertex(inner_box_gen);
    ZXVert inner_spid =
        inner.add_vertex(ZXType::ZSpider, QuantumType::Classical);
    inner.add_wire(
        inner_box, inner_outs[0], ZXWireType::H, QuantumType::Quantum, 0);
    inner.add_wire(
        inner_spid, inner_box, ZXWireType::Basic, QuantumType::Classical,
        std::nullopt, 1);
    inner.add_wire(
        inner_box, inner_spid, ZXWireType::Basic, QuantumType::Classical, 2);
    inner.add_wire(inner_spid, inner_outs[1]);
    ZXGen_ptr box_gen = std::make_shared<const ZXBox>(inner);

    ZXDiagram diag(1, 1, 0, 0);
    ZXVertVec b = diag.get_boundary();
    ZXVert box = diag.add_vertex(box_gen);
    ZXVert spid = diag.add_vertex(ZXType::ZSpider, 1., QuantumType::Quantum);
    diag.add_wire(box, b[0], ZXWireType::Basic, QuantumType::Quantum, 0);
    diag.add_wire(box, spid, ZXWireType::Basic, QuantumType::Quantum, 1);
    diag.add_wire(spid, b[1]);

    REQUIRE_NOTHROW(diag.check_validity());

    CHECK(Rewrite::decompose_boxes().apply(diag));

    CHECK(diag.count_vertices(ZXType::ZXBox) == 0);
    CHECK(diag.count_vertices(ZXType::ZSpider) == 2);
    CHECK(diag.count_vertices(ZXType::XSpider) == 1);

    CHECK(Rewrite::parallel_h_removal().apply(diag));
    CHECK(Rewrite::spider_fusion().apply(diag));

    REQUIRE_NOTHROW(diag.check_validity());
  }
}

SCENARIO("Mapping Hadamard edges to basic edges") {
  GIVEN("A diagram with a mixture of edge types") {
    ZXDiagram diag(1, 1, 1, 1);
    ZXVertVec ins = diag.get_boundary(ZXType::Input);
    ZXVertVec outs = diag.get_boundary(ZXType::Output);
    ZXVert z = diag.add_vertex(ZXType::ZSpider, QuantumType::Classical);
    ZXVert x = diag.add_vertex(ZXType::XSpider);
    diag.add_wire(ins[0], x, ZXWireType::H);
    diag.add_wire(ins[1], z, ZXWireType::H, QuantumType::Classical);
    diag.add_wire(outs[0], z, ZXWireType::Basic);
    diag.add_wire(outs[1], z, ZXWireType::Basic, QuantumType::Classical);

    REQUIRE_NOTHROW(diag.check_validity());

    CHECK(Rewrite::basic_wires().apply(diag));

    CHECK(diag.count_wires(ZXWireType::H) == 0);
    CHECK(diag.count_wires(ZXWireType::Basic) == 6);
    CHECK(diag.count_vertices(ZXType::Hbox) == 2);

    CHECK(diag.get_qtype(diag.neighbours(ins[0])[0]) == QuantumType::Quantum);
    CHECK(diag.get_qtype(diag.neighbours(ins[1])[0]) == QuantumType::Classical);

    CHECK_FALSE(Rewrite::basic_wires().apply(diag));
  }
}

SCENARIO("Check rewrite combinators act the same as individual rewrites") {
  ZXDiagram diag(1, 1, 0, 0);
  ZXVertVec ins = diag.get_boundary(ZXType::Input);
  ZXVertVec outs = diag.get_boundary(ZXType::Output);
  ZXVert z = diag.add_vertex(ZXType::ZSpider);
  ZXVert x = diag.add_vertex(ZXType::XSpider);
  diag.add_wire(ins[0], z);
  diag.add_wire(z, z);
  diag.add_wire(z, x, ZXWireType::H);
  diag.add_wire(z, x, ZXWireType::Basic);
  diag.add_wire(x, x, ZXWireType::H);
  diag.add_wire(x, outs[0]);
  GIVEN("A diagram undergoing a sequence of rewrites") {
    ZXDiagram copy = diag;
    Rewrite seq = Rewrite::sequence(
        {Rewrite::self_loop_removal(), Rewrite::spider_fusion()});
    CHECK(seq.apply(copy));        // Both should make changes
    CHECK(seq.apply(copy));        // More self loops created by fusion
    CHECK_FALSE(seq.apply(copy));  // No more changes
  }
  GIVEN("Iterated rewrites on a diagram") {
    ZXDiagram copy = diag;
    Rewrite seq = Rewrite::sequence(
        {Rewrite::self_loop_removal(), Rewrite::spider_fusion()});
    Rewrite loop = Rewrite::repeat(seq);
    CHECK(loop.apply(copy));        // Should iterate until completion
    CHECK_FALSE(loop.apply(copy));  // Check for completion
  }
  GIVEN("Iterated (with metric) rewrites on a diagram") {
    ZXDiagram copy = diag;
    Rewrite seq = Rewrite::sequence(
        {Rewrite::self_loop_removal(), Rewrite::spider_fusion()});
    Rewrite loop = Rewrite::repeat_with_metric(seq, &ZXDiagram::n_vertices);
    CHECK(loop.apply(copy));        // Should iterate until completion
    CHECK_FALSE(loop.apply(copy));  // Check for completion
  }
  GIVEN("Iterated (while) rewrites on a diagram") {
    ZXDiagram copy = diag;
    Rewrite loop = Rewrite::repeat_while(
        Rewrite::self_loop_removal(), Rewrite::spider_fusion());
    CHECK(loop.apply(copy));        // Should iterate until completion
    CHECK_FALSE(loop.apply(copy));  // Check for completion
  }
}

}  // namespace test_ZXAxioms
}  // namespace zx
}  // namespace tket
