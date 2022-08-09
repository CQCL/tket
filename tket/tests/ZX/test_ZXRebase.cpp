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
#include "ZX/ZXDiagram.hpp"

namespace tket::zx::test_ZXRebase {

SCENARIO("Decompose box with an identity wire") {
  GIVEN("Identity Box") {
    ZXDiagram iboxd(1, 1, 0, 0);
    ZXVertVec ib = iboxd.get_boundary();
    iboxd.add_wire(ib[0], ib[1]);
    REQUIRE_NOTHROW(iboxd.check_validity());
    ZXGen_ptr ibox_gen = std::make_shared<const ZXBox>(iboxd);
    ZXDiagram diag(1, 1, 0, 0);
    ZXVert box = diag.add_vertex(ibox_gen);
    ZXVertVec b = diag.get_boundary();
    diag.add_wire(box, b[0], ZXWireType::Basic, QuantumType::Quantum, 0);
    diag.add_wire(box, b[1], ZXWireType::Basic, QuantumType::Quantum, 1);

    REQUIRE_NOTHROW(diag.check_validity());

    CHECK(Rewrite::decompose_boxes().apply(diag));
    REQUIRE_NOTHROW(diag.check_validity());
    REQUIRE(diag.n_wires() == 1);
    REQUIRE(diag.n_vertices() == 2);
    ZXVertVec d = diag.get_boundary();
    REQUIRE(diag.neighbours(d[0])[0] == d[1]);
  }
}

SCENARIO("Take a generic diagram and apply each rebase to it") {
  ZXDiagram diag(2, 1, 0, 1);
  ZXVertVec ins = diag.get_boundary(ZXType::Input);
  ZXVertVec outs = diag.get_boundary(ZXType::Output);
  ZXVert h0 = diag.add_vertex(ZXType::Hbox);
  ZXVert h1 = diag.add_vertex(ZXType::Hbox, -3.7, QuantumType::Classical);
  ZXVert xy = diag.add_vertex(ZXType::XY, 0.4);
  ZXVert xz = diag.add_vertex(ZXType::XZ, 0.7, QuantumType::Classical);
  ZXVert yz = diag.add_vertex(ZXType::YZ, 1.2);
  ZXVert px = diag.add_vertex(
      ZXGen::create_gen(ZXType::PX, false, QuantumType::Classical));
  ZXVert py = diag.add_vertex(
      ZXGen::create_gen(ZXType::PY, true, QuantumType::Classical));
  ZXVert pz = diag.add_vertex(ZXGen::create_gen(ZXType::PZ, true));
  ZXVert zspid = diag.add_vertex(ZXType::ZSpider, 0.9);
  ZXVert xspid = diag.add_vertex(ZXType::XSpider, 1.8, QuantumType::Classical);
  diag.add_wire(ins.at(0), h0);
  diag.add_wire(h0, xy);
  diag.add_wire(xy, yz);
  diag.add_wire(xy, pz);
  diag.add_wire(pz, outs.at(0));
  diag.add_wire(yz, zspid);
  diag.add_wire(ins.at(1), xspid);
  diag.add_wire(xspid, xz, ZXWireType::Basic, QuantumType::Classical);
  diag.add_wire(xz, px, ZXWireType::Basic, QuantumType::Classical);
  diag.add_wire(xz, py, ZXWireType::Basic, QuantumType::Classical);
  diag.add_wire(py, h1, ZXWireType::Basic, QuantumType::Classical);
  diag.add_wire(h1, h1);
  diag.add_wire(h1, outs.at(1), ZXWireType::Basic, QuantumType::Classical);
  REQUIRE_NOTHROW(diag.check_validity());

  // Just check for gate counts here; check for semantic preservation in python
  GIVEN("Rebasing to ZX") {
    Rewrite::rebase_to_zx().apply(diag);
    REQUIRE_NOTHROW(diag.check_validity());
    CHECK(diag.count_vertices(ZXType::Hbox) == 0);
    CHECK(diag.count_vertices(ZXType::XY) == 0);
    CHECK(diag.count_vertices(ZXType::XZ) == 0);
    CHECK(diag.count_vertices(ZXType::YZ) == 0);
    CHECK(diag.count_vertices(ZXType::PX) == 0);
    CHECK(diag.count_vertices(ZXType::PY) == 0);
    CHECK(diag.count_vertices(ZXType::PZ) == 0);
    CHECK(diag.count_vertices(ZXType::Triangle) == 0);
    CHECK(diag.count_vertices(ZXType::ZXBox) == 0);
  }
  GIVEN("Rebasing to MBQC") {
    Rewrite::rebase_to_mbqc().apply(diag);
    REQUIRE_NOTHROW(diag.check_validity());
    CHECK(diag.count_vertices(ZXType::Hbox) == 0);
    CHECK(diag.count_vertices(ZXType::ZSpider) == 0);
    CHECK(diag.count_vertices(ZXType::XSpider) == 0);
    CHECK(diag.count_vertices(ZXType::Triangle) == 0);
    CHECK(diag.count_vertices(ZXType::ZXBox) == 0);
  }
}

}  // namespace tket::zx::test_ZXRebase
