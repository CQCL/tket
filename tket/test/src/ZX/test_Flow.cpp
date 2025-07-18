// Copyright Quantinuum
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
#include <catch2/matchers/catch_matchers.hpp>

#include "tket/ZX/Flow.hpp"

namespace tket {

namespace zx {

namespace test_flow {

SCENARIO("Testing flow verification") {
  // Diagram combines Ex. 2.43, "There and back again: a circuit extraction
  // tale", Backens et al. 2021 and Ex. C.13, "Relating measurement patterns to
  // circuits via Pauli flow", Simmons 2021
  ZXDiagram diag(1, 3, 0, 0);
  ZXVertVec ins = diag.get_boundary(ZXType::Input);
  ZXVertVec outs = diag.get_boundary(ZXType::Output);
  // Gflow example from Backens et al.
  ZXVert ga = diag.add_vertex(ZXType::XY, 0.3);
  ZXVert gb = diag.add_vertex(ZXType::XY, 0.7);
  ZXVert gc = diag.add_vertex(ZXType::XZ, 1.4);
  ZXVert gd = diag.add_vertex(ZXType::YZ, 0.9);
  ZXVert o0 = diag.add_vertex(ZXType::PX);
  ZXVert o1 = diag.add_vertex(ZXType::PX);
  ZXVert o2 = diag.add_vertex(ZXType::PX);
  diag.add_wire(ins.at(0), ga);
  diag.add_wire(ga, gb, ZXWireType::H);
  diag.add_wire(gb, gc, ZXWireType::H);
  diag.add_wire(gb, gd, ZXWireType::H);
  diag.add_wire(gc, gd, ZXWireType::H);
  diag.add_wire(gb, o0, ZXWireType::H);
  diag.add_wire(o0, outs.at(0));
  // Pauli flow example from Simmons (angles cut to Paulis)
  ZXVert pi = diag.add_vertex(ZXType::XY, 0.9);
  ZXVert pa = diag.add_vertex(ZXType::PZ);
  ZXVert pb = diag.add_vertex(ZXType::PX);
  ZXVert pc = diag.add_vertex(ZXType::XY, 0.2);
  ZXVert pd = diag.add_vertex(ZXGen::create_gen(ZXType::PY, true));
  diag.add_wire(gc, pi, ZXWireType::H);
  diag.add_wire(pi, pb, ZXWireType::H);
  diag.add_wire(pa, pb, ZXWireType::H);
  diag.add_wire(pa, pc, ZXWireType::H);
  diag.add_wire(pa, pd, ZXWireType::H);
  diag.add_wire(pb, pd, ZXWireType::H);
  diag.add_wire(pc, pd, ZXWireType::H);
  diag.add_wire(pc, o1, ZXWireType::H);
  diag.add_wire(pd, o2, ZXWireType::H);
  diag.add_wire(o1, outs.at(1));
  diag.add_wire(o2, outs.at(2));

  // Give a valid Pauli flow
  std::map<ZXVert, ZXVertSeqSet> c{
      {ga, {gb}},              // Odd = {ga, gc, gd, o0}
      {gb, {gc}},              // Odd = {gb, gc, pi}
      {gc, {gc, gd}},          // Odd = {gc, gd, pi}
      {gd, {gd, o0, pi}},      // Odd = {pb}
      {pi, {pb, o2}},          // Odd = {pi, pa}
      {pa, {pa, pc, pd, o2}},  // Odd = {pd, o1, o2}
      {pb, {pc, pd, o1}},      // Odd = {pb, pd, o1, o2}
      {pc, {o1}},              // Odd = {pc}
      {pd, {o2}},              // Odd = {pd}
  };
  std::map<ZXVert, unsigned> d{
      {ga, 7}, {gb, 6}, {gc, 5}, {gd, 4}, {pi, 3}, {pa, 2},
      {pb, 2}, {pc, 1}, {pd, 1}, {o0, 0}, {o1, 0}, {o2, 0},
  };

  Flow fl{c, d};
  REQUIRE_NOTHROW(fl.verify(diag));

  // Check for ordering of corrections
  d.at(ga) = 4;
  fl = {c, d};
  REQUIRE_THROWS_WITH(
      fl.verify(diag), "A qubit has an X correction in its past");
  d.at(gb) = 3;
  fl = {c, d};
  REQUIRE_THROWS_WITH(
      fl.verify(diag), "A qubit has a Z correction in its past");
  // Revert to valid flow
  d.at(ga) = 7;
  d.at(gb) = 6;

  // Check history Y measurements have Y corrections
  diag.set_vertex_ZXGen_ptr(pb, ZXGen::create_gen(ZXType::PY));
  c.at(pa) = {pa};
  fl = {c, d};
  REQUIRE_THROWS_WITH(
      fl.verify(diag), "A past Y vertex receives a Z correction");
  c.at(pa) = {pa, pc, pd};
  d.at(pd) = 2;
  fl = {c, d};
  REQUIRE_THROWS_WITH(
      fl.verify(diag), "A past Y vertex receives an X correction");
  // Revert to valid flow
  diag.set_vertex_ZXGen_ptr(pb, ZXGen::create_gen(ZXType::PX));
  c.at(pa) = {pa, pc, pd, o2};
  d.at(pd) = 1;

  // Check all basis corrections are ok
  // Correct XY with I, X, Y
  std::vector<ZXVertSeqSet> cs{{}, {pc, o2}, {pc, o1, o2}};
  for (const ZXVertSeqSet& cc : cs) {
    c.at(pc) = cc;
    fl = {c, d};
    REQUIRE_THROWS_WITH(
        fl.verify(diag), "XY vertex must be corrected with a Z");
  }
  c.at(pc) = {o1};
  // Correct XZ with I, X, Z
  cs = {{}, {gc, o0}, {pi}};
  for (const ZXVertSeqSet& cc : cs) {
    c.at(gc) = cc;
    fl = {c, d};
    REQUIRE_THROWS_WITH(
        fl.verify(diag), "XZ vertex must be corrected with a Y");
  }
  c.at(gc) = {gc, gd};
  // Correct YZ with I, Y, Z
  diag.set_vertex_ZXGen_ptr(pa, ZXGen::create_gen(ZXType::YZ, Expr(1.2)));
  cs = {{}, {pa, pd}, {pc}};
  for (const ZXVertSeqSet& cc : cs) {
    c.at(pa) = cc;
    fl = {c, d};
    REQUIRE_THROWS_WITH(
        fl.verify(diag), "YZ vertex must be corrected with an X");
  }
  diag.set_vertex_ZXGen_ptr(pa, ZXGen::create_gen(ZXType::PZ));
  c.at(pa) = {pa, pc, pd, o2};
  // Correct PX with I, X
  diag.set_vertex_ZXGen_ptr(pc, ZXGen::create_gen(ZXType::PX));
  cs = {{}, {pc, o2}};
  for (const ZXVertSeqSet& cc : cs) {
    c.at(pc) = cc;
    fl = {c, d};
    REQUIRE_THROWS_WITH(
        fl.verify(diag), "PX vertex must be corrected with a Y or Z");
  }
  diag.set_vertex_ZXGen_ptr(pc, ZXGen::create_gen(ZXType::XY, Expr(0.2)));
  c.at(pc) = {o1};
  // Correct PY with I, Y
  diag.set_vertex_ZXGen_ptr(pc, ZXGen::create_gen(ZXType::PY));
  cs = {{}, {pc, o1, o2}};
  for (const ZXVertSeqSet& cc : cs) {
    c.at(pc) = cc;
    fl = {c, d};
    REQUIRE_THROWS_WITH(
        fl.verify(diag), "PY vertex must be corrected with an X or Z");
  }
  diag.set_vertex_ZXGen_ptr(pc, ZXGen::create_gen(ZXType::XY, Expr(0.2)));
  c.at(pc) = {o1};
  // Correct PZ with I, Z
  cs = {{}, {pc, o2}};
  for (const ZXVertSeqSet& cc : cs) {
    c.at(pa) = cc;
    fl = {c, d};
    REQUIRE_THROWS_WITH(
        fl.verify(diag), "PZ vertex must be corrected with an X or Y");
  }
}

SCENARIO("Testing causal flow identification and focussing") {
  // Diagram based on Fig. 8, "Determinism in the one-way model",
  // Danos & Kashefi 2006
  ZXDiagram diag(2, 2, 0, 0);
  ZXVertVec ins = diag.get_boundary(ZXType::Input);
  ZXVertVec outs = diag.get_boundary(ZXType::Output);
  // Input measurements
  ZXVert i0 = diag.add_vertex(ZXType::XY, 0.3);
  ZXVert i1 = diag.add_vertex(ZXType::XY, 0.7);
  diag.add_wire(ins.at(0), i0);
  diag.add_wire(ins.at(1), i1);
  // Chain on qubit 0
  ZXVert v0 = diag.add_vertex(ZXType::XY, 1.4);
  ZXVert o0 = diag.add_vertex(ZXType::PX);
  diag.add_wire(i0, v0, ZXWireType::H);
  diag.add_wire(v0, o0, ZXWireType::H);
  diag.add_wire(o0, outs.at(0));
  // Chain on qubit 1
  ZXVert v1a = diag.add_vertex(ZXType::XY, 0.9);
  ZXVert v1b = diag.add_vertex(ZXType::XY, 0.2);
  ZXVert v1c = diag.add_vertex(ZXType::XY, 1.2);
  ZXVert v1d = diag.add_vertex(ZXType::XY, 1.6);
  ZXVert v1e = diag.add_vertex(ZXType::XY, 0.4);
  ZXVert o1 = diag.add_vertex(ZXType::PX);
  diag.add_wire(i1, v1a, ZXWireType::H);
  diag.add_wire(v1a, v1b, ZXWireType::H);
  diag.add_wire(v1b, v1c, ZXWireType::H);
  diag.add_wire(v1c, v1d, ZXWireType::H);
  diag.add_wire(v1d, v1e, ZXWireType::H);
  diag.add_wire(v1e, o1, ZXWireType::H);
  diag.add_wire(o1, outs.at(1));
  // Cross-chain links
  diag.add_wire(i0, v1a, ZXWireType::H);
  diag.add_wire(i0, v1d, ZXWireType::H);

  Flow f = Flow::identify_causal_flow(diag);

  CHECK(f.c(i0) == ZXVertSeqSet{v0});
  CHECK(f.c(v0) == ZXVertSeqSet{o0});
  CHECK(f.c(i1) == ZXVertSeqSet{v1a});
  CHECK(f.c(v1a) == ZXVertSeqSet{v1b});
  CHECK(f.c(v1b) == ZXVertSeqSet{v1c});
  CHECK(f.c(v1c) == ZXVertSeqSet{v1d});
  CHECK(f.c(v1d) == ZXVertSeqSet{v1e});
  CHECK(f.c(v1e) == ZXVertSeqSet{o1});
  REQUIRE_NOTHROW(f.verify(diag));

  REQUIRE_NOTHROW(f.focus(diag));
  CHECK(f.c(i0) == ZXVertSeqSet{v0});
  CHECK(f.c(v0) == ZXVertSeqSet{o0});
  CHECK(f.c(i1) == ZXVertSeqSet{v1a, v0, v1c, v1e});
  CHECK(f.c(v1a) == ZXVertSeqSet{v1b, v1d, v0, o1});
  CHECK(f.c(v1b) == ZXVertSeqSet{v1c, v1e});
  CHECK(f.c(v1c) == ZXVertSeqSet{v1d, v0, o1});
  CHECK(f.c(v1d) == ZXVertSeqSet{v1e});
  CHECK(f.c(v1e) == ZXVertSeqSet{o1});
  REQUIRE_NOTHROW(f.verify(diag));
}

SCENARIO("Testing Pauli flow identification and focussing") {
  // Diagram combines Ex. 2.43, "There and back again: a circuit extraction
  // tale", Backens et al. 2021 and Ex. C.13, "Relating measurement patterns to
  // circuits via Pauli flow", Simmons 2021
  ZXDiagram diag(1, 3, 0, 0);
  ZXVertVec ins = diag.get_boundary(ZXType::Input);
  ZXVertVec outs = diag.get_boundary(ZXType::Output);
  // Gflow example from Backens et al.
  ZXVert ga = diag.add_vertex(ZXType::XY, 0.3);
  ZXVert gb = diag.add_vertex(ZXType::XY, 0.7);
  ZXVert gc = diag.add_vertex(ZXType::XZ, 1.4);
  ZXVert gd = diag.add_vertex(ZXType::YZ, 0.9);
  ZXVert o0 = diag.add_vertex(ZXType::PX);
  diag.add_wire(ins.at(0), ga);
  diag.add_wire(ga, gb, ZXWireType::H);
  diag.add_wire(gb, gc, ZXWireType::H);
  diag.add_wire(gb, gd, ZXWireType::H);
  diag.add_wire(gc, gd, ZXWireType::H);
  diag.add_wire(gb, o0, ZXWireType::H);
  diag.add_wire(o0, outs.at(0));
  // Pauli flow example from Simmons (angles cut to Paulis)
  ZXVert pi = diag.add_vertex(ZXType::XY, 0.9);
  ZXVert pa = diag.add_vertex(ZXType::PZ);
  ZXVert pb = diag.add_vertex(ZXType::PX);
  ZXVert pc = diag.add_vertex(ZXType::XY, 0.2);
  ZXVert pd = diag.add_vertex(ZXGen::create_gen(ZXType::PY, true));
  ZXVert o1 = diag.add_vertex(ZXType::PX);
  ZXVert o2 = diag.add_vertex(ZXType::PX);
  diag.add_wire(gc, pi, ZXWireType::H);
  diag.add_wire(pi, pb, ZXWireType::H);
  diag.add_wire(pa, pb, ZXWireType::H);
  diag.add_wire(pa, pc, ZXWireType::H);
  diag.add_wire(pa, pd, ZXWireType::H);
  diag.add_wire(pb, pd, ZXWireType::H);
  diag.add_wire(pc, pd, ZXWireType::H);
  diag.add_wire(pc, o1, ZXWireType::H);
  diag.add_wire(pd, o2, ZXWireType::H);
  diag.add_wire(o1, outs.at(1));
  diag.add_wire(o2, outs.at(2));

  Flow f = Flow::identify_pauli_flow(diag);

  // Check exact sets for deterministic output
  CHECK(f.d(ga) == 3);
  CHECK(f.d(gb) == 1);
  CHECK(f.d(gc) == 3);
  CHECK(f.d(gd) == 2);
  CHECK(f.d(o0) == 0);
  CHECK(f.d(pi) == 1);
  CHECK(f.d(pa) == 1);
  CHECK(f.d(pb) == 1);
  CHECK(f.d(pc) == 1);
  CHECK(f.d(pd) == 1);
  CHECK(f.d(o1) == 0);
  CHECK(f.d(o2) == 0);
  CHECK(f.c(ga).get<TagSeq>() == ZXVertSeqSet{gb, pi}.get<TagSeq>());
  CHECK(f.c(gb).get<TagSeq>() == ZXVertSeqSet{o0}.get<TagSeq>());
  CHECK(f.c(gc).get<TagSeq>() == ZXVertSeqSet{pi, gc}.get<TagSeq>());
  CHECK(f.c(gd).get<TagSeq>() == ZXVertSeqSet{pi, gd}.get<TagSeq>());
  CHECK(f.c(o0).get<TagSeq>() == ZXVertSeqSet{}.get<TagSeq>());
  CHECK(f.c(pi).get<TagSeq>() == ZXVertSeqSet{pb, o2}.get<TagSeq>());
  CHECK(f.c(pa).get<TagSeq>() == ZXVertSeqSet{pd, pa}.get<TagSeq>());
  CHECK(f.c(pb).get<TagSeq>() == ZXVertSeqSet{pd, o1, o2}.get<TagSeq>());
  CHECK(f.c(pc).get<TagSeq>() == ZXVertSeqSet{o1}.get<TagSeq>());
  CHECK(f.c(pd).get<TagSeq>() == ZXVertSeqSet{o2}.get<TagSeq>());
  CHECK(f.c(o1).get<TagSeq>() == ZXVertSeqSet{}.get<TagSeq>());
  CHECK(f.c(o2).get<TagSeq>() == ZXVertSeqSet{}.get<TagSeq>());

  REQUIRE_NOTHROW(f.verify(diag));
  REQUIRE_NOTHROW(f.focus(diag));
  REQUIRE_NOTHROW(f.verify(diag));
}

SCENARIO("Test focussed set identification") {
  // Diagram combines Ex. 2.43, "There and back again: a circuit extraction
  // tale", Backens et al. 2021 and Ex. C.13, "Relating measurement patterns to
  // circuits via Pauli flow", Simmons 2021
  ZXDiagram diag(1, 3, 0, 0);
  ZXVertVec ins = diag.get_boundary(ZXType::Input);
  ZXVertVec outs = diag.get_boundary(ZXType::Output);
  // Gflow example from Backens et al.
  ZXVert ga = diag.add_vertex(ZXType::XY, 0.3);
  ZXVert gb = diag.add_vertex(ZXType::XY, 0.7);
  ZXVert gc = diag.add_vertex(ZXType::XZ, 1.4);
  ZXVert gd = diag.add_vertex(ZXType::YZ, 0.9);
  ZXVert o0 = diag.add_vertex(ZXType::PX);
  diag.add_wire(ins.at(0), ga);
  diag.add_wire(ga, gb, ZXWireType::H);
  diag.add_wire(gb, gc, ZXWireType::H);
  diag.add_wire(gb, gd, ZXWireType::H);
  diag.add_wire(gc, gd, ZXWireType::H);
  diag.add_wire(gb, o0, ZXWireType::H);
  diag.add_wire(o0, outs.at(0));
  // Pauli flow example from Simmons (angles cut to Paulis)
  ZXVert pi = diag.add_vertex(ZXType::XY, 0.9);
  ZXVert pa = diag.add_vertex(ZXType::PZ);
  ZXVert pb = diag.add_vertex(ZXType::PX);
  ZXVert pc = diag.add_vertex(ZXType::XY, 0.2);
  ZXVert pd = diag.add_vertex(ZXGen::create_gen(ZXType::PY, true));
  ZXVert o1 = diag.add_vertex(ZXType::PX);
  ZXVert o2 = diag.add_vertex(ZXType::PX);
  diag.add_wire(gc, pi, ZXWireType::H);
  diag.add_wire(pi, pb, ZXWireType::H);
  diag.add_wire(pa, pb, ZXWireType::H);
  diag.add_wire(pa, pc, ZXWireType::H);
  diag.add_wire(pa, pd, ZXWireType::H);
  diag.add_wire(pb, pd, ZXWireType::H);
  diag.add_wire(pc, pd, ZXWireType::H);
  diag.add_wire(pc, o1, ZXWireType::H);
  diag.add_wire(pd, o2, ZXWireType::H);
  diag.add_wire(o1, outs.at(1));
  diag.add_wire(o2, outs.at(2));
  std::set<ZXVert> output_set = {o0, o1, o2};

  std::set<ZXVertSeqSet> focussed = Flow::identify_focussed_sets(diag);

  REQUIRE(focussed.size() == 2);
  for (const ZXVertSeqSet& fset : focussed) {
    std::map<ZXVert, unsigned> parities;
    for (const ZXVert& v : fset.get<TagSeq>()) {
      ZXType vtype = diag.get_zxtype(v);
      REQUIRE(
          (vtype == ZXType::XY || vtype == ZXType::PX || vtype == ZXType::PY));
      for (const ZXVert& n : fset.get<TagSeq>()) {
        auto inserted = parities.insert({n, 1});
        if (!inserted.second) {
          ++(inserted.first->second);
        }
      }
    }
    for (const std::pair<const ZXVert, unsigned>& p : parities) {
      if (p.second % 2 == 1) {
        ZXType vtype = diag.get_zxtype(p.first);
        REQUIRE(
            (vtype == ZXType::XZ || vtype == ZXType::YZ ||
             vtype == ZXType::PY || vtype == ZXType::PZ ||
             output_set.find(p.first) != output_set.end()));
        REQUIRE(
            (vtype != ZXType::PY ||
             fset.get<TagKey>().find(p.first) != fset.get<TagKey>().end()));
      }
    }
  }
}

SCENARIO("Test example sent by Korbinian Staudacher") {
  ZXDiagram diag(4, 4, 0, 0);
  ZXVertVec ins = diag.get_boundary(ZXType::Input);
  ZXVertVec outs = diag.get_boundary(ZXType::Output);
  ZXVert v0 = diag.add_vertex(ZXType::XY, 1.25);
  ZXVert v1 = diag.add_clifford_vertex(ZXType::PX, false);
  ZXVert v2 = diag.add_vertex(ZXType::XY, 1.25);
  ZXVert v3 = diag.add_vertex(ZXType::XY, 1.75);
  ZXVert v4 = diag.add_clifford_vertex(ZXType::PY, false);
  ZXVert v5 = diag.add_clifford_vertex(ZXType::PY, true);
  ZXVert v6 = diag.add_clifford_vertex(ZXType::PX, true);
  ZXVert v7 = diag.add_clifford_vertex(ZXType::PY, true);
  ZXVert v8 = diag.add_vertex(ZXType::XY, 0.75);
  ZXVert v9 = diag.add_clifford_vertex(ZXType::PX, false);
  ZXVert v10 = diag.add_vertex(ZXType::XY, 0.75);
  ZXVert v11 = diag.add_clifford_vertex(ZXType::PY, true);
  ZXVert v12 = diag.add_clifford_vertex(ZXType::PX, false);
  ZXVert v13 = diag.add_clifford_vertex(ZXType::PX, false);
  ZXVert v14 = diag.add_clifford_vertex(ZXType::PX, false);
  ZXVert v15 = diag.add_clifford_vertex(ZXType::PX, false);
  ZXVert v16 = diag.add_clifford_vertex(ZXType::PX, false);

  std::map<ZXVert, unsigned> vmap{
      {v0, 0},   {v1, 1},   {v2, 2},   {v3, 3},   {v4, 4},   {v5, 5},
      {v6, 6},   {v7, 7},   {v8, 8},   {v9, 9},   {v10, 10}, {v11, 11},
      {v12, 12}, {v13, 13}, {v14, 14}, {v15, 15}, {v16, 16},
  };

  diag.add_wire(ins[0], v1);
  diag.add_wire(ins[1], v2);
  diag.add_wire(ins[2], v3);
  diag.add_wire(ins[3], v0);
  diag.add_wire(v13, outs[0]);
  diag.add_wire(v16, outs[1]);
  diag.add_wire(v15, outs[2]);
  diag.add_wire(v14, outs[3]);

  diag.add_wire(v0, v5, ZXWireType::H);
  diag.add_wire(v1, v2, ZXWireType::H);
  diag.add_wire(v1, v6, ZXWireType::H);
  diag.add_wire(v2, v4, ZXWireType::H);
  diag.add_wire(v3, v6, ZXWireType::H);
  diag.add_wire(v3, v9, ZXWireType::H);
  diag.add_wire(v3, v12, ZXWireType::H);
  diag.add_wire(v4, v7, ZXWireType::H);
  diag.add_wire(v5, v6, ZXWireType::H);
  diag.add_wire(v5, v8, ZXWireType::H);
  diag.add_wire(v6, v12, ZXWireType::H);
  diag.add_wire(v7, v9, ZXWireType::H);
  diag.add_wire(v7, v10, ZXWireType::H);
  diag.add_wire(v7, v12, ZXWireType::H);
  diag.add_wire(v7, v16, ZXWireType::H);
  diag.add_wire(v8, v11, ZXWireType::H);
  diag.add_wire(v9, v10, ZXWireType::H);
  diag.add_wire(v10, v15, ZXWireType::H);
  diag.add_wire(v11, v14, ZXWireType::H);
  diag.add_wire(v12, v13, ZXWireType::H);

  Flow f = Flow::identify_pauli_flow(diag);

  // Check exact sets for deterministic output
  CHECK(f.d(v0) == 2);
  CHECK(f.d(v1) == 2);
  CHECK(f.d(v2) == 2);
  CHECK(f.d(v3) == 1);
  CHECK(f.d(v4) == 2);
  CHECK(f.d(v5) == 2);
  CHECK(f.d(v6) == 1);
  CHECK(f.d(v7) == 1);
  CHECK(f.d(v8) == 1);
  CHECK(f.d(v9) == 2);
  CHECK(f.d(v10) == 1);
  CHECK(f.d(v11) == 1);
  CHECK(f.d(v12) == 1);
  CHECK(f.d(v13) == 0);
  CHECK(f.d(v14) == 0);
  CHECK(f.d(v15) == 0);
  CHECK(f.d(v16) == 0);
  CHECK(f.c(v0).get<TagSeq>() == ZXVertSeqSet{v5, v8}.get<TagSeq>());
  CHECK(f.c(v1).get<TagSeq>() == ZXVertSeqSet{v6, v8}.get<TagSeq>());
  CHECK(f.c(v2).get<TagSeq>() == ZXVertSeqSet{v4, v7, v10}.get<TagSeq>());
  CHECK(f.c(v3).get<TagSeq>() == ZXVertSeqSet{v9, v15, v16}.get<TagSeq>());
  CHECK(f.c(v4).get<TagSeq>() == ZXVertSeqSet{v7, v10}.get<TagSeq>());
  CHECK(f.c(v5).get<TagSeq>() == ZXVertSeqSet{v8}.get<TagSeq>());
  CHECK(f.c(v6).get<TagSeq>() == ZXVertSeqSet{v9, v12, v15}.get<TagSeq>());
  CHECK(f.c(v7).get<TagSeq>() == ZXVertSeqSet{v16}.get<TagSeq>());
  CHECK(f.c(v8).get<TagSeq>() == ZXVertSeqSet{v11, v14}.get<TagSeq>());
  CHECK(f.c(v9).get<TagSeq>() == ZXVertSeqSet{v10}.get<TagSeq>());
  CHECK(f.c(v10).get<TagSeq>() == ZXVertSeqSet{v15}.get<TagSeq>());
  CHECK(f.c(v11).get<TagSeq>() == ZXVertSeqSet{v14}.get<TagSeq>());
  CHECK(f.c(v12).get<TagSeq>() == ZXVertSeqSet{v13}.get<TagSeq>());
  CHECK(f.c(v13).get<TagSeq>() == ZXVertSeqSet{}.get<TagSeq>());
  CHECK(f.c(v14).get<TagSeq>() == ZXVertSeqSet{}.get<TagSeq>());
  CHECK(f.c(v15).get<TagSeq>() == ZXVertSeqSet{}.get<TagSeq>());
  CHECK(f.c(v16).get<TagSeq>() == ZXVertSeqSet{}.get<TagSeq>());

  REQUIRE_NOTHROW(f.verify(diag));
  REQUIRE_NOTHROW(f.focus(diag));
  REQUIRE_NOTHROW(f.verify(diag));
}

}  // namespace test_flow

}  // namespace zx

}  // namespace tket
