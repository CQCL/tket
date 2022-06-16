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

#include "Converters/Converters.hpp"
#include "ZX/Flow.hpp"

namespace tket::zx::test_ZXExtraction {

// Just check for completion and validity here
// Check for semantic preservation in python

SCENARIO("Extracting a circuit with XY-gflow") {
  // Based on taking multiple instances of Fig. 2, "Generalised flow and
  // determinism in measurement-based quantum computation", Dan Browne et al.
  // 2007
  ZXDiagram diag(3, 3, 0, 0);
  ZXVertVec ins = diag.get_boundary(ZXType::Input);
  ZXVertVec outs = diag.get_boundary(ZXType::Output);
  ZXVert v00 = diag.add_vertex(ZXType::XY, 0.7);
  ZXVert v01 = diag.add_vertex(ZXType::XY, 0.2);
  ZXVert v02 = diag.add_vertex(ZXType::XY, 1.9);
  ZXVert v10 = diag.add_vertex(ZXType::XY, 0.56);
  ZXVert v11 = diag.add_vertex(ZXType::XY, 1.2);
  ZXVert v12 = diag.add_vertex(ZXType::XY, 0.9);
  ZXVert o0 = diag.add_vertex(ZXType::PX);
  ZXVert o1 = diag.add_vertex(ZXType::PX);
  ZXVert o2 = diag.add_vertex(ZXType::PX);
  diag.add_wire(ins.at(0), v00);
  diag.add_wire(ins.at(1), v01);
  diag.add_wire(ins.at(2), v02);
  diag.add_wire(v00, v10, ZXWireType::H);
  diag.add_wire(v00, v12, ZXWireType::H);
  diag.add_wire(v01, v10, ZXWireType::H);
  diag.add_wire(v01, v11, ZXWireType::H);
  diag.add_wire(v01, v12, ZXWireType::H);
  diag.add_wire(v02, v11, ZXWireType::H);
  diag.add_wire(v02, v12, ZXWireType::H);
  diag.add_wire(v10, o0, ZXWireType::H);
  diag.add_wire(v10, o2, ZXWireType::H);
  diag.add_wire(v11, o0, ZXWireType::H);
  diag.add_wire(v11, o1, ZXWireType::H);
  diag.add_wire(v11, o2, ZXWireType::H);
  diag.add_wire(v12, o1, ZXWireType::H);
  diag.add_wire(v12, o2, ZXWireType::H);
  diag.add_wire(o0, outs.at(0));
  diag.add_wire(o1, outs.at(1));
  diag.add_wire(o2, outs.at(2));
  REQUIRE_NOTHROW(Flow::identify_pauli_flow(diag));
  Circuit c = zx_to_circuit(diag);
  REQUIRE_NOTHROW(c.assert_valid());
  CHECK(c.n_qubits() == 3);
}

SCENARIO("Extracting a circuit from XY/YZ diagram with gflow") {
  // Diagram from Fig. 1(c), "Reducing the number of non-Clifford gates in
  // quantum circuits", Aleks Kissinger & John Van de Wetering, 2020
  ZXDiagram diag(5, 5, 0, 0);
  ZXVertVec ins = diag.get_boundary(ZXType::Input);
  ZXVertVec outs = diag.get_boundary(ZXType::Output);
  ZXVert i0 = diag.add_vertex(ZXType::XY);
  ZXVert i1 = diag.add_vertex(ZXType::XY);
  ZXVert i2 = diag.add_vertex(ZXType::XY, 0.25);
  ZXVert i3ext = diag.add_vertex(ZXType::XY);
  ZXVert i3 = diag.add_vertex(ZXType::XY, 0.25);
  ZXVert i4 = diag.add_vertex(ZXType::XY);
  ZXVert inter0 = diag.add_vertex(ZXType::XY, -0.25);
  ZXVert inter1 = diag.add_vertex(ZXType::XY, -0.25);
  ZXVert o0 = diag.add_vertex(ZXType::XY);
  ZXVert o0ext = diag.add_vertex(ZXType::PX);
  ZXVert o1 = diag.add_vertex(ZXType::XY);
  ZXVert o1ext = diag.add_vertex(ZXType::PX);
  ZXVert o2 = diag.add_vertex(ZXType::XY);
  ZXVert o2ext = diag.add_vertex(ZXType::PX);
  ZXVert o3 = diag.add_vertex(ZXType::PX);
  ZXVert o4 = diag.add_vertex(ZXType::XY, 0.25);
  ZXVert o4ext = diag.add_vertex(ZXType::PX);
  ZXVert g0 = diag.add_vertex(ZXType::YZ, -0.25);
  ZXVert g1 = diag.add_vertex(ZXType::YZ, 0.25);
  ZXVert g2 = diag.add_vertex(ZXType::YZ, 0.25);
  ZXVert g3 = diag.add_vertex(ZXType::YZ, 0.25);
  ZXVert g4 = diag.add_vertex(ZXType::YZ, 0.25);
  ZXVert g5 = diag.add_vertex(ZXType::YZ, 0.25);
  ZXVert g6 = diag.add_vertex(ZXType::YZ, -0.25);
  ZXVert g7 = diag.add_vertex(ZXType::YZ, -0.25);
  ZXVert g8 = diag.add_vertex(ZXType::YZ, -0.25);
  ZXVert g9 = diag.add_vertex(ZXType::YZ, -0.25);
  // Input wires
  diag.add_wire(ins.at(0), i0);
  diag.add_wire(ins.at(1), i1);
  diag.add_wire(ins.at(2), i2);
  diag.add_wire(ins.at(3), i3ext);
  diag.add_wire(ins.at(4), i4);
  diag.add_wire(i3ext, i3, ZXWireType::H);

  // Interior wires
  diag.add_wire(i0, i1, ZXWireType::H);
  diag.add_wire(i0, i3, ZXWireType::H);
  diag.add_wire(i0, i4, ZXWireType::H);
  diag.add_wire(i0, inter1, ZXWireType::H);
  diag.add_wire(i0, o0, ZXWireType::H);
  diag.add_wire(i1, o1, ZXWireType::H);
  diag.add_wire(i2, o2, ZXWireType::H);
  diag.add_wire(i3, inter0, ZXWireType::H);
  diag.add_wire(i3, o3, ZXWireType::H);
  diag.add_wire(i3, o4, ZXWireType::H);
  diag.add_wire(i4, inter0, ZXWireType::H);
  diag.add_wire(inter0, inter1, ZXWireType::H);
  diag.add_wire(inter1, o4, ZXWireType::H);

  // Gadget wires
  diag.add_wire(g0, i0, ZXWireType::H);
  diag.add_wire(g0, i1, ZXWireType::H);
  diag.add_wire(g0, inter0, ZXWireType::H);
  diag.add_wire(g1, i1, ZXWireType::H);
  diag.add_wire(g1, inter0, ZXWireType::H);
  diag.add_wire(g2, i0, ZXWireType::H);
  diag.add_wire(g2, inter0, ZXWireType::H);
  diag.add_wire(g3, i0, ZXWireType::H);
  diag.add_wire(g3, i1, ZXWireType::H);
  diag.add_wire(g3, o4, ZXWireType::H);
  diag.add_wire(g4, i3, ZXWireType::H);
  diag.add_wire(g4, inter1, ZXWireType::H);
  diag.add_wire(g5, i2, ZXWireType::H);
  diag.add_wire(g5, inter1, ZXWireType::H);
  diag.add_wire(g6, i2, ZXWireType::H);
  diag.add_wire(g6, i3, ZXWireType::H);
  diag.add_wire(g6, inter1, ZXWireType::H);
  diag.add_wire(g7, i1, ZXWireType::H);
  diag.add_wire(g7, o4, ZXWireType::H);
  diag.add_wire(g8, i0, ZXWireType::H);
  diag.add_wire(g8, o4, ZXWireType::H);
  diag.add_wire(g9, i2, ZXWireType::H);
  diag.add_wire(g9, i3, ZXWireType::H);

  // Output wires
  diag.add_wire(o0, o0ext, ZXWireType::H);
  diag.add_wire(o1, o1ext, ZXWireType::H);
  diag.add_wire(o2, o2ext, ZXWireType::H);
  diag.add_wire(o4, o4ext, ZXWireType::H);
  diag.add_wire(o0ext, outs.at(0));
  diag.add_wire(o1ext, outs.at(1));
  diag.add_wire(o2ext, outs.at(2));
  diag.add_wire(o3, outs.at(3));
  diag.add_wire(o4ext, outs.at(4));
  REQUIRE_NOTHROW(Flow::identify_pauli_flow(diag));
  Circuit c = zx_to_circuit(diag);
  REQUIRE_NOTHROW(c.assert_valid());
  CHECK(c.n_qubits() == 5);
}

}  // namespace tket::zx::test_ZXExtraction
