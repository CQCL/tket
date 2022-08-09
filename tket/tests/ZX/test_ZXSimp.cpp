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
namespace test_ZXSimp {

SCENARIO("Testing graph state simplification") {
  /**
   * Diagram 1: testing simplification on graph states
   * This diagram follows from section A of:
   *  \ref https://arxiv.org/pdf/1902.03178.pdf
   **/
  ZXDiagram diag1(4, 4, 0, 0);
  ZXVertVec d1_in = diag1.get_boundary(ZXType::Input);
  ZXVertVec d1_out = diag1.get_boundary(ZXType::Output);

  ZXVert c11 = diag1.add_vertex(ZXType::ZSpider, 1.5);
  ZXVert c12 = diag1.add_vertex(ZXType::ZSpider, 0.5);
  ZXVert c13 = diag1.add_vertex(ZXType::ZSpider);
  ZXVert c14 = diag1.add_vertex(ZXType::XSpider);
  ZXVert c15 = diag1.add_vertex(ZXType::ZSpider, 0.25);
  ZXVert c21 = diag1.add_vertex(ZXType::ZSpider, 0.5);
  ZXVert c22 = diag1.add_vertex(ZXType::ZSpider);
  ZXVert c23 = diag1.add_vertex(ZXType::ZSpider);
  ZXVert c24 = diag1.add_vertex(ZXType::ZSpider, 0.25);
  ZXVert c25 = diag1.add_vertex(ZXType::ZSpider);
  ZXVert c31 = diag1.add_vertex(ZXType::XSpider);
  ZXVert c32 = diag1.add_vertex(ZXType::XSpider);
  ZXVert c33 = diag1.add_vertex(ZXType::ZSpider, 0.5);
  ZXVert c34 = diag1.add_vertex(ZXType::ZSpider, 0.5);
  ZXVert c35 = diag1.add_vertex(ZXType::XSpider);
  ZXVert c41 = diag1.add_vertex(ZXType::ZSpider);
  ZXVert c42 = diag1.add_vertex(ZXType::ZSpider);
  ZXVert c43 = diag1.add_vertex(ZXType::ZSpider, 1.5);
  ZXVert c44 = diag1.add_vertex(ZXType::XSpider, 1.0);
  ZXVert c45 = diag1.add_vertex(ZXType::ZSpider, 0.5);
  ZXVert c46 = diag1.add_vertex(ZXType::XSpider, 1.0);

  Wire cw11 = diag1.add_wire(d1_in[0], c11, ZXWireType::Basic);
  Wire cw12 = diag1.add_wire(c11, c12, ZXWireType::H);
  Wire cw13 = diag1.add_wire(c12, c13, ZXWireType::Basic);
  Wire cw13_ = diag1.add_wire(c13, c41, ZXWireType::H);
  Wire cw14 = diag1.add_wire(c13, c14, ZXWireType::Basic);
  Wire cw14_ = diag1.add_wire(c14, c42, ZXWireType::Basic);
  Wire cw15 = diag1.add_wire(c14, c15, ZXWireType::H);
  Wire cw16 = diag1.add_wire(c15, d1_out[0], ZXWireType::H);

  Wire cw21 = diag1.add_wire(d1_in[1], c21, ZXWireType::Basic);
  Wire cw22 = diag1.add_wire(c21, c22, ZXWireType::Basic);
  Wire cw22_ = diag1.add_wire(c22, c31, ZXWireType::Basic);
  Wire cw23 = diag1.add_wire(c22, c23, ZXWireType::H);
  Wire cw23_ = diag1.add_wire(c23, c32, ZXWireType::Basic);
  Wire cw24 = diag1.add_wire(c23, c24, ZXWireType::Basic);
  Wire cw25 = diag1.add_wire(c24, c25, ZXWireType::H);
  Wire cw25_ = diag1.add_wire(c25, c35, ZXWireType::Basic);
  Wire cw26 = diag1.add_wire(d1_out[1], c25, ZXWireType::Basic);

  Wire cw31 = diag1.add_wire(d1_in[2], c31, ZXWireType::Basic);
  Wire cw32 = diag1.add_wire(c31, c32, ZXWireType::Basic);
  Wire cw33 = diag1.add_wire(c32, c33, ZXWireType::Basic);
  Wire cw34 = diag1.add_wire(c33, c34, ZXWireType::H);
  Wire cw35 = diag1.add_wire(c34, c35, ZXWireType::Basic);
  Wire cw36 = diag1.add_wire(c35, d1_out[2], ZXWireType::Basic);

  Wire cw41 = diag1.add_wire(d1_in[3], c41, ZXWireType::H);
  Wire cw42 = diag1.add_wire(c41, c42, ZXWireType::Basic);
  Wire cw43 = diag1.add_wire(c42, c43, ZXWireType::H);
  Wire cw44 = diag1.add_wire(c43, c44, ZXWireType::Basic);
  Wire cw45 = diag1.add_wire(c44, c45, ZXWireType::Basic);
  Wire cw46 = diag1.add_wire(c45, c46, ZXWireType::Basic);
  Wire cw47 = diag1.add_wire(c46, d1_out[3], ZXWireType::Basic);

  REQUIRE_NOTHROW(diag1.check_validity());

  /** Apply rewrites to diagram 1 to turn into a graph-like form **/

  Rewrite::red_to_green().apply(diag1);
  Rewrite::spider_fusion().apply(diag1);
  Rewrite::parallel_h_removal().apply(diag1);
  Rewrite::io_extension().apply(diag1);
  Rewrite::separate_boundaries().apply(diag1);

  /**
   * Graph simplification via Pauli & Clifford removal
   * We perform the full simpliciation procedure described in theorem 5.4
   * of \ref https://arxiv.org/pdf/1902.03178.pdf
   **/

  CHECK(Rewrite::remove_interior_cliffords().apply(diag1));
  CHECK_FALSE(Rewrite::extend_at_boundary_paulis().apply(
      diag1));  // If remove_interior_cliffords is exhaustive, this should not
                // need to be applied
  CHECK(Rewrite::remove_interior_paulis().apply(diag1));

  CHECK_FALSE(Rewrite::parallel_h_removal().apply(diag1));
}

}  // namespace test_ZXSimp
}  // namespace zx
}  // namespace tket
