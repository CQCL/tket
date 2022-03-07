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

#include "Utils/GraphHeaders.hpp"
#include "ZX/Rewrite.hpp"

namespace tket {

namespace zx {

bool Rewrite::decompose_boxes_fun(ZXDiagram& diag) {
  bool success = false;
  std::list<ZXVert> to_decompose;
  BGL_FORALL_VERTICES(v, *diag.graph, ZXGraph) {
    if (diag.get_zxtype(v) == ZXType::ZXBox) to_decompose.push_back(v);
  }
  for (const ZXVert& box : to_decompose) {
    const ZXBox& zxb = diag.get_vertex_ZXGen<ZXBox>(box);
    ZXDiagram inner(*zxb.get_diagram());
    decompose_boxes_fun(inner);
    std::vector<std::pair<Wire, WireEnd>> sub_boundary(zxb.n_ports());
    for (const Wire& w : diag.adj_wires(box)) {
      if (diag.source(w) == box)
        sub_boundary.at(*diag.source_port(w)) = {w, WireEnd::Source};
      if (diag.target(w) == box)
        sub_boundary.at(*diag.target_port(w)) = {w, WireEnd::Target};
    }
    ZXDiagram::Subdiagram sub{sub_boundary, {box}};
    diag.substitute(inner, sub);
    success = true;
  }
  return success;
}

Rewrite Rewrite::decompose_boxes() { return Rewrite(decompose_boxes_fun); }

bool Rewrite::basic_wires_fun(ZXDiagram& diag) {
  ZXGen_ptr qhad =
      std::make_shared<const PhasedGen>(ZXType::Hbox, -1, QuantumType::Quantum);
  ZXGen_ptr chad = std::make_shared<const PhasedGen>(
      ZXType::Hbox, -1, QuantumType::Classical);
  WireVec targets;
  BGL_FORALL_EDGES(w, *diag.graph, ZXGraph) {
    if (diag.get_wire_type(w) == ZXWireType::H) targets.push_back(w);
  }
  for (const Wire& w : targets) {
    WireProperties wp = diag.get_wire_info(w);
    ZXGen_ptr had = (wp.qtype == QuantumType::Quantum) ? qhad : chad;
    ZXVert h = diag.add_vertex(had);
    ZXVert s = diag.source(w);
    ZXVert t = diag.target(w);
    wp.type = ZXWireType::Basic;
    WireProperties wp2 = wp;
    wp.target_port = std::nullopt;
    wp2.source_port = std::nullopt;
    diag.add_wire(s, h, wp);
    diag.add_wire(h, t, wp);
    diag.remove_wire(w);
  }
  return !targets.empty();
}

Rewrite Rewrite::basic_wires() { return Rewrite(basic_wires_fun); }

bool Rewrite::rebase_to_zx_fun(ZXDiagram& diag) {
  ZXVertVec verts;
  BGL_FORALL_VERTICES(v, *diag.graph, ZXGraph) {
    ZXType t = diag.get_zxtype(v);
    if (!is_boundary_type(t) && !is_spider_type(t)) verts.push_back(v);
  }
  for (const ZXVert& v : verts) {
    switch (diag.get_zxtype(v)) {
      case ZXType::Hbox: {
        // TODO:: Find a paper with a decomposition
        break;
      }
      case ZXType::XY: {
        ZXGen_ptr vgen = diag.get_vertex_ZXGen_ptr(v);
        const PhasedGen& vg = static_cast<const PhasedGen&>(*vgen);
        ZXGen_ptr new_gen = ZXGen::create_gen(
            ZXType::ZSpider, -vg.get_param(), *vg.get_qtype());
        diag.set_vertex_ZXGen_ptr(v, new_gen);
        break;
      }
      case ZXType::XZ: {
        ZXGen_ptr vgen = diag.get_vertex_ZXGen_ptr(v);
        const PhasedGen& vg = static_cast<const PhasedGen&>(*vgen);
        ZXVert phase_v =
            diag.add_vertex(ZXType::XSpider, vg.get_param(), *vg.get_qtype());
        ZXGen_ptr new_gen =
            ZXGen::create_gen(ZXType::ZSpider, Expr(0.5), *vg.get_qtype());
        diag.set_vertex_ZXGen_ptr(v, new_gen);
        diag.add_wire(v, phase_v, ZXWireType::Basic, *vg.get_qtype());
        break;
      }
      case ZXType::YZ: {
        ZXGen_ptr vgen = diag.get_vertex_ZXGen_ptr(v);
        const PhasedGen& vg = static_cast<const PhasedGen&>(*vgen);
        ZXVert phase_v =
            diag.add_vertex(ZXType::XSpider, vg.get_param(), *vg.get_qtype());
        ZXGen_ptr new_gen = ZXGen::create_gen(ZXType::ZSpider, *vg.get_qtype());
        diag.set_vertex_ZXGen_ptr(v, new_gen);
        diag.add_wire(v, phase_v, ZXWireType::Basic, *vg.get_qtype());
        break;
      }
      case ZXType::PX: {
        ZXGen_ptr vgen = diag.get_vertex_ZXGen_ptr(v);
        const CliffordGen& vg = static_cast<const CliffordGen&>(*vgen);
        ZXGen_ptr new_gen = ZXGen::create_gen(
            ZXType::ZSpider, vg.get_param() ? Expr(1.) : Expr(0.),
            *vg.get_qtype());
        diag.set_vertex_ZXGen_ptr(v, new_gen);
        break;
      }
      case ZXType::PY: {
        ZXGen_ptr vgen = diag.get_vertex_ZXGen_ptr(v);
        const CliffordGen& vg = static_cast<const CliffordGen&>(*vgen);
        ZXGen_ptr new_gen = ZXGen::create_gen(
            ZXType::ZSpider, vg.get_param() ? Expr(0.5) : Expr(-0.5),
            *vg.get_qtype());
        diag.set_vertex_ZXGen_ptr(v, new_gen);
        break;
      }
      case ZXType::PZ: {
        ZXGen_ptr vgen = diag.get_vertex_ZXGen_ptr(v);
        const CliffordGen& vg = static_cast<const CliffordGen&>(*vgen);
        ZXVert phase_v = diag.add_vertex(
            ZXType::XSpider, vg.get_param() ? Expr(1.) : Expr(0.),
            *vg.get_qtype());
        ZXGen_ptr new_gen = ZXGen::create_gen(ZXType::ZSpider, *vg.get_qtype());
        diag.set_vertex_ZXGen_ptr(v, new_gen);
        diag.add_wire(v, phase_v, ZXWireType::Basic, *vg.get_qtype());
        break;
      }
      case ZXType::Triangle: {
        // Decomposition of the triangle given in Lemma 3.3, "Completeness of
        // the ZX-calculus for Pure Qubit Clifford+T Quantum Mechanics", K. Feng
        // Ng, Q. Wang, 2018
        ZXDiagram tri(0, 0, 0, 0);
        QuantumType qt = *diag.get_qtype(v);
        ZXVert in = tri.add_vertex(ZXType::Input, qt);
        ZXVert out = tri.add_vertex(ZXType::Output, qt);
        tri.multiply_scalar(
            (qt == QuantumType::Quantum) ? Expr(2.)
                                         : Expr(SymEngine::sqrt(Expr(2.))));
        ZXVert split = tri.add_vertex(ZXType::XSpider, qt);
        ZXVert lrz = tri.add_vertex(ZXType::ZSpider, -0.25, qt);
        ZXVert rrz = tri.add_vertex(ZXType::ZSpider, 0.25, qt);
        ZXVert laxis = tri.add_vertex(ZXType::XSpider, qt);
        ZXVert raxis = tri.add_vertex(ZXType::XSpider, qt);
        ZXVert lph = tri.add_vertex(ZXType::ZSpider, -0.25, qt);
        ZXVert rph = tri.add_vertex(ZXType::ZSpider, 0.25, qt);
        ZXVert merge = tri.add_vertex(ZXType::ZSpider, qt);
        tri.add_wire(in, split, ZXWireType::Basic, qt);
        tri.add_wire(split, lrz, ZXWireType::Basic, qt);
        tri.add_wire(split, rrz, ZXWireType::Basic, qt);
        tri.add_wire(lrz, laxis, ZXWireType::Basic, qt);
        tri.add_wire(rrz, raxis, ZXWireType::Basic, qt);
        tri.add_wire(laxis, lph, ZXWireType::Basic, qt);
        tri.add_wire(raxis, rph, ZXWireType::Basic, qt);
        tri.add_wire(laxis, merge, ZXWireType::Basic, qt);
        tri.add_wire(raxis, merge, ZXWireType::Basic, qt);
        tri.add_wire(merge, out, ZXWireType::Basic, qt);

        Wire w0 = diag.wire_at_port(v, 0);
        Wire w1 = diag.wire_at_port(v, 1);
        WireEnd we0, we1;
        if (w0 == w1) {
          if (diag.source_port(w0) == 0) {
            we0 = WireEnd::Source;
            we1 = WireEnd::Target;
          } else {
            we0 = WireEnd::Target;
            we1 = WireEnd::Source;
          }
        } else {
          we0 = diag.end_of(w0, v);
          we1 = diag.end_of(w1, v);
        }
        ZXDiagram::Subdiagram sub{{{w0, we0}, {w1, we1}}, {v}};
        diag.substitute(tri, sub);
        break;
      }
      default:
        break;
    }
  }
  return !verts.empty();
}

Rewrite Rewrite::rebase_to_zx() { return Rewrite(rebase_to_zx_fun); }

bool Rewrite::rebase_to_mbqc_fun(ZXDiagram& diag) {
  ZXVertVec verts;
  BGL_FORALL_VERTICES(v, *diag.graph, ZXGraph) {
    ZXType t = diag.get_zxtype(v);
    if (!is_boundary_type(t) && !is_MBQC_type(t)) verts.push_back(v);
  }
  for (const ZXVert& v : verts) {
    switch (diag.get_zxtype(v)) {
      case ZXType::ZSpider: {
        ZXGen_ptr vgen = diag.get_vertex_ZXGen_ptr(v);
        const PhasedGen& vg = static_cast<const PhasedGen&>(*vgen);
        ZXGen_ptr new_gen =
            ZXGen::create_gen(ZXType::XY, -vg.get_param(), *vg.get_qtype());
        diag.set_vertex_ZXGen_ptr(v, new_gen);
        break;
      }
      case ZXType::XSpider: {
        ZXGen_ptr vgen = diag.get_vertex_ZXGen_ptr(v);
        const PhasedGen& vg = static_cast<const PhasedGen&>(*vgen);
        ZXGen_ptr new_gen =
            ZXGen::create_gen(ZXType::XY, -vg.get_param(), *vg.get_qtype());
        diag.set_vertex_ZXGen_ptr(v, new_gen);
        for (const Wire& w : diag.adj_wires(v)) {
          ZXWireType wt = diag.get_wire_type(w);
          diag.set_wire_type(
              w, (wt == ZXWireType::Basic) ? ZXWireType::H : ZXWireType::Basic);
        }
        break;
      }
      case ZXType::Hbox: {
        // TODO:: Find a paper with a decomposition
        break;
      }
      case ZXType::Triangle: {
        ZXDiagram tri(0, 0, 0, 0);
        QuantumType qt = *diag.get_qtype(v);
        ZXVert in = tri.add_vertex(ZXType::Input, qt);
        ZXVert out = tri.add_vertex(ZXType::Output, qt);
        ZXVert t = tri.add_vertex(ZXType::Triangle, qt);
        tri.add_wire(in, t, ZXWireType::Basic, qt, std::nullopt, 0);
        tri.add_wire(out, t, ZXWireType::Basic, qt, std::nullopt, 1);
        rebase_to_zx_fun(tri);
        rebase_to_mbqc_fun(tri);

        Wire w0 = diag.wire_at_port(v, 0);
        Wire w1 = diag.wire_at_port(v, 1);
        WireEnd we0, we1;
        if (w0 == w1) {
          if (diag.source_port(w0) == 0) {
            we0 = WireEnd::Source;
            we1 = WireEnd::Target;
          } else {
            we0 = WireEnd::Target;
            we1 = WireEnd::Source;
          }
        } else {
          we0 = diag.end_of(w0, v);
          we1 = diag.end_of(w1, v);
        }
        ZXDiagram::Subdiagram sub{{{w0, we0}, {w1, we1}}, {v}};
        diag.substitute(tri, sub);
        break;
      }
      default:
        break;
    }
  }
  return !verts.empty();
}

Rewrite Rewrite::rebase_to_mbqc() { return Rewrite(rebase_to_mbqc_fun); }

}  // namespace zx

}  // namespace tket
