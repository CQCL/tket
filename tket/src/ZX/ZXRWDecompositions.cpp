// Copyright 2019-2021 Cambridge Quantum Computing
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
  std::list<ZXVert> to_decompose;
  BGL_FORALL_VERTICES(v, *diag.graph, ZXGraph) {
    if (diag.get_zxtype(v) == ZXType::ZXBox) to_decompose.push_back(v);
  }
  bool success = !to_decompose.empty();
  while (!to_decompose.empty()) {
    ZXVert box = to_decompose.front();
    to_decompose.pop_front();
    const ZXBox& zxb = diag.get_vertex_ZXGen<ZXBox>(box);
    const ZXDiagram& inner = *zxb.get_diagram();
    // Substitute diagram in place of box
    auto iso = diag.copy_graph(inner, false);
    std::map<Wire, Wire> incident_map;
    for (const Wire& w : diag.adj_wires(box)) {
      WireProperties wp = diag.get_wire_info(w);
      ZXVert s = diag.source(w);
      ZXVert t = diag.target(w);
      if (s == box) {
        ZXVert b = iso.first.at(inner.boundary.at(*wp.source_port));
        Wire bw = diag.adj_wires(b).at(0);
        if (diag.end_of(bw, b) == WireEnd::Source) {
          s = diag.target(bw);
          wp.source_port = wp.target_port;
        } else {
          s = diag.source(bw);
          wp.source_port = wp.source_port;
        }
      }
      if (t == box) {
        ZXVert b = iso.first.at(inner.boundary.at(*wp.target_port));
        Wire bw = diag.adj_wires(b).at(0);
        if (diag.end_of(bw, b) == WireEnd::Source) {
          t = diag.target(bw);
          wp.target_port = wp.target_port;
        } else {
          t = diag.source(bw);
          wp.target_port = wp.source_port;
        }
      }
      Wire new_w = diag.add_wire(s, t, wp);
      incident_map.emplace(w, new_w);
    }
    // Remove box and boundaries, and add new boxes to recursively decompose
    diag.remove_vertex(box);
    for (const std::pair<const ZXVert, ZXVert>& pair : iso.first) {
      ZXType type = diag.get_zxtype(pair.second);
      if (is_boundary_type(type))
        diag.remove_vertex(pair.second);
      else if (type == ZXType::ZXBox)
        to_decompose.push_back(pair.second);
    }
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

}  // namespace zx

}  // namespace tket
