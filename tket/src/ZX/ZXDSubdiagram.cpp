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
#include "ZX/ZXDiagram.hpp"

namespace tket {

namespace zx {

ZXDiagram::Subdiagram::Subdiagram() : boundary_(), verts_() {}

ZXDiagram::Subdiagram::Subdiagram(
    const std::vector<std::pair<Wire, WireEnd>>& cut, const ZXVertSeqSet& verts)
    : boundary_(cut), verts_(verts) {}

void ZXDiagram::Subdiagram::check_validity(const ZXDiagram& diag) const {
  std::set<std::pair<Wire, WireEnd>> boundary_lookup;
  for (const std::pair<Wire, WireEnd>& w : boundary_) {
    auto inserted = boundary_lookup.insert(w);
    if (!inserted.second)
      throw ZXError(
          "Malformed ZX Subdiagram: Wire appears multiple times in boundary");
    if (verts_.find(diag.vertex_at_end(w.first, w.second)) == verts_.end())
      throw ZXError(
          "Malformed ZX Subdiagram: Vertex adjacent to boundary is not in "
          "vertex set");
  }
  for (const ZXVert& v : verts_) {
    if (is_boundary_type(diag.get_zxtype(v)))
      throw ZXError("Malformed ZX Subdiagram: Contains a boundary vertex");
    for (const Wire& w : diag.adj_wires(v)) {
      ZXVert n = diag.other_end(w, v);
      if (verts_.find(n) == verts_.end()) {
        if (boundary_lookup.find({w, diag.end_of(w, v)}) ==
            boundary_lookup.end())
          throw ZXError("Malformed ZX Subdiagram: subdiagram is not closed");
      } else {
        if ((boundary_lookup.find({w, WireEnd::Source}) ==
             boundary_lookup.end()) ^
            (boundary_lookup.find({w, WireEnd::Target}) ==
             boundary_lookup.end()))
          throw ZXError(
              "Malformed ZX Subdiagram: wire between two interior vertices "
              "contains one boundary");
      }
    }
  }
}

ZXDiagram ZXDiagram::Subdiagram::to_diagram(const ZXDiagram& orig) const {
  ZXDiagram diag;
  std::map<ZXVert, ZXVert> vert_iso;
  std::map<std::pair<Wire, WireEnd>, ZXVert> bound_iso;
  for (const std::pair<Wire, WireEnd>& bw : boundary_) {
    ZXVert bv = diag.add_vertex(ZXType::Open, orig.get_qtype(bw.first));
    diag.boundary.push_back(bv);
    bound_iso.insert({bw, bv});
  }
  for (const ZXVert& ov : verts_) {
    ZXVert v = diag.add_vertex(orig.get_vertex_ZXGen_ptr(ov));
    vert_iso.insert({ov, v});
    for (const Wire& w : orig.adj_wires(ov)) {
      bool internal = true;
      if (orig.source(w) == ov) {
        auto found_bound = bound_iso.find({w, WireEnd::Source});
        if (found_bound != bound_iso.end()) {
          internal = false;
          WireProperties wp = orig.get_wire_info(w);
          diag.add_wire(
              v, found_bound->second, ZXWireType::Basic, wp.qtype,
              wp.source_port, std::nullopt);
        }
      }
      if (orig.target(w) == ov) {
        auto found_bound = bound_iso.find({w, WireEnd::Target});
        if (found_bound != bound_iso.end()) {
          internal = false;
          WireProperties wp = orig.get_wire_info(w);
          diag.add_wire(
              found_bound->second, v, ZXWireType::Basic, wp.qtype, std::nullopt,
              wp.target_port);
        }
      }
      if (internal) {
        ZXVert other = orig.other_end(w, ov);
        auto found_other = vert_iso.find(other);
        if (found_other != vert_iso.end()) {
          WireProperties wp = orig.get_wire_info(w);
          if (orig.source(w) == ov)
            diag.add_wire(v, found_other->second, wp);
          else
            diag.add_wire(found_other->second, v, wp);
        }
      }
    }
  }
  return diag;
}

void ZXDiagram::substitute(
    const ZXDiagram& to_insert, const Subdiagram& to_replace) {
  unsigned n_bounds = to_insert.boundary.size();
  if (n_bounds != to_replace.boundary_.size())
    throw ZXError(
        "ZXDiagram substitution error: boundary size of replacement does not "
        "fit size of subdiagram");

  std::pair<std::map<ZXVert, ZXVert>, std::map<Wire, Wire>> iso =
      copy_graph(to_insert, false);

  std::map<Wire, ZXVert> loops;
  for (unsigned i = 0; i < n_bounds; ++i) {
    std::pair<Wire, WireEnd> w = to_replace.boundary_.at(i);
    WireProperties wp = get_wire_info(w.first);
    ZXVert new_b = iso.first.at(to_insert.boundary.at(i));
    if (wp.qtype != get_qtype(new_b))
      throw ZXError(
          "ZXDiagram substitution error: QuantumType mismatch at a boundary");
    ZXVert to_connect =
        (w.second == WireEnd::Source) ? target(w.first) : source(w.first);
    if (to_replace.verts_.find(to_connect) != to_replace.verts_.end()) {
      // Connect two boundaries together
      auto inserted = loops.insert({w.first, new_b});
      if (!inserted.second) {
        // Already seen the other end, so actually connect the loop
        ZXVert other_b = inserted.first->second;
        Wire wb = adj_wires(new_b).at(0);
        ZXVert adj = other_end(wb, new_b);
        std::optional<unsigned> adj_port =
            (source(wb) == adj) ? source_port(wb) : target_port(wb);
        if (adj == other_b) {
          // Creates a wire-loop, resolve to a scalar
          bool hadamard =
              (wp.type == ZXWireType::H) ^ (get_wire_type(wb) == ZXWireType::H);
          if (hadamard)
            multiply_scalar(0.);
          else if (wp.qtype == QuantumType::Quantum)
            multiply_scalar(4.);
          else
            multiply_scalar(2.);
        } else {
          // Connect the adjacent vertices
          Wire wob = adj_wires(other_b).at(0);
          ZXVert other_adj = other_end(wob, other_b);
          std::optional<unsigned> other_adj_port =
              (source(wob) == other_adj) ? source_port(wob) : target_port(wob);
          bool hadamard = (wp.type == ZXWireType::H) ^
                          (get_wire_type(wb) == ZXWireType::H) ^
                          (get_wire_type(wob) == ZXWireType::H);
          add_wire(
              adj, other_adj, hadamard ? ZXWireType::H : ZXWireType::Basic,
              wp.qtype, adj_port, other_adj_port);
        }
        remove_vertex(new_b);
        remove_vertex(other_b);
      }
    } else {
      Wire wb = adj_wires(new_b).at(0);
      if (get_wire_type(wb) == ZXWireType::H) {
        if (wp.type == ZXWireType::Basic)
          wp.type = ZXWireType::H;
        else
          wp.type = ZXWireType::Basic;
      }
      ZXVert adj = other_end(wb, new_b);
      std::optional<unsigned> adj_port =
          (source(wb) == adj) ? source_port(wb) : target_port(wb);
      if (w.second == WireEnd::Source) {
        wp.source_port = adj_port;
        add_wire(adj, target(w.first), wp);
      } else {
        wp.target_port = adj_port;
        add_wire(source(w.first), adj, wp);
      }
      remove_vertex(new_b);
    }
  }
  for (const ZXVert& v : to_replace.verts_) {
    remove_vertex(v);
  }
}

}  // namespace zx

}  // namespace tket
