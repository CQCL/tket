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

/**
 * Helper method for local complementation and pivoting.
 * Checks that all neighbours of some vertex v are also ZSpiders.
 * If v is Classical, all neighbours must also be Classical.
 */
static bool can_complement_neighbourhood(
    const ZXDiagram& diag, QuantumType vqtype, const ZXVertVec& neighbours) {
  for (const ZXVert& n : neighbours) {
    if (diag.get_zxtype(n) != ZXType::ZSpider ||
        (vqtype == QuantumType::Classical &&
         *diag.get_qtype(n) == QuantumType::Quantum))
      return false;
  }
  return true;
}

bool Rewrite::remove_interior_cliffords_fun(ZXDiagram& diag) {
  if (!diag.is_graphlike()) return false;
  bool success = false;
  ZXVertSeqSet candidates;
  BGL_FORALL_VERTICES(v, *diag.graph, ZXGraph) { candidates.insert(v); }
  auto& view = candidates.get<TagSeq>();
  while (!candidates.empty()) {
    auto it = view.begin();
    ZXVert v = *it;
    view.erase(it);
    if (!diag.is_proper_clifford_spider(v)) continue;
    const PhasedGen& spid = diag.get_vertex_ZXGen<PhasedGen>(v);
    QuantumType vqtype = *spid.get_qtype();
    ZXVertVec neighbours = diag.neighbours(v);
    if (!can_complement_neighbourhood(diag, vqtype, neighbours)) continue;
    // Found an internal proper clifford spider on which we can perform local
    // complementation
    /**
     * Complement the neighbourhoods' edges and modify the phase information
     * on the neighbours.
     **/
    auto xi = neighbours.begin(), x_end = neighbours.end();
    for (; xi != x_end; ++xi) {
      for (auto yi = xi + 1; yi != x_end; ++yi) {
        // Don't add a doubled edge between classicals to preserve graph-like
        if (!(vqtype == QuantumType::Quantum &&
              *diag.get_qtype(*xi) == QuantumType::Classical &&
              *diag.get_qtype(*yi) == QuantumType::Classical)) {
          std::optional<Wire> wire = diag.wire_between(*xi, *yi);
          if (wire)
            diag.remove_wire(*wire);
          else
            diag.add_wire(*xi, *yi, ZXWireType::H, vqtype);
        }
      }
      const PhasedGen& xi_op = diag.get_vertex_ZXGen<PhasedGen>(*xi);
      // If `v` is Quantum, Classical neighbours will pick up both the +theta
      // and -theta phases, cancelling out
      if (vqtype == QuantumType::Quantum &&
          *xi_op.get_qtype() == QuantumType::Classical)
        continue;
      // Update phase information
      ZXGen_ptr xi_new_op = std::make_shared<const PhasedGen>(
          ZXType::ZSpider, xi_op.get_param() - spid.get_param(),
          *xi_op.get_qtype());
      diag.set_vertex_ZXGen_ptr(*xi, xi_new_op);
      candidates.insert(
          *xi);  // Changing the phase could introduce a new proper Clifford
    }
    diag.remove_vertex(v);
    success = true;
  }
  return success;
}

Rewrite Rewrite::remove_interior_cliffords() {
  return Rewrite(remove_interior_cliffords_fun);
}

static void add_phase_to_vertices(
    ZXDiagram& diag, const ZXVertSeqSet& verts, const Expr& phase) {
  for (const ZXVert& v : verts) {
    const PhasedGen& old_spid = diag.get_vertex_ZXGen<PhasedGen>(v);
    ZXGen_ptr new_spid = std::make_shared<const PhasedGen>(
        ZXType::ZSpider, old_spid.get_param() + phase, *old_spid.get_qtype());
    diag.set_vertex_ZXGen_ptr(v, new_spid);
  }
}

static void bipartite_complementation(
    ZXDiagram& diag, const ZXVertSeqSet& sa, const ZXVertSeqSet& sb,
    QuantumType qtype) {
  for (const ZXVert& a : sa.get<TagSeq>()) {
    for (const ZXVert& b : sb.get<TagSeq>()) {
      // Don't add a doubled edge between classicals to preserve graph-like
      if (!(qtype == QuantumType::Quantum &&
            *diag.get_qtype(a) == QuantumType::Classical &&
            *diag.get_qtype(b) == QuantumType::Classical)) {
        std::optional<Wire> wire = diag.wire_between(a, b);
        if (wire)
          diag.remove_wire(*wire);
        else
          diag.add_wire(a, b, ZXWireType::H, qtype);
      }
    }
  }
}

bool Rewrite::remove_interior_paulis_fun(ZXDiagram& diag) {
  if (!diag.is_graphlike()) return false;
  bool success = false;
  ZXVertSeqSet candidates;  // Need an indirect iterator as BGL_FORALL_VERTICES
                            // breaks when removing the current vertex
  BGL_FORALL_VERTICES(v, *diag.graph, ZXGraph) { candidates.insert(v); }
  auto& view = candidates.get<TagSeq>();
  while (!candidates.empty()) {
    auto it = view.begin();
    ZXVert v = *it;
    view.erase(it);
    // Check `v` is an interior Pauli
    if (!diag.is_pauli_spider(v)) continue;
    ZXVertVec v_ns = diag.neighbours(v);
    QuantumType vqtype = *diag.get_qtype(v);
    if (!can_complement_neighbourhood(diag, vqtype, v_ns)) continue;
    // Look for an interior Pauli neighbour
    bool pair_found = false;
    ZXVert u;
    ZXVertVec u_ns;
    for (const ZXVert& n : v_ns) {
      if (!diag.is_pauli_spider(n)) continue;
      ZXVertVec n_ns = diag.neighbours(n);
      QuantumType nqtype = *diag.get_qtype(n);
      if (can_complement_neighbourhood(diag, nqtype, n_ns)) {
        pair_found = true;
        u = n;
        u_ns = n_ns;
        break;
      }
    }
    if (!pair_found) continue;
    // Found a valid pair
    // Identify the three sets from the neighbourhoods of `u` and `v`
    ZXVertSeqSet excl_v{v_ns.begin(), v_ns.end()};
    excl_v.erase(u);
    ZXVertSeqSet excl_u, joint;
    auto& lookup_v = excl_v.get<TagKey>();
    for (const ZXVert& nu : u_ns) {
      if (lookup_v.find(nu) != lookup_v.end())
        joint.insert(nu);
      else
        excl_u.insert(nu);
    }
    excl_u.erase(v);
    excl_v.erase(joint.begin(), joint.end());
    const PhasedGen& v_spid = diag.get_vertex_ZXGen<PhasedGen>(v);
    const PhasedGen& u_spid = diag.get_vertex_ZXGen<PhasedGen>(u);

    add_phase_to_vertices(
        diag, joint, v_spid.get_param() + u_spid.get_param() + 1.);
    add_phase_to_vertices(diag, excl_u, v_spid.get_param());
    add_phase_to_vertices(diag, excl_v, u_spid.get_param());

    // Because `can_complement_neighbourhood` checks all neighbours,
    // v and u have the same QuantumType
    bipartite_complementation(diag, joint, excl_u, vqtype);
    bipartite_complementation(diag, joint, excl_v, vqtype);
    bipartite_complementation(diag, excl_u, excl_v, vqtype);

    diag.remove_vertex(u);
    diag.remove_vertex(v);
    candidates.erase(u);
    success = true;
  }
  return success;
}

Rewrite Rewrite::remove_interior_paulis() {
  return Rewrite(remove_interior_paulis_fun);
}

bool Rewrite::extend_at_boundary_paulis_fun(ZXDiagram& diag) {
  if (!diag.is_graphlike()) return false;
  bool success = false;
  for (const ZXVert& b : diag.get_boundary()) {
    // Valid ZX graph requires boundaries to have a unique neighbour
    Wire bw = diag.adj_wires(b).at(0);
    ZXVert u = diag.other_end(bw, b);
    if (!diag.is_pauli_spider(u)) continue;
    bool has_internal_pauli = false;
    for (const ZXVert& w : diag.neighbours(u)) {
      if (!diag.is_pauli_spider(w)) continue;
      bool interior = true;
      for (const ZXVert& wn : diag.neighbours(w)) {
        if (is_boundary_type(diag.get_zxtype(wn))) {
          interior = false;
          break;
        }
      }
      if (interior) {
        has_internal_pauli = true;
        break;
      }
    }
    if (!has_internal_pauli) continue;
    // We would like to pivot about (`u`, `w`) but `u` is by a boundary, so we
    // extend it
    ZXGen_ptr u_op = diag.get_vertex_ZXGen_ptr(u);
    QuantumType qtype = *u_op->get_qtype();
    ZXGen_ptr id =
        std::make_shared<const PhasedGen>(ZXType::ZSpider, 0., qtype);
    ZXVert z1 = diag.add_vertex(id);
    ZXVert z2 = diag.add_vertex(u_op);
    diag.add_wire(u, z1, ZXWireType::H, qtype);
    diag.add_wire(z1, z2, ZXWireType::H, qtype);
    diag.add_wire(z2, b, ZXWireType::Basic, qtype);
    diag.remove_wire(bw);
    diag.set_vertex_ZXGen_ptr(u, id);
    success = true;
  }
  return success;
}

Rewrite Rewrite::extend_at_boundary_paulis() {
  return Rewrite(extend_at_boundary_paulis_fun);
}

}  // namespace zx

}  // namespace tket
