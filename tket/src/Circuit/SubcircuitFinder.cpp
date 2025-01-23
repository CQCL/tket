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

#include <algorithm>
#include <boost/graph/transitive_closure.hpp>
#include <cstddef>
#include <iterator>
#include <optional>
#include <tkassert/Assert.hpp>
#include <vector>

#include "tket/Circuit/Circuit.hpp"
#include "tket/Circuit/DAGDefs.hpp"
#include "tket/OpType/EdgeType.hpp"

namespace tket {

// Representation of a connected convex subcircuit.
typedef struct {
  VertexSet verts;  // all vertices in the subcircuit
  VertexSet preds;  // all predecessors of vertices in the subcircuit that are
                    // not themselves in it
  VertexSet succs;  // all successors of vertices in the subcircuit that are not
                    // themselves in it
} subcircuit_info_t;

// Test whether the union of two connected convex subcircuits is connected.
static bool union_is_connected(
    const subcircuit_info_t &subcircuit_info0,
    const subcircuit_info_t &subcircuit_info1) {
  const VertexSet &verts0 = subcircuit_info0.verts;
  const VertexSet &verts1 = subcircuit_info1.verts;
  for (const Vertex &v : subcircuit_info0.succs) {
    if (verts1.contains(v)) {
      return true;
    }
  }
  for (const Vertex &v : subcircuit_info1.succs) {
    if (verts0.contains(v)) {
      return true;
    }
  }
  return false;
}

static VertexSet set_diff(const VertexSet &A, const VertexSet &B) {
  VertexSet C;
  for (const Vertex &v : A) {
    if (!B.contains(v)) {
      C.insert(v);
    }
  }
  return C;
}

// Return the union of two disjoint convex connected subcircuits, assuming that
// the union is convex and connected.
static subcircuit_info_t convex_union(
    const subcircuit_info_t &subcircuit_info0,
    const subcircuit_info_t &subcircuit_info1) {
  const VertexSet &verts0 = subcircuit_info0.verts;
  const VertexSet &preds0 = subcircuit_info0.preds;
  const VertexSet &succs0 = subcircuit_info0.succs;
  const VertexSet &verts1 = subcircuit_info1.verts;
  const VertexSet &preds1 = subcircuit_info1.preds;
  const VertexSet &succs1 = subcircuit_info1.succs;

  // verts = verts0 ∪ verts1
  VertexSet verts = verts0;
  verts.insert(verts1.begin(), verts1.end());

  // preds = (preds0 ∪ preds1) ∖ verts
  VertexSet preds01 = preds0;
  preds01.insert(preds1.begin(), preds1.end());
  VertexSet preds = set_diff(preds01, verts);

  // succs = (succs0 ∪ succs1) ∖ verts
  VertexSet succs01 = succs0;
  succs01.insert(succs1.begin(), succs1.end());
  VertexSet succs = set_diff(succs01, verts);

  return {verts, preds, succs};
}

static std::set<std::pair<Vertex, Vertex>> order_relations(
    Circuit *circ /*const*/) {
  // Put the vertices in reverse topological order:
  VertexVec verts;
  boost::topological_sort(circ->dag, std::back_inserter(verts));
  // Construct a map v --> {all vertices in the future of v}:
  std::map<Vertex, VertexSet> futures;
  for (const Vertex &v : verts) {
    VertexSet v_futures = {v};
    for (const Vertex &w : circ->get_successors(v)) {
      const VertexSet &w_futures = futures[w];
      v_futures.insert(w_futures.begin(), w_futures.end());
    }
    futures[v] = v_futures;
  }
  // Construct the set of order relations:
  std::set<std::pair<Vertex, Vertex>> rels;
  for (const auto &v_ws : futures) {
    const Vertex &v = v_ws.first;
    for (const Vertex &w : v_ws.second) {
      rels.insert({v, w});
    }
  }
  return rels;
}

// Helper class for finding connected convex subcircuits.
class SubcircuitFinder {
 public:
  SubcircuitFinder(Circuit *circ /*const*/)
      : circ_(circ), order_relations_(order_relations(circ)) {}
  std::vector<VertexSet> find_subcircuits(
      std::function<bool(Op_ptr)> criterion) {
    // Find a maximal partition of the vertices satisfying the criterion into
    // connected convex subcircuits. We use a greedy algorithm, beginning with
    // the trivial partition into singletons, and then repeatedly looking for a
    // pair that we can merge.
    std::vector<subcircuit_info_t> subcircuit_infos;
    // Form a trivial partition, one node per subcircuit:
    BGL_FORALL_VERTICES(v, circ_->dag, DAG) {
      if (criterion(circ_->get_Op_ptr_from_Vertex(v))) {
        VertexVec preds = circ_->get_predecessors(v);
        VertexVec succs = circ_->get_successors(v);
        subcircuit_infos.push_back(
            {{v}, {preds.begin(), preds.end()}, {succs.begin(), succs.end()}});
      }
    }
    while (true) {
      // If possible, find a mergeable pair of subcircuits and merge them.
      std::optional<std::pair<std::size_t, std::size_t>> indices =
          find_mergeable_pair(subcircuit_infos);
      if (!indices.has_value()) {
        break;
      } else {
        std::size_t i0 = indices->first;
        std::size_t i1 = indices->second;
        TKET_ASSERT(i0 < i1);
        subcircuit_info_t subcircuit_info0 = subcircuit_infos[i0];
        subcircuit_info_t subcircuit_info1 = subcircuit_infos[i1];
        // Must erase in this order, since i0 < i1:
        subcircuit_infos.erase(subcircuit_infos.begin() + i1);
        subcircuit_infos.erase(subcircuit_infos.begin() + i0);
        subcircuit_infos.push_back(
            convex_union(subcircuit_info0, subcircuit_info1));
      }
    }
    std::vector<VertexSet> subcircuits;
    for (const subcircuit_info_t &subcircuit_info : subcircuit_infos) {
      subcircuits.push_back(subcircuit_info.verts);
    }
    return subcircuits;
  }

 private:
  // Test whether the union of two disjoint connected convex subcircuits is
  // convex.
  bool union_is_convex(
      const subcircuit_info_t &subcircuit_info0,
      const subcircuit_info_t &subcircuit_info1) {
    const VertexSet &preds0 = subcircuit_info0.preds;
    const VertexSet &succs0 = subcircuit_info0.succs;
    const VertexSet &preds1 = subcircuit_info1.preds;
    const VertexSet &succs1 = subcircuit_info1.succs;
    for (const Vertex &v0 : succs0) {
      for (const Vertex &v1 : preds1) {
        if (order_relations_.contains({v0, v1})) {
          return false;
        }
      }
    }
    for (const Vertex &v1 : succs1) {
      for (const Vertex &v0 : preds0) {
        if (order_relations_.contains({v1, v0})) {
          return false;
        }
      }
    }
    return true;
  }

  // Given a vector of disjoint connected convex subcircuits, look for a pair
  // whose union is connected and convex, and return the indices of such a pair
  // if it exists.
  std::optional<std::pair<std::size_t, std::size_t>> find_mergeable_pair(
      const std::vector<subcircuit_info_t> &subcircuit_infos) {
    std::size_t n_subcircuits = subcircuit_infos.size();
    for (std::size_t i0 = 0; i0 < n_subcircuits; i0++) {
      for (std::size_t i1 = i0 + 1; i1 < n_subcircuits; i1++) {
        const subcircuit_info_t &subcircuit_info0 = subcircuit_infos[i0];
        const subcircuit_info_t &subcircuit_info1 = subcircuit_infos[i1];
        if (union_is_connected(subcircuit_info0, subcircuit_info1) &&
            union_is_convex(subcircuit_info0, subcircuit_info1)) {
          return std::pair<std::size_t, std::size_t>{i0, i1};
        }
      }
    }
    return std::nullopt;
  }
  const Circuit *circ_;
  std::set<std::pair<Vertex, Vertex>> order_relations_;
};

std::vector<VertexSet> Circuit::get_subcircuits(
    std::function<bool(Op_ptr)> criterion) /*const*/ {
  index_vertices();
  SubcircuitFinder finder(this);
  return finder.find_subcircuits(criterion);
}

Subcircuit Circuit::make_subcircuit(VertexSet verts) const {
  const std::map<Edge, UnitID> unitmap = edge_unit_map();
  EdgeSet q_ins, q_outs, c_ins, c_outs;
  for (const Vertex &v : verts) {
    for (const Edge &v_in : get_in_edges(v)) {
      if (!verts.contains(source(v_in))) {
        switch (get_edgetype(v_in)) {
          case EdgeType::Quantum:
            q_ins.insert(v_in);
            break;
          case EdgeType::Classical:
            c_ins.insert(v_in);
            break;
          default:
            // Boolean and WASM edges ignored.
            // TODO What to do? Throw?
            break;
        }
      }
    }
    for (const Edge &v_out : get_all_out_edges(v)) {
      if (!verts.contains(target(v_out))) {
        switch (get_edgetype(v_out)) {
          case EdgeType::Quantum:
            q_outs.insert(v_out);
            break;
          case EdgeType::Classical:
            c_outs.insert(v_out);
            break;
          default:
            // Boolean and WASM edges ignored.
            // TODO What to do? Throw?
            break;
        }
      }
    }
  }
  EdgeSet crs;
  for (const Edge &e : c_outs) {
    Vertex v = source(e);
    for (const Edge &b_out : get_out_edges_of_type(v, EdgeType::Boolean)) {
      crs.insert(b_out);
    }
  }
  // Convert to vectors, taking care of order:
  TKET_ASSERT(q_ins.size() == q_outs.size());
  TKET_ASSERT(c_ins.size() == c_outs.size());
  EdgeVec q_ins_vec{q_ins.begin(), q_ins.end()};
  EdgeVec q_outs_vec{q_outs.begin(), q_outs.end()};
  EdgeVec c_ins_vec{c_ins.begin(), c_ins.end()};
  EdgeVec c_outs_vec{c_outs.begin(), c_outs.end()};
  EdgeVec crs_vec{crs.begin(), crs.end()};
  std::sort(
      q_ins_vec.begin(), q_ins_vec.end(),
      [&unitmap](const Edge &e0, const Edge &e1) {
        return unitmap.at(e0) < unitmap.at(e1);
      });
  std::sort(
      q_outs_vec.begin(), q_outs_vec.end(),
      [&unitmap](const Edge &e0, const Edge &e1) {
        return unitmap.at(e0) < unitmap.at(e1);
      });
  std::sort(
      c_ins_vec.begin(), c_ins_vec.end(),
      [&unitmap](const Edge &e0, const Edge &e1) {
        return unitmap.at(e0) < unitmap.at(e1);
      });
  std::sort(
      c_outs_vec.begin(), c_outs_vec.end(),
      [&unitmap](const Edge &e0, const Edge &e1) {
        return unitmap.at(e0) < unitmap.at(e1);
      });
  return {q_ins_vec, q_outs_vec, c_ins_vec, c_outs_vec, crs_vec, verts};
}

}  // namespace tket
