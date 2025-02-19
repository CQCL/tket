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
  std::map<UnitID, Edge> out_lookup;
  EdgeSet ins, outs, b_future;
  for (const Vertex &v : verts) {
    for (const Edge &e : get_in_edges(v)) {
      if (!verts.contains(source(e))) {
        ins.insert(e);
      }
    }
    for (const Edge &e : get_all_out_edges(v)) {
      if (!verts.contains(target(e))) {
        if (get_edgetype(e) == EdgeType::Boolean) {
          b_future.insert(e);
        } else {
          outs.insert(e);
          out_lookup.insert({unitmap.at(e), e});
        }
      }
    }
  }
  // Convert to vectors, sort by UnitID, making UnitIDs of corresponding ins and
  // outs match position
  EdgeVec ins_vec{ins.begin(), ins.end()};
  auto compare = [&unitmap](const Edge &e0, const Edge &e1) {
    return unitmap.at(e0) < unitmap.at(e1);
  };
  std::sort(ins_vec.begin(), ins_vec.end(), compare);
  std::vector<std::optional<Edge>> outs_vec;
  for (const Edge &e : ins_vec) {
    if (get_edgetype(e) == EdgeType::Boolean) {
      outs_vec.push_back(std::nullopt);
    } else {
      outs_vec.push_back(out_lookup.at(unitmap.at(e)));
    }
  }
  // b_future can go in any order
  EdgeVec b_future_vec{b_future.begin(), b_future.end()};
  return {ins_vec, outs_vec, b_future_vec, verts};
}

}  // namespace tket
