#include "QInteraction.hpp"

#include <algorithm>
#include <array>
#include <functional>
#include <iostream>
#include <iterator>
#include <memory>
#include <numeric>
#include <set>
#include <string>

#include "Circuit/CircUtils.hpp"
#include "Circuit/Circuit.hpp"
#include "Circuit/DAGDefs.hpp"
#include "OpType/EdgeType.hpp"
#include "OpType/OpType.hpp"
#include "Utils/Assert.hpp"
#include "Utils/GraphHeaders.hpp"

namespace tket {

QInteraction::QInteraction(const Circuit &circ, const Edge &e) : circ_(circ) {
  in_edges_.push_back(e);
  out_edges_.push_back(e);
  n_wires_ = 1;
}

// Combine with another subcircuit disjoint from this. Disjointness is assumed.
void QInteraction::combine(const QInteraction &other) {
  in_edges_.insert(
      in_edges_.end(), other.in_edges_.begin(), other.in_edges_.end());
  out_edges_.insert(
      out_edges_.end(), other.out_edges_.begin(), other.out_edges_.end());
  n_wires_ += other.n_wires_;
  vertices_.insert(other.vertices_.begin(), other.vertices_.end());
}

EdgeVec QInteraction::out_edges() const { return out_edges_; }

VertexSet QInteraction::vertices() const { return vertices_; }

unsigned QInteraction::n_wires() const { return n_wires_; }

unsigned QInteraction::n_vertices() const { return vertices_.size(); }

Subcircuit QInteraction::subcircuit() const {
  return {in_edges_, out_edges_, vertices_};
}

// Append a vertex following the subcircuit. It is assumed that every input
// edge of the vertex matches exactly one output edge of the existing
// subcircuit.
void QInteraction::append(const Vertex &v) {
  EdgeVec v_ins = circ_.get_in_edges_of_type(v, EdgeType::Quantum);
  EdgeVec v_outs = circ_.get_out_edges_of_type(v, EdgeType::Quantum);
  unsigned n = v_ins.size();
  TKET_ASSERT(n == v_outs.size());
  TKET_ASSERT(n <= n_wires_);
  for (unsigned i = 0; i < n; i++) {
    Edge e_new_0 = v_ins[i];
    Edge e_new_1 = v_outs[i];
    bool matched = false;
    for (unsigned j = 0; j < n_wires_; j++) {
      if (out_edges_[j] == e_new_0) {
        TKET_ASSERT(!matched);
        out_edges_[j] = e_new_1;
        matched = true;
      }
    }
    TKET_ASSERT(matched);
  }
  vertices_.insert(v);
}

// Construct an empty system.
QISystem::QISystem(
    Circuit &circ, std::function<Circuit(Circuit)> replacement_func)
    : circ_(circ),
      bin_(),
      interactions_(),
      idx_(0),
      replacement_func_(replacement_func) {}

// Add a new interaction to the system consisting of a single edge, and
// return its index.
int QISystem::create_new_interaction_from_edge(const Edge &e) {
  TKET_ASSERT(!interactions_.contains(idx_));
  interactions_[idx_] = std::make_unique<QInteraction>(circ_, e);
  return idx_++;
}

// Return the set of (indices of) interactions in the system that have v as a
// direct successor.
std::vector<int> QISystem::interactions_feeding_vertex(const Vertex &v) const {
  EdgeVec edges = circ_.get_in_edges_of_type(v, EdgeType::Quantum);
  std::vector<int> meets;
  for (const auto &[i, I] : interactions_) {
    EdgeVec I_outs = I->out_edges();
    EdgeSet I_outset = {I_outs.begin(), I_outs.end()};
    if (std::any_of(edges.begin(), edges.end(), [&I_outset](const Edge &e) {
          return I_outset.contains(e);
        })) {
      meets.push_back(i);
    }
  }
  return meets;
}

// The total width (number of wires) of a subset of the interactions.
unsigned QISystem::total_n_wires(const std::vector<int> &S) const {
  unsigned total = 0;
  for (int i : S) {
    total += interactions_.at(i)->n_wires();
  }
  return total;
}

// From a set of indices, choose the one indexing the largest interaction,
// in terms of vertex count.
int QISystem::largest_interaction(const std::vector<int> &S) const {
  return *std::max_element(
      S.begin(), S.end(), [this](const int &i0, const int &i1) {
        return interactions_.at(i0)->n_vertices() <
               interactions_.at(i1)->n_vertices();
      });
}

// Combine a set of existing interactions into one and append the vertex v.
// It is assumed that the interactions are combinable and the vertex
// appendable.
void QISystem::combine_and_append(const std::vector<int> &S, const Vertex &v) {
  unsigned N = S.size();
  TKET_ASSERT(N > 0);
  iptr &I = interactions_.at(S[0]);
  for (unsigned i = 1; i < N; i++) {
    I->combine(*interactions_.at(S[i]));
    interactions_.erase(S[i]);
  }
  I->append(v);
}

// Close an interaction, squashing it if possible, and erasing it from the
// set. Return true iff any substitution was made, and the (possibly new)
// vector of outgoing edges from the region of the interaction.
std::pair<bool, EdgeVec> QISystem::close_interaction(
    int i, const bool &replace) {
  iptr &I = interactions_.at(i);
  bool changed = false;
  EdgeVec outs = I->out_edges();
  if (I->n_wires() > 1 && replace) {
    Subcircuit sub = I->subcircuit();
    Circuit subc = circ_.subcircuit(sub);
    Circuit replacement = replacement_func_(subc);
    if (replacement.count_gates(OpType::CX) < subc.count_gates(OpType::CX)) {
      // 1. Collect data needed later to reconstruct the list of out-edges:
      std::vector<std::pair<Vertex, port_t>> out_vertex_ports;
      for (const Edge &e : outs) {
        out_vertex_ports.push_back({circ_.target(e), circ_.get_target_port(e)});
      }
      // 2. Do the substitution:
      VertexSet old_vertices = I->vertices();
      bin_.insert(bin_.end(), old_vertices.begin(), old_vertices.end());
      circ_.substitute(replacement, sub, Circuit::VertexDeletion::No);
      // 3. Construct the new list of out-edges:
      EdgeVec new_outs;
      for (const auto &[v, p] : out_vertex_ports) {
        new_outs.push_back(circ_.get_nth_in_edge(v, p));
      }
      changed = true;
      outs = std::move(new_outs);
    }
  }

  interactions_.erase(i);
  return {changed, outs};
}

// Close an interaction and spawn new ones on its outgoing edges. Return true
// iff any substitution is made.
bool QISystem::close_interaction_and_spawn(int i, const bool &replace) {
  auto [changed, outs] = close_interaction(i, replace);
  for (const Edge &e : outs) {
    create_new_interaction_from_edge(e);
  }
  return changed;
}

// Close all interactions that have v as a direct successor, and start new
// ones following them. Return true iff any substitution is made.
bool QISystem::close_interactions_feeding_vertex(
    const Vertex &v, const bool &replace) {
  bool changed = false;
  std::vector<int> v_interactions = interactions_feeding_vertex(v);

  for (int i : v_interactions) {
    auto [change, outs] = close_interaction(i, replace);
    changed |= change;
    for (const Edge &e : outs) {
      if (circ_.target(e) != v) {
        create_new_interaction_from_edge(e);
      }
    }
  }

  for (const Edge &e : circ_.get_out_edges_of_type(v, EdgeType::Quantum)) {
    create_new_interaction_from_edge(e);
  }

  return changed;
}

// Close all interactions. Return true iff any substitution is made.
bool QISystem::close_all_interactions(const bool &replace) {
  bool changed = false;
  // Form set of keys.
  std::set<int> indices;
  for (const auto &pair : interactions_) {
    indices.insert(pair.first);
  }
  // Close each one.
  for (int i : indices) {
    auto [change, outs] = close_interaction(i, replace);
    changed |= change;
  }
  return changed;
}

// Delete all vertices marked for deletion.
void QISystem::destroy_bin() {
  circ_.remove_vertices(
      bin_, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
}

std::map<int, iptr> *QISystem::get_interactions() { return &interactions_; }
}  // namespace tket