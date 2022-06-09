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

#include "ThreeQubitSquash.hpp"

#include <algorithm>
#include <array>
#include <iostream>
#include <iterator>
#include <memory>
#include <numeric>
#include <set>
#include <string>

#include "Circuit/CircUtils.hpp"
#include "Circuit/Circuit.hpp"
#include "Circuit/DAGDefs.hpp"
#include "Circuit/ThreeQubitConversion.hpp"
#include "Decomposition.hpp"
#include "OpType/EdgeType.hpp"
#include "OpType/OpType.hpp"
#include "OptimisationPass.hpp"
#include "Transform.hpp"
#include "Utils/Assert.hpp"
#include "Utils/GraphHeaders.hpp"

namespace tket {

namespace Transforms {

// Helper class defining a pure-quantum subcircuit of up to 3 qubits.
class QInteraction {
 public:
  // Construct an empty subcircuit with a single edge that is both input and
  // output, and no vertices.
  QInteraction(const Circuit &circ, const Edge &e) : circ_(circ) {
    in_edges_.push_back(e);
    out_edges_.push_back(e);
    n_wires_ = 1;
  }

  // Combine with another subcircuit disjoint from this. Disjointness, and the
  // fact that the combined subcircuit has at most three wires, is assumed.
  void combine(const QInteraction &other) {
    in_edges_.insert(
        in_edges_.end(), other.in_edges_.begin(), other.in_edges_.end());
    out_edges_.insert(
        out_edges_.end(), other.out_edges_.begin(), other.out_edges_.end());
    n_wires_ += other.n_wires_;
    vertices_.insert(other.vertices_.begin(), other.vertices_.end());
  }

  EdgeVec out_edges() const { return out_edges_; }

  VertexSet vertices() const { return vertices_; }

  unsigned n_wires() const { return n_wires_; }

  unsigned n_vertices() const { return vertices_.size(); }

  Subcircuit subcircuit() const { return {in_edges_, out_edges_, vertices_}; }

  // Append a vertex following the subcircuit. It is assumed that every input
  // edge of the vertex matches exactly one output edge of the existing
  // subcircuit.
  void append(const Vertex &v) {
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

 private:
  const Circuit &circ_;
  EdgeVec in_edges_;
  EdgeVec out_edges_;
  unsigned n_wires_;    // number of in/out edge pairs
  VertexSet vertices_;  // all internal vertices
};

typedef std::unique_ptr<QInteraction> iptr;

// Candidate substitution for a 2- or 3-qubit circuit.
static Circuit candidate_sub(const Circuit &circ) {
  unsigned n_qb = circ.n_qubits();
  if (n_qb == 2) {
    Circuit repl = two_qubit_canonical(get_matrix_from_2qb_circ(circ));
    clifford_simp(false).apply(repl);
    return repl;
  } else {
    TKET_ASSERT(n_qb == 3);
    Circuit repl = three_qubit_synthesis(get_3q_unitary(circ));
    // TODO: for now we decompose all the way to CX. In the future, it's worth
    // considering keeping TK2, and decomposing to CX (or other gates) later
    // when necessary.
    decompose_TK2().apply(repl);
    clifford_simp(false).apply(repl);
    return repl;
  }
}

// Helper class representing a system of disjoint interactions, each with at
// most three qubits. The interactions are represented by integer labels.
class QISystem {
 public:
  // Construct an empty system.
  explicit QISystem(Circuit &circ)
      : circ_(circ), bin_(), interactions_(), idx_(0) {}

  // Add a new interaction to the system consisting of a single edge, and
  // return its index.
  int create_new_interaction_from_edge(const Edge &e) {
    TKET_ASSERT(!interactions_.contains(idx_));
    interactions_[idx_] = std::make_unique<QInteraction>(circ_, e);
    return idx_++;
  }

  // Return the set of (indices of) interactions in the system that have v as a
  // direct successor.
  std::vector<int> interactions_feeding_vertex(const Vertex &v) const {
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
  unsigned total_n_wires(const std::vector<int> &S) const {
    unsigned total = 0;
    for (int i : S) {
      total += interactions_.at(i)->n_wires();
    }
    return total;
  }

  // From a set of indices, choose the one indexing the largest interaction,
  // in terms of vertex count.
  int largest_interaction(const std::vector<int> &S) const {
    return *std::max_element(
        S.begin(), S.end(), [this](const int &i0, const int &i1) {
          return interactions_.at(i0)->n_vertices() <
                 interactions_.at(i1)->n_vertices();
        });
  }

  // Combine a set of existing interactions into one and append the vertex v.
  // It is assumed that the interactions are combinable and the vertex
  // appendable.
  void combine_and_append(const std::vector<int> &S, const Vertex &v) {
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
  std::pair<bool, EdgeVec> close_interaction(int i) {
    iptr &I = interactions_.at(i);
    bool changed = false;
    EdgeVec outs = I->out_edges();
    switch (I->n_wires()) {
      case 1:
        break;
      case 2:
      case 3: {
        Subcircuit sub = I->subcircuit();
        Circuit subc = circ_.subcircuit(sub);
        Circuit replacement = candidate_sub(subc);
        if (replacement.count_gates(OpType::CX) <
            subc.count_gates(OpType::CX)) {
          // 1. Collect data needed later to reconstruct the list of out-edges:
          std::vector<std::pair<Vertex, port_t>> out_vertex_ports;
          for (const Edge &e : outs) {
            out_vertex_ports.push_back(
                {circ_.target(e), circ_.get_target_port(e)});
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
        break;
      }
      default:
        TKET_ASSERT(!"Interaction with invalid number of wires");
    }
    interactions_.erase(i);
    return {changed, outs};
  }

  // Close an interaction and spawn new ones on its outgoing edges. Return true
  // iff any substitution is made.
  bool close_interaction_and_spawn(int i) {
    auto [changed, outs] = close_interaction(i);
    for (const Edge &e : outs) {
      create_new_interaction_from_edge(e);
    }
    return changed;
  }

  // Close all interactions that have v as a direct successor, and start new
  // ones following them. Return true iff any substitution is made.
  bool close_interactions_feeding_vertex(const Vertex &v) {
    bool changed = false;
    std::vector<int> v_interactions = interactions_feeding_vertex(v);

    for (int i : v_interactions) {
      auto [change, outs] = close_interaction(i);
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
  bool close_all_interactions() {
    bool changed = false;
    // Form set of keys.
    std::set<int> indices;
    for (const auto &pair : interactions_) {
      indices.insert(pair.first);
    }
    // Close each one.
    for (int i : indices) {
      auto [change, outs] = close_interaction(i);
      changed |= change;
    }
    return changed;
  }

  // Delete all vertices marked for deletion.
  void destroy_bin() {
    circ_.remove_vertices(
        bin_, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
  }

 private:
  Circuit &circ_;
  VertexList bin_;
  std::map<int, iptr> interactions_;
  int idx_;
};

Transform three_qubit_squash() {
  return Transform([](Circuit &circ) {
    bool changed = false;

    // Step through the vertices in topological order.
    QISystem Is(circ);  // set of "live" interactions
    for (const Vertex &v : circ.vertices_in_order()) {
      const EdgeVec v_q_ins = circ.get_in_edges_of_type(v, EdgeType::Quantum);
      const EdgeVec v_q_outs = circ.get_out_edges_of_type(v, EdgeType::Quantum);
      unsigned n_q_ins = v_q_ins.size();
      unsigned n_q_outs = v_q_outs.size();

      // If v has no quantum wires, ignore it and move on.
      if (n_q_ins == 0 && n_q_outs == 0) continue;

      // If v is initial, create an interaction from its out-edge, and move on.
      if (n_q_ins == 0) {
        TKET_ASSERT(n_q_outs == 1);
        Is.create_new_interaction_from_edge(v_q_outs[0]);
        continue;
      }

      // If v is final, ignore it and move on.
      if (n_q_outs == 0) continue;

      // It's an internal operation with >0 quantum wires.
      TKET_ASSERT(n_q_ins == n_q_outs);

      Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
      OpType optype = op->get_type();

      // If there are any incoming classical wires, or if this is a Barrier or
      // Reset or Collapse operation, or if the operation contains symbols,
      // close all existing interactions meeting v and create new ones, then
      // move on.
      if (!circ.get_in_edges_of_type(v, EdgeType::Classical).empty() ||
          !circ.get_in_edges_of_type(v, EdgeType::Boolean).empty() ||
          optype == OpType::Barrier || optype == OpType::Reset ||
          optype == OpType::Collapse || !op->free_symbols().empty()) {
        changed |= Is.close_interactions_feeding_vertex(v);
        continue;
      }

      // The circuit should contain only 1-qubit and CX gates.
      if ((n_q_ins == 2 && optype != OpType::CX) || (n_q_ins > 2)) {
        throw std::invalid_argument(
            "Three-qubit squash requires circuits with 1q and CX gates only");
      }

      // Absorb v into existing interactions, closing or merging as necessary.
      bool done_with_v = false;
      while (!done_with_v) {
        std::vector<int> v_Is = Is.interactions_feeding_vertex(v);
        unsigned total_n_wires = Is.total_n_wires(v_Is);
        if (total_n_wires <= 3) {
          Is.combine_and_append(v_Is, v);
          done_with_v = true;
        } else {
          // Close one of the interactions meeting v.
          int i = Is.largest_interaction(v_Is);
          changed |= Is.close_interaction_and_spawn(i);
        }
      }
    }

    // Close all remaining interactions.
    changed |= Is.close_all_interactions();

    // Delete removed vertices.
    Is.destroy_bin();

    return changed;
  });
}

}  // namespace Transforms

}  // namespace tket
