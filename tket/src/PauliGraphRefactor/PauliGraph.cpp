// Copyright 2019-2024 Cambridge Quantum Computing
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

#include "tket/PauliGraphRefactor/PauliGraph.hpp"

#include "tket/Utils/SequencedContainers.hpp"

namespace tket {
namespace pg {

/**
 * PauliGraph Implementation
 */

PauliGraph::PauliGraph()
    : pauli_ac_(0, 0),
      pauli_index_(),
      c_graph_(),
      qubits_(),
      bits_(),
      last_writes_(),
      last_reads_(),
      initial_tableau_(std::nullopt),
      final_tableau_(std::nullopt) {}

PauliGraph::PauliGraph(const std::set<Qubit>& qubits, const std::set<Bit>& bits)
    : pauli_ac_(0, 0),
      pauli_index_(),
      c_graph_(),
      qubits_(qubits),
      bits_(bits),
      last_writes_(),
      last_reads_(),
      initial_tableau_(std::nullopt),
      final_tableau_(std::nullopt) {}

const std::set<Qubit>& PauliGraph::get_qubits() const { return qubits_; }

const std::set<Bit>& PauliGraph::get_bits() const { return bits_; }

std::optional<PGVert> PauliGraph::get_input_tableau() const {
  return initial_tableau_;
}

std::optional<PGVert> PauliGraph::get_output_tableau() const {
  return final_tableau_;
}

PGOp_ptr PauliGraph::get_vertex_PGOp_ptr(const PGVert& v) const {
  return c_graph_[v];
}

void PauliGraph::to_graphviz(std::ostream& out) const {
  out << "digraph G {\ncompound = true;\n";

  std::map<PGVert, unsigned> v_map;
  unsigned i = 0;
  BGL_FORALL_VERTICES(v, c_graph_, PGClassicalGraph) {
    v_map.insert({v, i});
    PGOp_ptr op = c_graph_[v];
    out << "subgraph cluster" << i << "{\nlabel = \"" << op->get_name()
        << "\";\n";
    if (op->n_paulis() == 0) {
      out << "p" << i << "_0;\n";
    } else {
      for (unsigned j = 0; j < op->n_paulis(); ++j) {
        out << "p" << i << "_" << j << ";\n";
      }
    }
    out << "}\n";
    ++i;
  }

  BGL_FORALL_EDGES(e, c_graph_, PGClassicalGraph) {
    PGVert vs = boost::source(e, c_graph_);
    PGVert vt = boost::target(e, c_graph_);
    unsigned vsi = v_map.at(vs);
    unsigned vti = v_map.at(vt);
    out << "p" << vsi << "_0 -> p" << vti << "_0 [ltail=cluster" << vsi
        << ",lhead=cluster" << vti << "];\n";
  }

  for (const PGPauli& r_pauli : pauli_index_.get<TagID>()) {
    for (const PGPauli& c_pauli : pauli_index_.get<TagID>()) {
      if (pauli_ac_(r_pauli.index, c_pauli.index)) {
        out << "p" << v_map.at(c_pauli.vert) << "_" << c_pauli.port << " -> p"
            << v_map.at(r_pauli.vert) << "_" << r_pauli.port << ";\n";
      }
    }
  }

  out << "}";
}

PGVert PauliGraph::add_vertex_at_end(PGOp_ptr op) {
  // Add vertex
  PGVert v = boost::add_vertex(c_graph_);
  c_graph_[v] = op;
  unsigned mat_offset = pauli_ac_.rows();
  pauli_ac_.conservativeResize(
      mat_offset + op->n_paulis(), mat_offset + op->n_paulis());
  for (unsigned i = 0; i < op->n_paulis(); ++i) {
    pauli_index_.insert({mat_offset + i, v, i});
    pauli_ac_.col(mat_offset + i).setZero();
    for (const std::pair<const Qubit, Pauli>& qp : op->port(i).string)
      qubits_.insert(qp.first);
  }
  for (const Bit& b : op->read_bits()) bits_.insert(b);
  for (const Bit& b : op->write_bits()) bits_.insert(b);
  if (op->get_type() == PGOpType::InputTableau) {
    if (boost::num_vertices(c_graph_) != 1)
      throw PGError(
          "Cannot insert InputTableau into PauliGraph - other operations "
          "already exist");
    initial_tableau_ = v;
  } else if (final_tableau_) {
    throw PGError(
        "Cannot insert additional operations to the end of the PauliGraph "
        "after the final tableau");
  } else if (op->get_type() == PGOpType::OutputTableau) {
    final_tableau_ = v;
  }
  // Find ancestors in the anticommutation matrix
  std::vector<SpPauliStabiliser> active = op->active_paulis();
  for (unsigned i = 0; i < active.size(); ++i) {
    for (const PGPauli& prev_pauli : pauli_index_.get<TagID>()) {
      if (prev_pauli.vert == v) continue;
      PGOp_ptr other_op = c_graph_[prev_pauli.vert];
      SpPauliStabiliser other_pauli = other_op->port(prev_pauli.port);
      pauli_ac_(mat_offset + i, prev_pauli.index) =
          !active.at(i).commutes_with(other_pauli);
    }
  }
  // Find classical predecessors
  std::unordered_set<PGVert> c_preds;
  for (const Bit& b : op->read_bits()) {
    auto prev = last_writes_.find(b);
    if (prev != last_writes_.end()) c_preds.insert(prev->second);
    last_reads_[b].insert(v);
  }
  for (const Bit& b : op->write_bits()) {
    auto prev = last_writes_.find(b);
    if (prev != last_writes_.end()) c_preds.insert(prev->second);
    last_writes_[b] = v;
    auto reads = last_reads_.find(b);
    if (reads != last_reads_.end()) {
      for (const PGVert& read : reads->second) c_preds.insert(read);
      last_reads_.erase(reads);
    }
  }
  // If op both reads and writes to the same bit, it may have been inserted into
  // c_preds
  auto self_loop = c_preds.find(v);
  if (self_loop != c_preds.end()) c_preds.erase(self_loop);
  for (const PGVert& pred : c_preds) boost::add_edge(pred, v, c_graph_);
  return v;
}

void PauliGraph::multiply_strings(
    unsigned source_r, unsigned target_r, quarter_turns_t coeff) {
  PGPauli source_pgp = *pauli_index_.get<TagID>().find(source_r);
  PGOp_ptr source_op = c_graph_[source_pgp.vert];
  PGPauli target_pgp = *pauli_index_.get<TagID>().find(target_r);
  PGOp_ptr target_op = c_graph_[target_pgp.vert];
  // Update strings in PGOps
  SpPauliStabiliser& target_port = target_op->port(target_pgp.port);
  target_port = source_op->port(source_pgp.port) * target_port;
  target_port.coeff = (target_port.coeff + coeff) % 4;
  // Update anticommutation matrix
  for (unsigned i = 0; i < pauli_ac_.rows(); ++i) {
    pauli_ac_(i, target_pgp.index) =
        pauli_ac_(i, source_pgp.index) ^ pauli_ac_(i, target_pgp.index);
    pauli_ac_(target_pgp.index) =
        pauli_ac_(source_pgp.index, i) ^ pauli_ac_(target_pgp.index, i);
  }
  pauli_ac_(target_pgp.index, target_pgp.index) = false;
}

void PauliGraph::verify() const {
  // Check validity of the graphs by finding a mutual topological ordering
  std::unordered_set<PGVert> consumed;
  std::map<Bit, PGVert> previous_write;
  std::map<Bit, std::unordered_set<PGVert>> previous_reads;
  bool found_more = true;
  while (found_more) {
    found_more = false;
    BGL_FORALL_VERTICES(v, c_graph_, PGClassicalGraph) {
      if (consumed.find(v) != consumed.end()) continue;
      bool initial = true;
      auto in_edge_range = boost::in_edges(v, c_graph_);
      for (auto it = in_edge_range.first; it != in_edge_range.second; ++it) {
        if (consumed.find(boost::source(*it, c_graph_)) == consumed.end()) {
          initial = false;
          break;
        }
      }
      if (!initial) continue;
      auto range = pauli_index_.get<TagOp>().equal_range(v);
      for (auto it = range.first; it != range.second; ++it) {
        for (const PGPauli& c_pauli : pauli_index_.get<TagID>()) {
          if (pauli_ac_(it->index, c_pauli.index) &&
              (consumed.find(c_pauli.vert) == consumed.end())) {
            initial = false;
            break;
          }
        }
      }
      if (!initial) continue;
      // Found a valid next vertex, so check its relative validity
      found_more = true;
      consumed.insert(v);
      // Check Classical history contains exactly predecessor hazards for all
      // classical bits and all active bits are registered
      PGOp_ptr op = c_graph_[v];
      std::unordered_set<PGVert> justified_preds;
      for (const Bit& b : op->read_bits()) {
        if (bits_.find(b) == bits_.end())
          throw PGError("PGOp reads from unregistered bit: " + op->get_name());
        auto w_found = previous_write.find(b);
        if (w_found != previous_write.end()) {
          if (!boost::edge(w_found->second, v, c_graph_).second)
            throw PGError(
                "No edge in PGClassicalGraph for RAW dependency on bit " +
                b.repr() + " between " + c_graph_[w_found->second]->get_name() +
                " and " + op->get_name());
          justified_preds.insert(w_found->second);
        }
        previous_reads[b].insert(v);
      }
      for (const Bit& b : op->write_bits()) {
        if (bits_.find(b) == bits_.end())
          throw PGError("PGOp writes to unregistered bit: " + op->get_name());
        auto w_found = previous_write.find(b);
        if (w_found != previous_write.end()) {
          if (!boost::edge(w_found->second, v, c_graph_).second)
            throw PGError(
                "No edge in PGClassicalGraph for WAW dependency on bit " +
                b.repr() + " between " + c_graph_[w_found->second]->get_name() +
                " and " + op->get_name());
          justified_preds.insert(w_found->second);
        }
        auto r_found = previous_reads.find(b);
        if (r_found != previous_reads.end()) {
          for (const PGVert& u : r_found->second) {
            if (u == v) continue;
            if (!boost::edge(u, v, c_graph_).second)
              throw PGError(
                  "No edge in PGClassicalGraph for WAR dependency on bit " +
                  b.repr() + " between " + c_graph_[u]->get_name() + " and " +
                  op->get_name());
            justified_preds.insert(u);
          }
          r_found->second.clear();
        }
        previous_write[b] = v;
      }
      for (auto it = in_edge_range.first; it != in_edge_range.second; ++it) {
        PGVert s = boost::source(*it, c_graph_);
        if (justified_preds.find(s) == justified_preds.end())
          throw PGError(
              "Edge in PGClassicalGraph despite no dependency between " +
              c_graph_[s]->get_name() + " and " + op->get_name());
      }
      // Check Pauli history contains all anti-commuting terms and all active
      // qubits are registered
      std::vector<SpPauliStabiliser> paulis = op->active_paulis();
      for (auto it = range.first; it != range.second; ++it) {
        const SpPauliStabiliser& tensor = paulis.at(it->port);
        for (const std::pair<const Qubit, Pauli>& qp : tensor.string) {
          if (qubits_.find(qp.first) == qubits_.end())
            throw PGError(
                "PGOp interacts with unregistered qubit: " + op->get_name());
        }
        for (const PGPauli& c_pauli : pauli_index_.get<TagID>()) {
          if (consumed.find(c_pauli.vert) == consumed.end()) continue;
          if ((c_pauli.vert != v) &&
              !tensor.commutes_with(
                  c_graph_[c_pauli.vert]->port(c_pauli.port))) {
            if (!pauli_ac_(it->index, c_pauli.index))
              throw PGError(
                  "PauliGraph anticommutation matrix is missing a link "
                  "between " +
                  c_graph_[c_pauli.vert]->get_name() + " and " +
                  op->get_name());
          } else if (pauli_ac_(it->index, c_pauli.index))
            throw PGError(
                "PauliGraph anticommutation matrix contains an invalid link "
                "between " +
                c_graph_[c_pauli.vert]->get_name() + " and " + op->get_name());
        }
      }

      // Tableau conditions
      if (op->get_type() == PGOpType::InputTableau) {
        if (initial_tableau_ != v)
          throw PGError("PauliGraph contains an untracked InputTableau");
        // No predecessors
        for (auto it = range.first; it != range.second; ++it) {
          for (unsigned c = 0; c < pauli_ac_.cols(); ++c) {
            if (pauli_ac_(it->index, c))
              throw PGError(
                  "PauliGraph InputTableau has predecessors in anticommutation "
                  "matrix");
          }
        }
        // Commuting rows
        const PGInputTableau& tab = dynamic_cast<const PGInputTableau&>(*op);
        for (unsigned i = 0; i < tab.n_paulis(); ++i) {
          for (unsigned j = i + 1; j < tab.n_paulis(); ++j) {
            if (tab.get_full_row(i).first.commutes_with(
                    tab.get_full_row(j).first) !=
                tab.get_full_row(i).second.commutes_with(
                    tab.get_full_row(j).second))
              throw PGError(
                  "InputTableau of PauliGraph does not have commuting rows");
          }
        }
      } else if (op->get_type() == PGOpType::OutputTableau) {
        if (final_tableau_ != v)
          throw PGError("PauliGraph contains an untracked OutputTableau");
        // No successors
        for (auto it = range.first; it != range.second; ++it) {
          for (unsigned r = 0; r < pauli_ac_.rows(); ++r) {
            if (pauli_ac_(r, it->index))
              throw PGError(
                  "PauliGraph OutputTableau has successors in anticommutation "
                  "matrix");
          }
        }
        // Commuting rows
        const PGOutputTableau& tab = dynamic_cast<const PGOutputTableau&>(*op);
        for (unsigned i = 0; i < tab.n_paulis(); ++i) {
          for (unsigned j = i + 1; j < tab.n_paulis(); ++j) {
            if (tab.get_full_row(i).first.commutes_with(
                    tab.get_full_row(j).first) !=
                tab.get_full_row(i).second.commutes_with(
                    tab.get_full_row(j).second))
              throw PGError(
                  "OutputTableau of PauliGraph does not have commuting rows");
          }
        }
      }
    }
  }
  if (consumed.size() != boost::num_vertices(c_graph_))
    throw PGError("Cannot obtain a topological ordering of PauliGraph");
}

std::list<PGOp_ptr> PauliGraph::pgop_sequence() const {
  std::list<PGOp_ptr> sequence;
  std::list<std::list<PGOp_ptr>> set_list = pgop_commuting_sets();
  for (const std::list<PGOp_ptr>& set : set_list)
    sequence.insert(sequence.end(), set.begin(), set.end());
  return sequence;
}

std::list<std::list<PGOp_ptr>> PauliGraph::pgop_commuting_sets() const {
  std::list<std::list<PGOp_ptr>> set_list;
  sequence_set_t<PGVert> remaining;
  BGL_FORALL_VERTICES(v, c_graph_, PGClassicalGraph) { remaining.insert(v); }
  while (!remaining.empty()) {
    std::list<PGVert> initials;
    for (const PGVert& v : remaining.get<TagSeq>()) {
      bool initial = true;
      auto in_edge_range = boost::in_edges(v, c_graph_);
      for (auto it = in_edge_range.first; it != in_edge_range.second; ++it) {
        if (remaining.find(boost::source(*it, c_graph_)) != remaining.end()) {
          initial = false;
          break;
        }
      }
      if (!initial) continue;
      auto range = pauli_index_.get<TagOp>().equal_range(v);
      for (auto it = range.first; it != range.second; ++it) {
        for (const PGPauli& c_pauli : pauli_index_.get<pg::TagID>()) {
          if (pauli_ac_(it->index, c_pauli.index) &&
              (remaining.find(c_pauli.vert) != remaining.end())) {
            initial = false;
            break;
          }
        }
      }
      if (initial) initials.push_back(v);
    }
    std::list<PGOp_ptr> initial_ops;
    for (const PGVert& v : initials) initial_ops.push_back(c_graph_[v]);
    set_list.push_back(initial_ops);
    auto& lookup = remaining.get<TagKey>();
    for (const PGVert& v : initials) {
      lookup.erase(lookup.find(v));
    }
  }
  return set_list;
}

void PauliGraph::symbol_substitution(
    const SymEngine::map_basic_basic& sub_map) {
  BGL_FORALL_VERTICES(v, c_graph_, PGClassicalGraph) {
    PGOp_ptr new_op = get_vertex_PGOp_ptr(v)->symbol_substitution(sub_map);
    if (new_op) c_graph_[v] = new_op;
  }
}

SymSet PauliGraph::free_symbols() const {
  SymSet symbols;
  BGL_FORALL_VERTICES(v, c_graph_, PGClassicalGraph) {
    SymSet s = get_vertex_PGOp_ptr(v)->free_symbols();
    symbols.insert(s.begin(), s.end());
  }
  return symbols;
}

bool PauliGraph::is_symbolic() const { return !free_symbols().empty(); }

}  // namespace pg
}  // namespace tket
