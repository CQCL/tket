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

#include "PauliGraph.hpp"

#include "Gate/Gate.hpp"
#include "Utils/GraphHeaders.hpp"

namespace tket {

bool operator<(
    const PauliGadgetProperties &pgp1, const PauliGadgetProperties &pgp2) {
  return (pgp1.tensor_.string < pgp2.tensor_.string);
}

PauliGraph::PauliGraph(unsigned n) : cliff_(n) {}

PauliGraph::PauliGraph(const qubit_vector_t &qbs, const bit_vector_t &bits)
    : cliff_(qbs), bits_(bits) {}

PauliVertSet PauliGraph::get_successors(const PauliVert &vert) const {
  PauliVertSet succs;
  for (auto iter = boost::adjacent_vertices(vert, graph_);
       iter.first != iter.second; iter.first++) {
    succs.insert(*iter.first);
  }
  return succs;
}

PauliVertSet PauliGraph::get_predecessors(const PauliVert &vert) const {
  PauliVertSet preds;
  for (auto iter = boost::inv_adjacent_vertices(vert, graph_);
       iter.first != iter.second; iter.first++) {
    preds.insert(*iter.first);
  }
  return preds;
}

PauliEdgeSet PauliGraph::get_in_edges(const PauliVert &vert) const {
  PauliEdgeSet ins;
  for (auto iter = boost::in_edges(vert, graph_); iter.first != iter.second;
       iter.first++) {
    ins.insert(*iter.first);
  }
  return ins;
}

PauliEdgeSet PauliGraph::get_out_edges(const PauliVert &vert) const {
  PauliEdgeSet outs;
  for (auto iter = boost::out_edges(vert, graph_); iter.first != iter.second;
       iter.first++) {
    outs.insert(*iter.first);
  }
  return outs;
}

PauliVert PauliGraph::source(const PauliEdge &edge) const {
  return boost::source(edge, graph_);
}

PauliVert PauliGraph::target(const PauliEdge &edge) const {
  return boost::target(edge, graph_);
}

void PauliGraph::apply_gate_at_end(
    const Gate &gate, const unit_vector_t &args) {
  for (const UnitID &arg : args) {
    if (arg.type() == UnitType::Qubit) {
      if (measures_.left.find(Qubit(arg)) != measures_.left.end()) {
        throw NotImplemented(
            "PauliGraph does not support mid-circuit measurements "
            "- cannot add gate after measure on qubit " +
            arg.repr());
      }
    } else if (measures_.right.find(Bit(arg)) != measures_.right.end()) {
      throw NotImplemented(
          "PauliGraph does not support mid-circuit measurements - "
          "cannot add gate after measure to bit " +
          arg.repr());
    }
  }
  OpType type = gate.get_type();
  if (type == OpType::Measure) {
    measures_.insert({Qubit(args.at(0)), Bit(args.at(1))});
    return;
  }
  qubit_vector_t qbs = {args.begin(), args.end()};
  switch (type) {
    case OpType::Z:
    case OpType::X:
    case OpType::Y:
    case OpType::S:
    case OpType::Sdg:
    case OpType::V:
    case OpType::Vdg:
    case OpType::H:
    case OpType::CX:
    case OpType::CY:
    case OpType::CZ:
    case OpType::SWAP: {
      cliff_.apply_gate_at_end(type, qbs);
      break;
    }
    case OpType::Rz: {
      QubitPauliTensor pauli = cliff_.get_zpauli(qbs.at(0));
      Expr angle = gate.get_params().at(0);
      std::optional<unsigned> cliff_angle = equiv_Clifford(angle);
      if (cliff_angle) {
        for (unsigned i = 0; i < cliff_angle.value(); i++) {
          cliff_.apply_gate_at_end(OpType::S, qbs);
        }
      } else
        apply_pauli_gadget_at_end(pauli, angle);
      break;
    }
    case OpType::Rx: {
      QubitPauliTensor pauli = cliff_.get_xpauli(qbs.at(0));
      Expr angle = gate.get_params().at(0);
      std::optional<unsigned> cliff_angle = equiv_Clifford(angle);
      if (cliff_angle) {
        for (unsigned i = 0; i < cliff_angle.value(); i++) {
          cliff_.apply_gate_at_end(OpType::V, qbs);
        }
      } else
        apply_pauli_gadget_at_end(pauli, angle);
      break;
    }
    case OpType::Ry: {
      Expr angle = gate.get_params().at(0);
      std::optional<unsigned> cliff_angle = equiv_Clifford(angle);
      if (cliff_angle) {
        if (cliff_angle.value() != 0) {
          cliff_.apply_gate_at_end(OpType::V, qbs);
          for (unsigned i = 0; i < cliff_angle.value(); i++) {
            cliff_.apply_gate_at_end(OpType::S, qbs);
          }
          cliff_.apply_gate_at_end(OpType::Vdg, qbs);
        }
      } else {
        QubitPauliTensor zpauli = cliff_.get_zpauli(qbs.at(0));
        QubitPauliTensor xpauli = cliff_.get_xpauli(qbs.at(0));
        QubitPauliTensor ypauli = i_ * xpauli * zpauli;
        apply_pauli_gadget_at_end(ypauli, angle);
      }
      break;
    }
    case OpType::T: {
      QubitPauliTensor pauli = cliff_.get_zpauli(qbs.at(0));
      apply_pauli_gadget_at_end(pauli, 0.25);
      break;
    }
    case OpType::Tdg: {
      QubitPauliTensor pauli = cliff_.get_zpauli(qbs.at(0));
      apply_pauli_gadget_at_end(pauli, 0.25);
      break;
    }
    case OpType::ZZMax: {
      cliff_.apply_gate_at_end(OpType::H, {qbs.at(1)});
      cliff_.apply_gate_at_end(OpType::CX, qbs);
      cliff_.apply_gate_at_end(OpType::Sdg, {qbs.at(0)});
      cliff_.apply_gate_at_end(OpType::S, {qbs.at(1)});
      cliff_.apply_gate_at_end(OpType::V, {qbs.at(1)});
      break;
    }
    case OpType::PhaseGadget:
    case OpType::ZZPhase: {
      Expr angle = gate.get_params().at(0);
      std::optional<unsigned> cliff_angle = equiv_Clifford(angle);
      if (cliff_angle) {
        if (cliff_angle.value() != 0) {
          for (unsigned i = 1; i < qbs.size(); i++) {
            cliff_.apply_gate_at_end(OpType::CX, {qbs.at(i - 1), qbs.at(i)});
          }
          for (unsigned i = 0; i < cliff_angle.value(); i++) {
            cliff_.apply_gate_at_end(OpType::S, {qbs.back()});
          }
          for (unsigned i = qbs.size() - 1; i > 0; i--) {
            cliff_.apply_gate_at_end(OpType::CX, {qbs.at(i - 1), qbs.at(i)});
          }
        }
      } else {
        QubitPauliTensor pauli;
        for (const Qubit &qb : qbs) {
          pauli = pauli * cliff_.get_zpauli(qb);
        }
        apply_pauli_gadget_at_end(pauli, angle);
      }
      break;
    }
    case OpType::XXPhase: {
      Expr angle = gate.get_params().at(0);
      std::optional<unsigned> cliff_angle = equiv_Clifford(angle);
      if (cliff_angle) {
        if (cliff_angle.value() != 0) {
          cliff_.apply_gate_at_end(OpType::CX, {qbs.at(1), qbs.at(0)});
          for (unsigned i = 0; i < cliff_angle.value(); i++) {
            cliff_.apply_gate_at_end(OpType::V, {qbs.back()});
          }
          cliff_.apply_gate_at_end(OpType::CX, {qbs.at(1), qbs.at(0)});
        }
      } else {
        QubitPauliTensor pauli =
            cliff_.get_xpauli(qbs.at(0)) * cliff_.get_xpauli(qbs.at(1));
        apply_pauli_gadget_at_end(pauli, angle);
      }
      break;
    }
    case OpType::YYPhase: {
      Expr angle = gate.get_params().at(0);
      std::optional<unsigned> cliff_angle = equiv_Clifford(angle);
      if (cliff_angle) {
        if (cliff_angle.value() != 0) {
          const Qubit &arg0 = qbs.at(0);
          const Qubit &arg1 = qbs.at(1);
          cliff_.apply_gate_at_end(OpType::V, {arg0});
          cliff_.apply_gate_at_end(OpType::V, {arg1});
          cliff_.apply_gate_at_end(OpType::CX, {arg1, arg0});
          for (unsigned i = 0; i < cliff_angle.value(); i++) {
            cliff_.apply_gate_at_end(OpType::V, {arg1});
          }
          cliff_.apply_gate_at_end(OpType::CX, {arg1, arg0});
          cliff_.apply_gate_at_end(OpType::Vdg, {arg0});
          cliff_.apply_gate_at_end(OpType::Vdg, {arg1});
        }
      } else {
        QubitPauliTensor pauli =
            -1. * cliff_.get_xpauli(qbs.at(0)) * cliff_.get_zpauli(qbs.at(0)) *
            cliff_.get_xpauli(qbs.at(1)) * cliff_.get_zpauli(qbs.at(1));
        apply_pauli_gadget_at_end(pauli, angle);
      }
      break;
    }
    default: {
      throw NotImplemented("Cannot add gate to PauliGraph: " + gate.get_name());
    }
  }
}

void PauliGraph::apply_pauli_gadget_at_end(
    const QubitPauliTensor &pauli, const Expr &angle) {
  PauliVertSet to_search = end_line_;
  PauliVertSet commuted;
  PauliVert new_vert = boost::add_vertex(graph_);
  graph_[new_vert] = {pauli, angle};
  while (!to_search.empty()) {
    // Get next candidate parent
    PauliVert to_compare = *to_search.begin();
    to_search.erase(to_search.begin());

    // Check that we have already commuted past all of its children
    bool ready = true;
    for (const PauliVert &child : get_successors(to_compare)) {
      if (commuted.get<TagKey>().find(child) == commuted.get<TagKey>().end()) {
        ready = false;
        break;
      }
    }
    if (!ready) continue;

    // Check if we can commute past it
    QubitPauliTensor compare_pauli = graph_[to_compare].tensor_;
    if (pauli.commutes_with(compare_pauli)) {
      if (pauli.string == compare_pauli.string) {
        // Identical strings - we can merge vertices
        if (pauli.coeff == compare_pauli.coeff) {
          graph_[to_compare].angle_ += angle;
        } else {
          graph_[to_compare].angle_ -= angle;
        }
        boost::clear_vertex(new_vert, graph_);
        boost::remove_vertex(new_vert, graph_);

        std::optional<unsigned> cl_ang =
            equiv_Clifford(graph_[to_compare].angle_);
        if (cl_ang) {
          cliff_.apply_pauli_at_front(graph_[to_compare].tensor_, *cl_ang);
          start_line_.erase(to_compare);
          for (const PauliVert &v : get_predecessors(to_compare)) {
            if (boost::out_degree(v, graph_) == 1) {
              end_line_.insert(v);
            }
          }
          end_line_.erase(to_compare);
          boost::clear_vertex(to_compare, graph_);
          boost::remove_vertex(to_compare, graph_);
        }
        return;
      } else {
        // Commute and continue searching
        PauliVertSet preds = get_predecessors(to_compare);
        to_search.insert(preds.begin(), preds.end());
        commuted.insert(to_compare);
      }
    } else {
      // Does not commute - add dependency edge
      boost::add_edge(to_compare, new_vert, graph_);
      end_line_.erase(to_compare);
    }
  }
  end_line_.insert(new_vert);
  if (get_predecessors(new_vert).empty()) start_line_.insert(new_vert);
}

PauliGraph::TopSortIterator::TopSortIterator()
    : pg_(nullptr),
      current_vert_(boost::graph_traits<PauliDAG>::null_vertex()) {}

PauliGraph::TopSortIterator::TopSortIterator(const PauliGraph &pg) {
  if (pg.start_line_.empty()) {
    current_vert_ = boost::graph_traits<PauliDAG>::null_vertex();
    return;
  }
  pg_ = &pg;
  for (const PauliVert &vert : pg.start_line_) {
    search_set_.insert({pg.graph_[vert].tensor_, vert});
  }
  current_vert_ = search_set_.begin()->second;
  search_set_.erase(search_set_.begin());
  visited_ = {current_vert_};
  for (const PauliVert &child : pg_->get_successors(current_vert_)) {
    search_set_.insert({pg.graph_[child].tensor_, child});
  }
}

const PauliVert &PauliGraph::TopSortIterator::operator*() const {
  return current_vert_;
}

const PauliVert *PauliGraph::TopSortIterator::operator->() const {
  return &current_vert_;
}

bool PauliGraph::TopSortIterator::operator==(
    const TopSortIterator &other) const {
  return this->current_vert_ == other.current_vert_;
}

bool PauliGraph::TopSortIterator::operator!=(
    const TopSortIterator &other) const {
  return !(*this == other);
}

PauliGraph::TopSortIterator PauliGraph::TopSortIterator::operator++(int) {
  PauliGraph::TopSortIterator it = *this;
  ++*this;
  return it;
}

PauliGraph::TopSortIterator &PauliGraph::TopSortIterator::operator++() {
  bool found_next = false;
  while (!found_next && !search_set_.empty()) {
    current_vert_ = search_set_.begin()->second;
    search_set_.erase(search_set_.begin());

    // Check that we have visited all parents
    found_next = true;
    for (const PauliVert &parent : pg_->get_predecessors(current_vert_)) {
      if (visited_.find(parent) == visited_.end()) {
        found_next = false;
        break;
      }
    }
  }
  if (found_next) {
    visited_.insert(current_vert_);
    for (const PauliVert &child : pg_->get_successors(current_vert_)) {
      search_set_.insert({pg_->graph_[child].tensor_, child});
    }
  } else {
    *this = TopSortIterator();
  }
  return *this;
}

PauliGraph::TopSortIterator PauliGraph::begin() const {
  return TopSortIterator(*this);
}

PauliGraph::TopSortIterator PauliGraph::end() const {
  return TopSortIterator();
}

void PauliGraph::to_graphviz_file(const std::string &filename) const {
  std::ofstream dot_file(filename);
  to_graphviz(dot_file);
}

void PauliGraph::to_graphviz(std::ostream &out) const {
  out << "digraph G {\n";

  std::map<PauliVert, unsigned> index_map;
  unsigned i = 0;
  BGL_FORALL_VERTICES(v, graph_, PauliDAG) {
    index_map.insert({v, i});
    out << i << " [label = \"" << graph_[v].tensor_.to_str() << ", "
        << graph_[v].angle_ << "\"];\n";
    ++i;
  }

  BGL_FORALL_EDGES(e, graph_, PauliDAG) {
    PauliVert v_so = source(e);
    PauliVert v_ta = target(e);
    out << index_map.at(v_so) << " -> " << index_map.at(v_ta) << ";\n";
  }

  out << "}";
}

}  // namespace tket
