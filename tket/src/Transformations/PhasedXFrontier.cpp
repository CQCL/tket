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

#include "PhasedXFrontier.hpp"

#include <string>

#include "Circuit/CircPool.hpp"
#include "Circuit/CircUtils.hpp"
#include "OpType/OpType.hpp"
#include "StandardSquash.hpp"
#include "Utils/Expression.hpp"

namespace tket {

namespace Transforms {

template <typename T>
static std::map<T, unsigned> count(std::vector<T> xs) {
  std::map<T, unsigned> cnts;
  for (const auto& x : xs) {
    auto it = cnts.find(x);
    if (it != cnts.end()) {
      ++(it->second);
    } else {
      cnts[x] = 1;
    }
  }
  return cnts;
}

bool PhasedXFrontier::is_finished() {
  for (unsigned i = 0; i < circ_.n_qubits(); ++i) {
    Vertex v = circ_.target(intervals_[i].second);
    if (!circ_.detect_final_Op(v)) {
      return false;
    }
  }
  return true;
}

std::set<unsigned> PhasedXFrontier::qubits_ending_in(const Vertex& v) const {
  std::set<unsigned> qubits;
  for (unsigned i = 0; i < circ_.n_qubits(); ++i) {
    if (circ_.target(intervals_[i].second) == v) {
      qubits.insert(i);
    }
  }
  return qubits;
}

void PhasedXFrontier::squash_intervals() {
  for (unsigned i = 0; i < circ_.n_qubits(); ++i) {
    squash_interval(i);
  }
}

OptEdge PhasedXFrontier::get_beta_edge(unsigned i) const {
  const auto& [start, end] = intervals_[i];
  Edge e = start;
  while (e != end) {
    Vertex v = circ_.target(e);
    OpType type = circ_.get_OpType_from_Vertex(v);
    if (type == OpType::PhasedX || type == OpType::NPhasedX) {
      return e;
    }
    e = circ_.get_next_edge(v, e);
  }
  return std::nullopt;
}

OptEdgeVec PhasedXFrontier::get_all_beta_edges() const {
  OptEdgeVec beta_edges;
  std::vector<Vertex> beta_vertices;
  for (unsigned i = 0; i < circ_.n_qubits(); ++i) {
    OptEdge beta_e = get_beta_edge(i);
    beta_edges.push_back(beta_e);
    if (beta_e) {
      beta_vertices.push_back(circ_.target(*beta_e));
    }
  }

  // for OpType::NPhasedX, we need to reset some betas to zero if they
  // are shadowed (i.e. if one NPhasedX is behind another one)
  std::map<Vertex, unsigned> arities = count(beta_vertices);
  for (auto& e : beta_edges) {
    if (e) {
      Vertex v = circ_.target(*e);
      Op_ptr op = circ_.get_Op_ptr_from_Vertex(v);
      if (arities[v] != op->n_qubits()) {
        e = std::nullopt;
      }
    }
  }
  return beta_edges;
}

OptVertexVec PhasedXFrontier::get_all_beta_vertices() const {
  OptVertexVec vertices;
  for (const auto& e : get_all_beta_edges()) {
    if (e) {
      Vertex v = circ_.target(*e);
      vertices.push_back(v);
    } else {
      OptVertex optv;
      vertices.push_back(optv);
    }
  }
  return vertices;
}

std::vector<Expr> PhasedXFrontier::get_all_betas() const {
  std::vector<Expr> betas;
  for (const auto& v : get_all_beta_vertices()) {
    Expr beta = 0.;
    if (v) {
      Op_ptr op = circ_.get_Op_ptr_from_Vertex(*v);
      beta = op->get_params().front();
    }
    betas.push_back(beta);
  }
  return betas;
}

void PhasedXFrontier::next_interval(unsigned i) {
  auto& [start, end] = intervals_[i];

  // new interval
  start = get_interval_start(end);
  end = get_interval_end(start);
}

void PhasedXFrontier::next_multiqb(const Vertex& v) {
  std::set<unsigned> curr_qubits = qubits_ending_in(v);
  // move forward
  for (unsigned i : curr_qubits) {
    next_interval(i);
  }
}

Edge PhasedXFrontier::get_interval_end(Edge e) const {
  Vertex v = circ_.target(e);
  while (!circ_.detect_final_Op(v) && !is_interval_boundary(v)) {
    auto p = circ_.get_next_pair(v, e);
    v = p.first;
    e = p.second;
  }
  return e;
}

Edge PhasedXFrontier::get_interval_start(Edge e) const {
  Vertex v = circ_.target(e);
  if (!circ_.detect_final_Op(v)) {
    e = circ_.get_next_edge(v, e);
  }
  return e;
}

bool PhasedXFrontier::is_interval_boundary(Op_ptr op) {
  OpType type = op->get_type();
  return is_gate_type(type) && as_gate_ptr(op)->n_qubits() > 1 &&
         type != OpType::NPhasedX;
}

bool PhasedXFrontier::is_interval_boundary(Vertex v) const {
  Op_ptr op = circ_.get_Op_ptr_from_Vertex(v);
  return is_interval_boundary(op);
}

/**
 * @brief Implements the AbstractSquasher interface for squashing to PhasedX+Rz
 */
class PhasedXSquasher : public StandardSquasher {
 public:
  PhasedXSquasher()
      : StandardSquasher(
            OpTypeSet{OpType::Rz, OpType::PhasedX},
            CircPool::tk1_to_PhasedXRz) {}

  // Accept any single-qubit gate that has TK1 angles to squash it.
  bool accepts(Gate_ptr gp) const override {
    OpType optype = gp->get_type();
    return is_single_qubit_unitary_type(optype);
  }
};

PhasedXFrontier::PhasedXFrontier(Circuit& circ)
    : intervals_(),
      circ_(circ),
      squasher_(std::make_unique<PhasedXSquasher>(), circ, false) {
  const unsigned n = circ_.n_qubits();
  intervals_.resize(n);

  // initialise intervals
  const qubit_vector_t all_qbs = circ_.all_qubits();
  for (unsigned i = 0; i < n; ++i) {
    Qubit q = all_qbs[i];
    Vertex v_in = circ_.get_in(q);
    EdgeVec e_vec = circ_.get_all_out_edges(v_in);
    TKET_ASSERT(e_vec.size() == 1);

    // interval [start, end]
    Edge start = e_vec.front();
    Edge end = get_interval_end(start);
    intervals_[i] = {start, end};
  }
}

void PhasedXFrontier::squash_interval(unsigned i) {
  auto& [start_e, end_e] = intervals_[i];

  // backup edges
  VertPort start{circ_.source(start_e), circ_.get_source_port(start_e)};
  VertPort end{circ_.target(end_e), circ_.get_target_port(end_e)};

  squasher_.squash_between(start_e, end_e);

  // restore interval edges
  start_e = circ_.get_nth_out_edge(start.first, start.second);
  end_e = circ_.get_nth_in_edge(end.first, end.second);
}

PhasedXFrontier::BackupIntervals PhasedXFrontier::backup_intervals() const {
  BackupIntervals ret;
  for (unsigned i = 0; i < circ_.n_qubits(); ++i) {
    const auto& [start, end] = intervals_[i];
    ret.start.push_back({circ_.source(start), circ_.get_source_port(start)});
    ret.end.push_back({circ_.target(end), circ_.get_target_port(end)});
  }
  return ret;
}

void PhasedXFrontier::restore_intervals(const BackupIntervals& b) {
  for (unsigned i = 0; i < circ_.n_qubits(); ++i) {
    Edge start = circ_.get_nth_out_edge(b.start[i].first, b.start[i].second);
    Edge end = circ_.get_nth_in_edge(b.end[i].first, b.end[i].second);
    intervals_[i].first = start;
    intervals_[i].second = end;
  }
}

void PhasedXFrontier::skip_global_gates(unsigned n) {
  for (unsigned i = 0; i < circ_.n_qubits(); ++i) {
    unsigned count = 0;
    auto& [e, end] = intervals_[i];
    while (e != end) {
      Vertex v = circ_.target(e);
      OpType type = circ_.get_OpType_from_Vertex(v);
      e = circ_.get_next_edge(v, e);
      if (type == OpType::NPhasedX ||
          (circ_.n_qubits() == 1 && type == OpType::PhasedX)) {
        unsigned in_edges = circ_.n_in_edges_of_type(v, EdgeType::Quantum);
        unsigned out_edges = circ_.n_out_edges_of_type(v, EdgeType::Quantum);
        if (in_edges != circ_.n_qubits() || out_edges != circ_.n_qubits()) {
          throw NotValid("Found non-global NPhasedX gate");
        }
        ++count;
        if (count == n) {
          break;
        }
      }
    }
    if (count < n) {
      throw NotValid("Did not find expected global gates");
    }
  }
}

bool PhasedXFrontier::are_phasedx_left() const {
  PhasedXFrontier copy = *this;
  const unsigned n = circ_.n_qubits();
  for (unsigned i = 0; i < n; ++i) {
    copy.next_interval(i);
  }
  OptVertexVec all_phasedx = copy.get_all_beta_vertices();

  return !all_nullopt(all_phasedx);
}

void PhasedXFrontier::insert_1_phasedx(unsigned i) {
  std::vector<Expr> betas = get_all_betas();
  OptEdgeVec edges = get_all_beta_edges();
  OptVertexVec vertices = get_all_beta_vertices();

  if (!vertices[i]) {
    throw NotValid("No PhasedX found on qubit " + std::to_string(i));
  }
  Expr beta = betas[i];

  VertexSet bin;
  EdgeVec in_hole, out_hole;
  Circuit sub1(circ_.n_qubits());
  Circuit sub2(circ_.n_qubits());

  for (unsigned i = 0; i < circ_.n_qubits(); ++i) {
    if (vertices[i]) {
      Vertex v = vertices[i].value();
      Edge e = edges[i].value();
      Op_ptr op = circ_.get_Op_ptr_from_Vertex(v);
      OpType type = op->get_type();
      Expr new_beta = betas[i] - beta;

      in_hole.push_back(e);
      out_hole.push_back(circ_.get_next_edge(v, e));

      if (!bin.contains(v)) {
        bin.insert(v);

        // v's qubits
        std::vector<unsigned> qubits;
        for (unsigned j = 0; j < circ_.n_qubits(); ++j) {
          if (vertices[j] == v) {
            qubits.push_back(j);
          }
        }

        if (type == OpType::NPhasedX || type == OpType::PhasedX) {
          Expr alpha = op->get_params()[1];
          if (!equiv_0(new_beta, 4)) {
            if (equiv_0(new_beta, 2)) {
              if (qubits.size() % 2) {
                sub1.add_phase(-1);
              }
            } else {
              sub2.add_op(type, {new_beta, 0.}, qubits);
            }
          }
          if (!equiv_0(alpha, 2)) {
            for (unsigned q : qubits) {
              sub1.add_op<unsigned>(OpType::Rz, {-alpha}, {q});
              sub2.add_op<unsigned>(OpType::Rz, {alpha}, {q});
            }
          }
        } else {
          throw NotValid("Encountered invalid beta angle OpType");
        }
      }
    } else {
      Edge interval_begin = intervals_[i].first;
      in_hole.push_back(interval_begin);
      out_hole.push_back(interval_begin);
      if (!equiv_0(beta, 2)) {
        sub2.add_op<unsigned>(OpType::PhasedX, {-beta, 0}, {i});
      } else if (!equiv_0(beta, 4)) {
        sub2.add_phase(-1);
      }
    }
  }

  Circuit sub(circ_.n_qubits());
  sub.append(sub1);
  sub.add_op<Qubit>(OpType::NPhasedX, {beta, 0}, sub.all_qubits());
  sub.append(sub2);

  Subcircuit hole{in_hole, out_hole, bin};

  // perform substitution
  BackupIntervals backup = backup_intervals();
  circ_.substitute(sub, hole, Circuit::VertexDeletion::Yes);
  restore_intervals(backup);

  skip_global_gates(1);
}

void PhasedXFrontier::insert_2_phasedx() {
  EdgeVec hole_in, hole_out;
  Circuit sub1(circ_.n_qubits());
  Circuit sub2(circ_.n_qubits());
  Circuit sub3(circ_.n_qubits());
  VertexSet bin;

  std::vector<Expr> betas = get_all_betas();
  OptEdgeVec edges = get_all_beta_edges();
  OptVertexVec vertices = get_all_beta_vertices();

  for (unsigned i = 0; i < circ_.n_qubits(); ++i) {
    if (vertices[i]) {
      Vertex v = vertices[i].value();
      Edge e = edges[i].value();
      Op_ptr op = circ_.get_Op_ptr_from_Vertex(v);
      OpType type = op->get_type();
      hole_in.push_back(e);
      hole_out.push_back(circ_.get_next_edge(v, e));
      bin.insert(v);
      Expr beta = betas[i];
      if (type == OpType::NPhasedX || type == OpType::PhasedX) {
        Expr alpha = op->get_params()[1];
        if (!equiv_0(alpha, 2)) {
          sub1.add_op<unsigned>(OpType::Rz, {-alpha}, {i});
          sub3.add_op<unsigned>(OpType::Rz, {alpha}, {i});
        }
      }
      if (!equiv_0(beta, 2)) {
        sub2.add_op<unsigned>(OpType::Rz, {beta}, {i});
      } else if (!equiv_0(beta, 4)) {
        sub2.add_phase(-1);
      }
    } else {
      Edge interval_begin = intervals_[i].first;
      hole_in.push_back(interval_begin);
      hole_out.push_back(interval_begin);
    }
  }

  Circuit sub(circ_.n_qubits());
  sub.append(sub1);
  sub.add_op(OpType::NPhasedX, {-0.5, 0.5}, sub.all_qubits());
  sub.append(sub2);
  sub.add_op(OpType::NPhasedX, {0.5, 0.5}, sub.all_qubits());
  sub.append(sub3);
  Subcircuit hole{hole_in, hole_out, bin};

  BackupIntervals backup = backup_intervals();
  circ_.substitute(sub, hole, Circuit::VertexDeletion::Yes);
  restore_intervals(backup);

  skip_global_gates(2);
}

bool all_nullopt(const OptVertexVec& vec) {
  for (const OptVertex& v : vec) {
    if (v) return false;
  }
  return true;
}

}  // namespace Transforms

}  // namespace tket
