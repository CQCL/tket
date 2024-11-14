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

#include "tket/Converters/Converters.hpp"

namespace tket {

APState circuit_to_apstate(const Circuit& circ) {
  APState aps(circ.n_qubits());
  std::map<UnitID, unsigned> qb_ordering;
  for (const Qubit& q : circ.all_qubits())
    qb_ordering.insert({q, (unsigned)qb_ordering.size()});
  for (const Command& com : circ) {
    auto args = com.get_args();
    std::vector<unsigned> qbs;
    for (const UnitID& q : args) qbs.push_back(qb_ordering.at(q));
    aps.apply_gate(com.get_op_ptr()->get_type(), qbs);
  }
  return aps;
}

Circuit apstate_to_circuit(const APState& ap) {
  unsigned n_qbs = ap.A.rows();
  Circuit circ(n_qbs);
  circ.qubit_create_all();
  for (unsigned q = 0; q < n_qbs; ++q) {
    if (!ap.A(q, q))
      circ.add_op<unsigned>(OpType::H, {q});
    else if (ap.B(q))
      circ.add_op<unsigned>(OpType::X, {q});
  }
  for (unsigned trgt = 0; trgt < n_qbs; ++trgt) {
    if (ap.A(trgt, trgt)) {
      for (unsigned ctrl = trgt + 1; ctrl < n_qbs; ++ctrl)
        if (ap.A(trgt, ctrl)) circ.add_op<unsigned>(OpType::CX, {ctrl, trgt});
    }
  }
  for (unsigned q1 = 0; q1 < n_qbs; ++q1) {
    if (ap.A(q1, q1)) continue;
    for (unsigned q2 = q1 + 1; q2 < n_qbs; ++q2) {
      if (ap.E(q1, q2)) circ.add_op<unsigned>(OpType::CZ, {q1, q2});
    }
    switch (ap.P(q1) % 4) {
      case 1: {
        circ.add_op<unsigned>(OpType::S, {q1});
        break;
      }
      case 2: {
        circ.add_op<unsigned>(OpType::Z, {q1});
        break;
      }
      case 3: {
        circ.add_op<unsigned>(OpType::Sdg, {q1});
        break;
      }
      default: {
        break;
      }
    }
  }
  circ.add_phase(ap.phase);
  return circ;
}

APState tableau_to_apstate(SymplecticTableau tab) {
  unsigned n_qbs = tab.get_n_qubits();
  if (tab.get_n_rows() != n_qbs)
    throw std::logic_error(
        "tableau_to_apstate requires a tableau with n commuting rows for n "
        "qubits");
  MatrixXb fullmat = MatrixXb::Zero(n_qbs, 2 * n_qbs);
  /**
   * Gaussian elimination by the x matrix first ensures the bottom rows are only
   * Zs, i.e. describing rows of A. Reversing the columns of the x matrix
   * guarantees that each row has an X on at most one free qubit, simplifying
   * the code for finding E and P
   */
  for (unsigned c = 0; c < n_qbs; ++c) {
    fullmat.col(c) = tab.xmat.col(n_qbs - 1 - c);
  }
  fullmat.block(0, n_qbs, n_qbs, n_qbs) = tab.zmat;
  std::vector<std::pair<unsigned, unsigned>> row_ops =
      gaussian_elimination_row_ops(fullmat);
  for (const std::pair<unsigned, unsigned>& op : row_ops)
    tab.row_mult(op.first, op.second);

  MatrixXb A = MatrixXb::Zero(n_qbs, n_qbs);
  VectorXb B = VectorXb::Zero(n_qbs);
  MatrixXb E = MatrixXb::Zero(n_qbs, n_qbs);
  Eigen::VectorXi P = Eigen::VectorXi::Zero(n_qbs);

  unsigned n_leading = 0;
  for (unsigned r = n_qbs; r > 0;) {
    --r;
    bool is_a_row = true;
    for (unsigned c = 0; c < n_qbs; ++c) {
      if (tab.xmat(r, c)) {
        is_a_row = false;
        break;
      }
    }
    if (!is_a_row) break;
    unsigned first = 0;
    for (unsigned c = 0; c < n_qbs; ++c) {
      if (tab.zmat(r, c)) {
        first = c;
        break;
      }
    }
    A.row(first) = tab.zmat.row(r);
    B(first) = tab.phase(r);
    ++n_leading;
  }
  /**
   * Each free qubit is after all leaders connected to it, by reduced
   * row-echelon of A. Therefore Gaussian elimination of the x matrix in reverse
   * order gives rows corresponding to the columns of A of free qubits (plus the
   * free qubit itself).
   *
   * Then the corresponding row of the z matrix is exactly the free qubit's row
   * of E, from which we take the diagonal and phase vector to determine P.
   */
  for (unsigned r = 0; r < n_qbs - n_leading; ++r) {
    unsigned free = 0;
    for (unsigned c = n_qbs; c > 0;) {
      --c;
      if (tab.xmat(r, c)) {
        free = c;
        break;
      }
    }
    E.row(free) = tab.zmat.row(r);
    if (tab.zmat(r, free)) P(free) += 1;
    if (tab.phase(r)) P(free) += 2;
    E(free, free) = false;
  }

  return APState(A, B, E, P, 0);
}

SymplecticTableau apstate_to_tableau(const APState& ap) {
  unsigned n_qbs = ap.A.rows();
  MatrixXb xmat = MatrixXb::Zero(n_qbs, n_qbs);
  MatrixXb zmat = MatrixXb::Zero(n_qbs, n_qbs);
  VectorXb phase = VectorXb::Zero(n_qbs);

  for (unsigned q = 0; q < n_qbs; ++q) {
    if (ap.A(q, q)) {
      // Leading qubit
      zmat.row(q) = ap.A.row(q);
      phase(q) = ap.B(q);
    } else {
      // Free qubit
      xmat.row(q) = ap.A.col(q);
      xmat(q, q) = true;
      zmat.row(q) = ap.E.row(q);
      zmat(q, q) = (ap.P(q) % 2 == 1);
      phase(q) = (ap.P(q) % 4 > 1);
    }
  }

  return SymplecticTableau(xmat, zmat, phase);
}

}  // namespace tket