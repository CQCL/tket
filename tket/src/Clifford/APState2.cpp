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

#include "tket/Clifford/APState2.hpp"

#include <boost/foreach.hpp>
#include <tkassert/Assert.hpp>

#include "tket/OpType/OpTypeInfo.hpp"

namespace tket {

/*************************** APState IMPLEMENTATION ***************************/

APState::APState(
    const MatrixXb& A_, const VectorXb& B_, const MatrixXb& C_,
    const MatrixXb& E_, const Eigen::VectorXi& P_, const Expr& phase_)
    : A(A_), B(B_), C(C_), E(E_), P(P_), phase(phase_) {
  verify();
}

APState::APState(unsigned n_qubits)
    : A(MatrixXb::Identity(n_qubits, n_qubits)),
      B(VectorXb::Zero(n_qubits)),
      C(MatrixXb::Zero(n_qubits, n_qubits)),
      E(MatrixXb::Zero(n_qubits, n_qubits)),
      P(Eigen::VectorXi::Zero(n_qubits)),
      phase(0) {}

unsigned clifford_phase(const Complex& c) {
  unsigned res = 0;
  unsigned n_possibles = 0;
  if (c.real() > EPS) {
    res = 0;
    ++n_possibles;
  }
  if (c.imag() > EPS) {
    res = 1;
    ++n_possibles;
  }
  if (c.real() < -EPS) {
    res = 2;
    ++n_possibles;
  }
  if (c.imag() < -EPS) {
    res = 3;
    ++n_possibles;
  }
  TKET_ASSERT(n_possibles == 1);
  return res;
}

APState::APState(const Eigen::VectorXcd& sv) {
  unsigned n_qbs = 0;
  while (sv.size() > (1 << n_qbs)) ++n_qbs;
  TKET_ASSERT(sv.size() == (1 << n_qbs));

  // Find non-zero entries as a vector space and offset
  unsigned z0 = 0;
  for (unsigned x = 0; x < sv.size(); ++x) {
    if (sv(x) != 0.) {
      z0 = x;
      break;
    }
  }
  std::vector<unsigned> offsets;
  unsigned n_non_zero = 0;
  for (unsigned x = 1; x < sv.size(); ++x) {
    if (sv(z0 ^ x) != 0.) {
      ++n_non_zero;
      if (n_non_zero == (1u << offsets.size())) offsets.push_back(x);
    }
  }

  // Find A as the dual space
  MatrixXb offset_mat(offsets.size(), n_qbs);
  for (unsigned r = 0; r < offsets.size(); ++r) {
    unsigned off = offsets.at(r);
    for (unsigned c = 0; c < n_qbs; ++c) {
      // Binary encoding of offsets in reverse order to guarantee free qubits
      // are the later ones, meaning we produce A in row echelon form
      offset_mat(r, c) = (off & (1 << c)) != 0;
    }
  }
  std::vector<std::pair<unsigned, unsigned>> row_ops =
      gaussian_elimination_row_ops(offset_mat);
  for (const std::pair<unsigned, unsigned>& op : row_ops) {
    for (unsigned j = 0; j < n_qbs; ++j) {
      offset_mat(op.second, j) ^= offset_mat(op.first, j);
    }
  }
  std::map<unsigned, unsigned> mat_leaders;
  A = MatrixXb::Zero(n_qbs, n_qbs);
  for (unsigned c = 0; c < n_qbs; ++c) {
    bool free_qubit = false;
    for (unsigned r = 0; r < offsets.size(); ++r) {
      if (offset_mat(r, c)) {
        auto inserted = mat_leaders.insert({r, c});
        if (inserted.second) {
          free_qubit = true;
          break;
        } else {
          // Reverse bit orderings back to normal
          A(n_qbs - 1 - c, n_qbs - 1 - inserted.first->second) = true;
        }
      }
    }
    A(n_qbs - 1 - c, n_qbs - 1 - c) = !free_qubit;
  }

  VectorXb z0_vec(n_qbs);
  for (unsigned i = 0; i < n_qbs; ++i) {
    z0_vec(i) = ((z0 & (1 << (n_qbs - 1 - i))) != 0);
  }
  B = A * z0_vec;

  E = MatrixXb::Zero(n_qbs, n_qbs);
  P = Eigen::VectorXi::Zero(n_qbs);
  unsigned neutral_z0 = z0;  // Index with 0 for all free qubits
  std::map<unsigned, unsigned> offset_for_free;
  for (const std::pair<const unsigned, unsigned>& row_qfree : mat_leaders) {
    unsigned offset = 0;
    for (unsigned i = 0; i < n_qbs; ++i)
      if (offset_mat(row_qfree.first, i)) offset += (1 << i);
    unsigned qfree = n_qbs - 1 - row_qfree.second;
    offset_for_free.insert({qfree, offset});
    if ((neutral_z0 & (1 << row_qfree.second)) != 0)
      neutral_z0 = neutral_z0 ^ offset;
  }
  for (const std::pair<const unsigned, unsigned>& qfree_offset :
       offset_for_free) {
    unsigned qfree = qfree_offset.first;
    unsigned offset = qfree_offset.second;
    auto complex_phase = sv(neutral_z0 ^ offset) / sv(neutral_z0);
    P(qfree) = clifford_phase(complex_phase);
    for (const std::pair<const unsigned, unsigned>& qfree_offset2 :
         offset_for_free) {
      if (qfree_offset == qfree_offset2) break;  // Only go up to solved phases
      unsigned qfree2 = qfree_offset2.first;
      unsigned offset2 = qfree_offset2.second;
      auto complex_phase_cz =
          sv(neutral_z0 ^ offset ^ offset2) / sv(neutral_z0);
      E(qfree, qfree2) = E(qfree2, qfree) =
          ((unsigned)(clifford_phase(complex_phase_cz) - P(qfree) - P(qfree2)) %
           4) == 2;
    }
  }

  phase = Expr(std::arg(sv(neutral_z0)) / PI);
}

APState::APState(const Eigen::MatrixXcd& dm) {
  Eigen::VectorXcd dm_as_vec = dm.reshaped();
  APState pure_double(dm_as_vec);
  unsigned n_qbs = pure_double.A.cols() / 2;
  A = MatrixXb::Zero(n_qbs, n_qbs);
  B = VectorXb::Zero(n_qbs);
  C = MatrixXb::Zero(n_qbs, n_qbs);

  MatrixXb full_mat(2 * n_qbs, 2 * n_qbs + 1);
  full_mat.block(0, 0, 2 * n_qbs, 2 * n_qbs) = pure_double.A;
  full_mat.col(2 * n_qbs) = pure_double.B;

  // This full matrix is generated by AI, IA, and CC. Gaussian elimination puts
  // zeros at the bottom, followed by IA.
  std::vector<std::pair<unsigned, unsigned>> row_ops =
      gaussian_elimination_row_ops(full_mat);
  for (const std::pair<unsigned, unsigned>& op : row_ops) {
    for (unsigned j = 0; j < 2 * n_qbs + 1; ++j) {
      full_mat(op.second, j) ^= full_mat(op.first, j);
    }
  }

  unsigned first_zero = 0;
  for (unsigned r = n_qbs; r > 0;) {
    --r;
    bool is_zero = true;
    // By symmetry of AI, IA, and CC, the bottom non-empty row must have some
    // component on the right side
    for (unsigned c = n_qbs; c < 2 * n_qbs; ++c) {
      if (full_mat(r, c)) {
        is_zero = false;
        break;
      }
    }
    if (!is_zero) {
      first_zero = r + 1;
      break;
    }
  }
  unsigned first_right = 0;
  for (unsigned r = first_zero; r > 0;) {
    --r;
    bool is_right_only = true;
    for (unsigned c = 0; c < n_qbs; ++c) {
      if (full_mat(r, c)) {
        is_right_only = false;
        break;
      }
    }
    if (!is_right_only) {
      first_right = r + 1;
      break;
    }
  }
  A.block(0, 0, first_zero - first_right, n_qbs) =
      full_mat.block(first_right, n_qbs, first_zero - first_right, n_qbs);
  B.head(first_zero - first_right) =
      full_mat.block(first_right, 2 * n_qbs, first_zero - first_right, n_qbs);

  // Flip the column order and reduce the remaining rows to obtain CC above AI;
  // get the row combinations from the reordered matrix but apply them to the
  // matrix with the correct column ordering
  MatrixXb remaining_rows(first_right, 2 * n_qbs);
  for (unsigned r = 0; r < first_right; ++r) {
    for (unsigned c = 0; c < 2 * n_qbs; ++c) {
      remaining_rows(r, c) = full_mat(r, 2 * n_qbs - 1 - c);
    }
  }
  row_ops = gaussian_elimination_row_ops(remaining_rows);
  for (const std::pair<unsigned, unsigned>& op : row_ops) {
    for (unsigned j = 0; j < 2 * n_qbs; ++j) {
      full_mat(op.second, j) ^= full_mat(op.first, j);
    }
  }

  unsigned first_left = 0;
  for (unsigned r = first_right; r > 0;) {
    --r;
    bool is_left_only = true;
    for (unsigned c = n_qbs; c < 2 * n_qbs; ++c) {
      if (full_mat(r, c)) {
        is_left_only = false;
        break;
      }
    }
    if (!is_left_only) {
      first_left = r + 1;
      break;
    }
  }
  C.block(0, 0, first_left, n_qbs) = full_mat.block(0, 0, first_left, n_qbs);

  /**
   * Recall that pure_double.A is generated by AI, IA, and CC. The statevector
   * constructor builds this in normal form, so it is in reduced row-echelon
   * form. The first block of qubits may have leaders for leaders or mixed
   * qubits in the density matrix, but the second block only has them for true
   * leaders in the density matrix. This means pure_double.E and pure_double.P
   * are not just the direct sum of two copies of the true E and P, but the
   * bottom segments contain them exactly.
   */
  E = pure_double.E.block(n_qbs, n_qbs, n_qbs, n_qbs);
  P = pure_double.P.tail(n_qbs);

  phase = 0;
}

void APState::verify() const {
  // Check dimensions all agree
  unsigned n_qubits = A.rows();
  TKET_ASSERT(A.cols() == n_qubits);
  TKET_ASSERT(B.size() == n_qubits);
  TKET_ASSERT(C.rows() == n_qubits);
  TKET_ASSERT(C.cols() == n_qubits);
  TKET_ASSERT(E.rows() == n_qubits);
  TKET_ASSERT(E.cols() == n_qubits);
  TKET_ASSERT(P.size() == n_qubits);

  for (unsigned r = 0; r < n_qubits; ++r) {
    for (unsigned c = 0; c < r; ++c) {
      // Check E is symmetric
      TKET_ASSERT(E(r, c) == E(c, r));
    }
    // Check diagonals of E are zero
    TKET_ASSERT(!E(r, r));
  }
}

VectorXb z2_mult(const MatrixXb& M, const VectorXb& v) {
  VectorXb res = VectorXb::Zero(M.rows());
  for (unsigned i = 0; i < M.cols(); ++i) {
    if (v(i)) {
      for (unsigned j = 0; j < M.rows(); ++j) res(j) ^= M(j, i);
    }
  }
  return res;
}

Eigen::VectorXcd APState::to_statevector() const {
  unsigned n_qubits = A.cols();
  Eigen::VectorXcd sv = Eigen::VectorXcd::Zero(1 << n_qubits);
  unsigned n_terms = 0;
  Complex g_phase = std::exp(i_ * PI * eval_expr(phase).value());
  for (unsigned x = 0; x < (1u << n_qubits); ++x) {
    VectorXb x_binary = VectorXb::Zero(n_qubits);
    for (unsigned i = 0; i < n_qubits; ++i) {
      unsigned mask = 1 << (n_qubits - 1 - i);
      x_binary(i) = ((x & mask) != 0);
    }

    if (z2_mult(A, x_binary) == B) {
      ++n_terms;
      unsigned i_phases = 0;
      for (unsigned q = 0; q < n_qubits; ++q) {
        if (x_binary(q)) {
          i_phases += P(q);
          for (unsigned q2 = q; q2 < n_qubits; ++q2) {
            if (E(q, q2) && x_binary(q2)) i_phases += 2;
          }
        }
      }
      switch (i_phases % 4) {
        case 0: {
          sv(x) = g_phase;
          break;
        }
        case 1: {
          sv(x) = i_ * g_phase;
          break;
        }
        case 2: {
          sv(x) = -g_phase;
          break;
        }
        case 3: {
          sv(x) = -i_ * g_phase;
          break;
        }
        default: {
          TKET_ASSERT(false);
        }
      }
    }
  }
  return pow(n_terms, -0.5) * sv;
}

Eigen::MatrixXcd APState::to_density_matrix() const {
  unsigned n_qubits = A.cols();
  Eigen::MatrixXcd dm = Eigen::MatrixXcd::Zero(1 << n_qubits, 1 << n_qubits);
  for (unsigned x = 0; x < (1u << n_qubits); ++x) {
    VectorXb x_binary = VectorXb::Zero(n_qubits);
    for (unsigned i = 0; i < n_qubits; ++i) {
      unsigned mask = 1 << (n_qubits - 1 - i);
      x_binary(i) = ((x & mask) != 0);
    }

    if (z2_mult(A, x_binary) == B) {
      MatrixXb Cx = z2_mult(C, x_binary);
      for (unsigned y = 0; y < (1u << n_qubits); ++y) {
        VectorXb y_binary = VectorXb::Zero(n_qubits);
        for (unsigned i = 0; i < n_qubits; ++i) {
          unsigned mask = 1 << (n_qubits - 1 - i);
          y_binary(i) = ((y & mask) != 0);
        }

        if ((z2_mult(A, y_binary) == B) && (z2_mult(C, y_binary) == Cx)) {
          unsigned i_phases = 0;
          for (unsigned q = 0; q < n_qubits; ++q) {
            if (x_binary(q)) {
              i_phases += P(q);
              for (unsigned q2 = q; q2 < n_qubits; ++q2) {
                if (E(q, q2) && x_binary(q2)) i_phases += 2;
              }
            }
            if (y_binary(q)) {
              i_phases -= P(q);
              for (unsigned q2 = q; q2 < n_qubits; ++q2) {
                if (E(q, q2) && y_binary(q2)) i_phases += 2;
              }
            }
          }
          switch (i_phases % 4) {
            case 0: {
              dm(x, y) = 1;
              break;
            }
            case 1: {
              dm(x, y) = i_;
              break;
            }
            case 2: {
              dm(x, y) = -1;
              break;
            }
            case 3: {
              dm(x, y) = -i_;
              break;
            }
            default: {
              TKET_ASSERT(false);
            }
          }
        }
      }
    }
  }
  return dm / dm.trace();
}

void APState::apply_CZ(unsigned ctrl, unsigned trgt) {
  E(ctrl, trgt) ^= true;
  E(trgt, ctrl) ^= true;
}

void APState::apply_S(unsigned q) { P(q) += 1; }

void APState::apply_V(unsigned q) {
  unsigned n_qbs = A.cols();
  std::list<unsigned> A_rows;
  std::list<unsigned> C_rows;
  for (unsigned r = 0; r < n_qbs; ++r) {
    if (A(r, q)) A_rows.push_back(r);
    if (C(r, q)) C_rows.push_back(r);
  }
  if ((unsigned)P(q) % 2 == 0) {
    bool a = ((unsigned)P(q) % 4) == 2;
    if (A_rows.empty()) {
      for (unsigned q2 = 0; q2 < n_qbs; ++q2) {
        if (E(q, q2)) {
          // Update local phase on neighbours
          P(q2) += a ? 3 : 1;
          // Local complementation between neighbours
          for (unsigned q3 = q2 + 1; q3 < n_qbs; ++q3) {
            if (E(q, q3)) {
              E(q2, q3) ^= true;
              E(q3, q2) ^= true;
            }
          }
          // Add connections to all C_rows
          for (unsigned r : C_rows) {
            C(r, q2) ^= true;
          }
        }
      }
      // Global phase
      phase += a ? .25 : -.25;
    } else {
      // Choose one of the red spiders to remove
      unsigned r = A_rows.back();
      A_rows.pop_back();
      // Stratify neighbourhoods of r (in A) and q (in E)
      std::list<unsigned> just_r;
      std::list<unsigned> just_q;
      std::list<unsigned> both;
      for (unsigned q2 = 0; q2 < n_qbs; ++q2) {
        if (q2 == q) continue;
        if (A(r, q2)) {
          if (E(q, q2))
            both.push_back(q2);
          else
            just_r.push_back(q2);
        } else if (E(q, q2))
          just_q.push_back(q2);
      }
      // Update A and B
      for (unsigned ar : A_rows) {
        A(ar, q) = false;
        for (unsigned q2 : just_r) A(ar, q2) ^= true;
        for (unsigned q2 : both) A(ar, q2) ^= true;
        B(ar) ^= B(r);
      }
      // Update C
      for (unsigned cr : C_rows) {
        C(cr, q) = false;
        for (unsigned q2 : just_r) C(cr, q2) ^= true;
        for (unsigned q2 : both) C(cr, q2) ^= true;
      }
      // Update E and P
      for (unsigned rn : just_r) {
        // Complementation within just_r
        for (unsigned rn2 : just_r) {
          E(rn, rn2) ^= true;
        }
        // Reset diagonal
        E(rn, rn) = false;
        // Complementation between just_r and just_q
        for (unsigned qn : just_q) {
          E(rn, qn) ^= true;
          E(qn, rn) ^= true;
        }
        // Connect to q
        E(rn, q) = true;
        E(q, rn) = true;
        // Local phases
        P(rn) += (a ^ B(r)) ? 1 : 3;
      }
      for (unsigned bn : both) {
        // Complementation within both
        for (unsigned bn2 : both) {
          E(bn, bn2) ^= true;
        }
        // Reset diagonal
        E(bn, bn) = false;
        // Complementation between both and just_q
        for (unsigned qn : just_q) {
          E(bn, qn) ^= true;
          E(qn, bn) ^= true;
        }
        // Local phases
        P(bn) += a ? 3 : 1;
      }
      for (unsigned qn : just_q) {
        // Remove connection to q
        E(q, qn) = false;
        E(qn, q) = false;
        // Local phases
        P(qn) += B(r) ? 2 : 0;
      }
      // Local phase on q
      P(q) = B(r) ? 1 : 3;
      // Global phase
      if (B(r)) phase += a ? .5 : 1.5;
      // Remove row r from A and B
      for (unsigned i = 0; i < n_qbs; ++i) {
        A(r, i) = false;
      }
      B(r) = false;
    }
  } else {
    bool a = ((unsigned)P(q) % 4) == 3;
    if (!A_rows.empty()) {
      // Choose one of the red spiders to remove
      unsigned r = A_rows.back();
      A_rows.pop_back();
      // Stratify neighbourhoods of r (in A) and q (in E)
      std::list<unsigned> just_r;
      std::list<unsigned> just_q;
      std::list<unsigned> both;
      for (unsigned q2 = 0; q2 < n_qbs; ++q2) {
        if (q2 == q) continue;
        if (A(r, q2)) {
          if (E(q, q2))
            both.push_back(q2);
          else
            just_r.push_back(q2);
        } else if (E(q, q2))
          just_q.push_back(q2);
      }
      // Update A and B
      for (unsigned ar : A_rows) {
        A(ar, q) = false;
        for (unsigned q2 : just_r) A(ar, q2) ^= true;
        for (unsigned q2 : both) A(ar, q2) ^= true;
        B(ar) ^= B(r);
      }
      // Update C
      for (unsigned cr : C_rows) {
        C(cr, q) = false;
        for (unsigned q2 : just_r) C(cr, q2) ^= true;
        for (unsigned q2 : both) C(cr, q2) ^= true;
      }
      // Update E and P
      for (unsigned rn : just_r) {
        // Complementation between just_r and just_q
        for (unsigned qn : just_q) {
          E(rn, qn) ^= true;
          E(qn, rn) ^= true;
        }
        // Complementation between just_r and both
        for (unsigned bn : both) {
          E(rn, bn) ^= true;
          E(bn, rn) ^= true;
        }
        // Connect to q
        E(rn, q) = true;
        E(q, rn) = true;
        // Local phases
        P(rn) += a ? 2 : 0;
      }
      for (unsigned bn : both) {
        // Complementation between both and just_q
        for (unsigned qn : just_q) {
          E(bn, qn) ^= true;
          E(qn, bn) ^= true;
        }
        // Local phases
        P(bn) += (a ^ B(r)) ? 0 : 2;
      }
      for (unsigned qn : just_q) {
        // Remove connection to q
        E(qn, q) = false;
        E(q, qn) = false;
        // Local phases
        P(qn) += B(r) ? 2 : 0;
      }
      // Local phase on q
      P(q) = B(r) ? 1 : 3;
      // Global phase
      if (a && B(r)) phase += 1;
      // Remove row r from A and B
      for (unsigned i = 0; i < n_qbs; ++i) {
        A(r, i) = false;
      }
      B(r) = false;
    } else if (!C_rows.empty()) {
      // Choose one of the discards to pivot around
      unsigned d = C_rows.back();
      C_rows.pop_back();
      // Stratify neighbourhoods of d (in C) and q (in E)
      std::list<unsigned> just_d;
      std::list<unsigned> just_q;
      std::list<unsigned> both;
      for (unsigned q2 = 0; q2 < n_qbs; ++q2) {
        if (q2 == q) continue;
        if (C(d, q2)) {
          if (E(q, q2))
            both.push_back(q2);
          else
            just_d.push_back(q2);
        } else if (E(q, q2))
          just_q.push_back(q2);
      }
      // Update C
      for (unsigned cr : C_rows) {
        C(cr, q) = false;
        for (unsigned q2 : just_d) C(cr, q2) ^= true;
        for (unsigned q2 : both) C(cr, q2) ^= true;
      }
      // Update E and P
      for (unsigned dn : just_d) {
        // Complementation between just_d and just_q
        for (unsigned qn : just_q) {
          E(dn, qn) ^= true;
          E(qn, dn) ^= true;
        }
        // Complementation between just_d and both
        for (unsigned bn : both) {
          E(dn, bn) ^= true;
          E(bn, dn) ^= true;
        }
        // Connect to q
        E(dn, q) = true;
        E(q, dn) = true;
        // Remove connection with d
        C(d, dn) = false;
        // Local phases
        P(dn) += a ? 2 : 0;
      }
      for (unsigned bn : both) {
        // Complementation between both and just_q
        for (unsigned qn : just_q) {
          E(bn, qn) ^= true;
          E(qn, bn) ^= true;
        }
        // Local phases
        P(bn) += a ? 0 : 2;
      }
      for (unsigned qn : just_q) {
        // Connect to d
        C(d, qn) = true;
        // Remove connection to q
        E(qn, q) = false;
        E(q, qn) = false;
        // No local phase change
      }
      // Local phase on q
      P(q) += 3;
      // No global phase change
    } else {
      VectorXb new_row = VectorXb::Zero(n_qbs);
      new_row(q) = true;
      for (unsigned q2 = 0; q2 < n_qbs; ++q2) {
        if (E(q, q2)) {
          // Connect to new red spider
          new_row(q2) = true;
          // Local phase on neighbours
          P(q2) += a ? 1 : 3;
          // Local complementation between neighbours
          for (unsigned q3 = q2 + 1; q3 < n_qbs; ++q3) {
            if (E(q, q3)) {
              E(q2, q3) ^= true;
              E(q3, q2) ^= true;
            }
          }
          // Remove connection with q
          E(q, q2) = false;
          E(q2, q) = false;
        }
      }
      // Reset local phase on q
      P(q) = 0;
      // Global phase
      phase += a ? 1.5 : 0;
      // Add new_row to A and B
      MatrixXb combined_mat = MatrixXb::Zero(n_qbs + 1, n_qbs + 1);
      combined_mat.block(0, 0, n_qbs, n_qbs) = A;
      combined_mat.block(0, n_qbs, n_qbs, 1) = B;
      combined_mat.block(n_qbs, 0, 1, n_qbs) = new_row.transpose();
      combined_mat(n_qbs, n_qbs) = a;
      std::vector<std::pair<unsigned, unsigned>> row_ops =
          gaussian_elimination_row_ops(combined_mat);
      for (const std::pair<unsigned, unsigned>& op : row_ops) {
        for (unsigned j = 0; j < n_qbs + 1; ++j) {
          combined_mat(op.second, j) ^= combined_mat(op.first, j);
        }
      }
      A = combined_mat.block(0, 0, n_qbs, n_qbs);
      B = combined_mat.block(0, n_qbs, n_qbs, 1);
    }
  }
}

void APState::apply_X(unsigned q) {
  // Push through the local phase
  phase += P(q) * .5;
  P(q) = -P(q);
  // Pushing through the CZs adds Zs on neighbours
  for (unsigned q2 = 0; q2 < E.cols(); ++q2) {
    if (E(q, q2)) P(q2) += 2;
  }
  // Pushing through adjacency matrix adds Xs onto connected reds
  for (unsigned r = 0; r < A.rows(); ++r) {
    if (A(r, q)) B(r) ^= true;
  }
}

void APState::apply_gate(OpType type, const std::vector<unsigned>& qbs) {
  switch (type) {
    case OpType::Z: {
      apply_S(qbs.at(0));
      apply_S(qbs.at(0));
      break;
    }
    case OpType::X: {
      apply_X(qbs.at(0));
      break;
    }
    case OpType::Y: {
      apply_S(qbs.at(0));
      apply_S(qbs.at(0));
      apply_X(qbs.at(0));
      phase += .5;
      break;
    }
    case OpType::S: {
      apply_S(qbs.at(0));
      break;
    }
    case OpType::Sdg: {
      apply_S(qbs.at(0));
      apply_S(qbs.at(0));
      apply_S(qbs.at(0));
      break;
    }
    case OpType::V: {
      apply_V(qbs.at(0));
      break;
    }
    case OpType::SX: {
      apply_V(qbs.at(0));
      phase += .25;
      break;
    }
    case OpType::Vdg: {
      apply_V(qbs.at(0));
      apply_X(qbs.at(0));
      phase += .5;
      break;
    }
    case OpType::SXdg: {
      apply_V(qbs.at(0));
      apply_X(qbs.at(0));
      phase += .25;
      break;
    }
    case OpType::H: {
      apply_S(qbs.at(0));
      apply_V(qbs.at(0));
      apply_S(qbs.at(0));
      break;
    }
    case OpType::CX: {
      apply_S(qbs.at(1));
      apply_V(qbs.at(1));
      apply_S(qbs.at(1));
      apply_CZ(qbs.at(0), qbs.at(1));
      apply_S(qbs.at(1));
      apply_V(qbs.at(1));
      apply_S(qbs.at(1));
      break;
    }
    case OpType::CY: {
      apply_V(qbs.at(1));
      apply_CZ(qbs.at(0), qbs.at(1));
      apply_V(qbs.at(1));
      apply_X(qbs.at(1));
      phase += .5;
      break;
    }
    case OpType::CZ: {
      apply_CZ(qbs.at(0), qbs.at(1));
      break;
    }
    case OpType::ZZMax: {
      apply_S(qbs.at(0));
      apply_S(qbs.at(1));
      apply_CZ(qbs.at(0), qbs.at(1));
      phase -= .25;
      break;
    }
    case OpType::ECR: {
      apply_S(qbs.at(0));
      apply_X(qbs.at(0));
      apply_S(qbs.at(1));
      apply_V(qbs.at(1));
      apply_CZ(qbs.at(0), qbs.at(1));
      apply_S(qbs.at(1));
      apply_V(qbs.at(1));
      apply_S(qbs.at(1));
      phase += -.25;
      break;
    }
    case OpType::ISWAPMax: {
      apply_V(qbs.at(0));
      apply_S(qbs.at(1));
      apply_V(qbs.at(1));
      apply_CZ(qbs.at(0), qbs.at(1));
      apply_V(qbs.at(0));
      apply_V(qbs.at(1));
      apply_CZ(qbs.at(0), qbs.at(1));
      apply_V(qbs.at(0));
      apply_V(qbs.at(1));
      apply_S(qbs.at(1));
      phase += 1.;
      break;
    }
    case OpType::SWAP: {
      apply_S(qbs.at(0));
      apply_V(qbs.at(0));
      apply_S(qbs.at(0));
      apply_CZ(qbs.at(0), qbs.at(1));
      apply_S(qbs.at(0));
      apply_V(qbs.at(0));
      apply_S(qbs.at(0));
      apply_S(qbs.at(1));
      apply_V(qbs.at(1));
      apply_S(qbs.at(1));
      apply_CZ(qbs.at(0), qbs.at(1));
      apply_S(qbs.at(0));
      apply_V(qbs.at(0));
      apply_S(qbs.at(0));
      apply_S(qbs.at(1));
      apply_V(qbs.at(1));
      apply_S(qbs.at(1));
      apply_CZ(qbs.at(0), qbs.at(1));
      apply_S(qbs.at(0));
      apply_V(qbs.at(0));
      apply_S(qbs.at(0));
      break;
    }
    case OpType::BRIDGE: {
      apply_S(qbs.at(2));
      apply_V(qbs.at(2));
      apply_S(qbs.at(2));
      apply_CZ(qbs.at(0), qbs.at(2));
      apply_S(qbs.at(2));
      apply_V(qbs.at(2));
      apply_S(qbs.at(2));
      break;
    }
    case OpType::noop: {
      break;
    }
    case OpType::Reset: {
      unsigned q = qbs.at(0);
      discard_qubit(q);
      unsigned new_qb = init_qubit();
      // Swap qubits to preserve the ordering
      if (q != new_qb) {
        unsigned n_qbs = A.cols();
        for (unsigned r = 0; r < n_qbs; ++r) {
          bool temp = A(r, q);
          A(r, q) = A(r, new_qb);
          A(r, new_qb) = temp;
        }
        for (unsigned r = 0; r < n_qbs; ++r) {
          bool temp = C(r, q);
          C(r, q) = C(r, new_qb);
          C(r, new_qb) = temp;
        }
        for (unsigned r = 0; r < n_qbs; ++r) {
          bool temp = E(r, q);
          E(r, q) = E(r, new_qb);
          E(r, new_qb) = temp;
        }
        for (unsigned c = 0; c < n_qbs; ++c) {
          bool temp = E(q, c);
          E(q, c) = E(new_qb, c);
          E(new_qb, c) = temp;
        }
        int temp = P(q);
        P(q) = P(new_qb);
        P(new_qb) = temp;
      }
      break;
    }
    case OpType::Collapse: {
      collapse_qubit(qbs.at(0));
      break;
    }
    case OpType::Phase: {
      throw std::logic_error(
          "OpType::Phase cannot be applied via APState::apply_gate");
    }
    default: {
      throw BadOpType(
          "Cannot be applied to a APState: not a Clifford gate", type);
    }
  }
}

unsigned APState::init_qubit() {
  unsigned n_qbs = A.cols();
  A.conservativeResize(n_qbs + 1, n_qbs + 1);
  A.col(n_qbs) = VectorXb::Zero(n_qbs + 1);
  A.row(n_qbs) = VectorXb::Zero(n_qbs + 1);
  A(n_qbs, n_qbs) = true;
  B.conservativeResize(n_qbs + 1);
  B(n_qbs) = false;
  C.conservativeResize(n_qbs + 1, n_qbs + 1);
  C.col(n_qbs) = VectorXb::Zero(n_qbs + 1);
  C.row(n_qbs) = VectorXb::Zero(n_qbs + 1);
  E.conservativeResize(n_qbs + 1, n_qbs + 1);
  E.col(n_qbs) = VectorXb::Zero(n_qbs + 1);
  E.row(n_qbs) = VectorXb::Zero(n_qbs + 1);
  P.conservativeResize(n_qbs + 1);
  P(n_qbs) = 0;
  return n_qbs;
}

unsigned APState::post_select(unsigned q) {
  unsigned n_qbs = A.cols();
  if (q >= n_qbs)
    throw std::invalid_argument("APState: post-selecting a non-existent qubit");
  MatrixXb AB = MatrixXb::Zero(n_qbs, n_qbs);
  AB.block(0, 0, n_qbs, n_qbs) = A;
  AB.col(q) = A.col(n_qbs - 1);
  AB.col(n_qbs - 1) = B;
  std::vector<std::pair<unsigned, unsigned>> row_ops =
      gaussian_elimination_row_ops(AB);
  for (const std::pair<unsigned, unsigned>& op : row_ops) {
    for (unsigned j = 0; j < n_qbs; ++j) {
      AB(op.second, j) ^= AB(op.first, j);
    }
  }
  A.resize(n_qbs - 1, n_qbs - 1);
  A = AB.block(0, 0, n_qbs - 1, n_qbs - 1);
  B.resize(n_qbs - 1);
  B = AB.block(0, n_qbs - 1, n_qbs - 1, 1);
  C.col(q) = C.col(n_qbs - 1);
  row_ops = gaussian_elimination_row_ops(C);
  for (const std::pair<unsigned, unsigned>& op : row_ops) {
    for (unsigned j = 0; j < n_qbs; ++j) {
      C(op.second, j) ^= C(op.first, j);
    }
  }
  C.conservativeResize(n_qbs - 1, n_qbs - 1);
  E.col(q) = E.col(n_qbs - 1);
  E.row(q) = E.row(n_qbs - 1);
  E.conservativeResize(n_qbs - 1, n_qbs - 1);
  P(q) = P(n_qbs - 1);
  P.conservativeResize(n_qbs - 1);
  return n_qbs - 1;
}

void APState::collapse_qubit(unsigned q) {
  unsigned n_qbs = A.cols();
  MatrixXb C_ext = MatrixXb::Zero(n_qbs + 1, n_qbs);
  C_ext.block(0, 0, n_qbs, n_qbs) = C;
  C_ext(n_qbs, q) = 1;
  std::vector<std::pair<unsigned, unsigned>> row_ops =
      gaussian_elimination_row_ops(C_ext);
  for (const std::pair<unsigned, unsigned>& op : row_ops) {
    for (unsigned j = 0; j < n_qbs; ++j) {
      C_ext(op.second, j) ^= C_ext(op.first, j);
    }
  }
  C = C_ext.block(0, 0, n_qbs, n_qbs);
}

unsigned APState::discard_qubit(unsigned q) {
  collapse_qubit(q);
  apply_V(q);
  collapse_qubit(q);
  return post_select(q);
}

void APState::normal_form() {
  unsigned n_qbs = A.cols();
  // Get A into reduced row-echelon form
  std::vector<std::pair<unsigned, unsigned>> row_ops =
      gaussian_elimination_row_ops(A);
  for (const std::pair<unsigned, unsigned>& op : row_ops) {
    for (unsigned j = 0; j < n_qbs; ++j) {
      A(op.second, j) ^= A(op.first, j);
    }
    B(op.second) ^= B(op.first);
  }
  // Identify leading qubits
  std::map<unsigned, unsigned> leader_to_row;
  for (unsigned r = 0; r < n_qbs; ++r) {
    bool leader_found = false;
    for (unsigned c = 0; c < n_qbs; ++c) {
      if (A(r, c)) {
        leader_found = true;
        leader_to_row.insert({c, r});
        break;
      }
    }
    if (!leader_found) break;
  }
  // Remove leaders from C
  for (unsigned r = 0; r < n_qbs; ++r) {
    for (const std::pair<const unsigned, unsigned>& lr : leader_to_row) {
      if (C(r, lr.first)) {
        for (unsigned c = 0; c < n_qbs; ++c) {
          C(r, c) ^= A(lr.second, c);
        }
      }
    }
  }
  // Remove leaders from E
  for (const std::pair<const unsigned, unsigned>& lr : leader_to_row) {
    for (unsigned q2 = lr.first; q2 < n_qbs; ++q2) {
      if (E(lr.first, q2)) {
        E(lr.first, q2) = false;
        E(q2, lr.first) = false;
        for (unsigned q3 = lr.first + 1; q3 < n_qbs; ++q3) {
          if (A(lr.second, q3)) {
            E(q2, q3) ^= true;
            E(q3, q2) ^= true;
          }
        }
        P(q2) += (B(lr.second) ^ A(lr.second, q2)) ? 2 : 0;
      }
    }
  }
  // Remove leaders from P
  for (const std::pair<const unsigned, unsigned>& lr : leader_to_row) {
    if ((unsigned)P(lr.first) % 2 == 1) {
      // Complementation within neighbours under A (excluding the leader)
      for (unsigned q2 = lr.first + 1; q2 < n_qbs; ++q2) {
        if (A(lr.second, q2)) {
          for (unsigned q3 = q2 + 1; q3 < n_qbs; ++q3) {
            if (A(lr.second, q3)) {
              E(q2, q3) ^= true;
              E(q3, q2) ^= true;
            }
          }
          // Local phases of neighbours
          P(q2) += (P(lr.first) % 4) + (B(lr.second) ? 2 : 0);
        }
      }
    } else if ((unsigned)P(lr.first) % 4 == 2) {
      // Local phases of neighbours
      for (unsigned q2 = lr.first + 1; q2 < n_qbs; ++q2) {
        if (A(lr.second, q2)) P(q2) += 2;
      }
    }
    // Global phase
    if (B(lr.second)) phase += P(lr.first) * .5;
    // Reset P
    P(lr.first) = 0;
  }
  // Get C into reduced row-echelon form
  row_ops = gaussian_elimination_row_ops(C);
  for (const std::pair<unsigned, unsigned>& op : row_ops) {
    for (unsigned j = 0; j < n_qbs; ++j) {
      C(op.second, j) ^= C(op.first, j);
    }
  }
  // Identify mixed qubits
  std::map<unsigned, unsigned> mixed_to_row;
  for (unsigned r = 0; r < n_qbs; ++r) {
    bool mixed_found = false;
    for (unsigned c = 0; c < n_qbs; ++c) {
      if (C(r, c)) {
        mixed_found = true;
        mixed_to_row.insert({c, r});
        break;
      }
    }
    if (!mixed_found) break;
  }
  // Remove E connections between mixed qubits
  for (auto mr1 = mixed_to_row.begin(); mr1 != mixed_to_row.end(); ++mr1) {
    for (auto mr2 = mixed_to_row.begin(); mr2 != mr1; ++mr2) {
      if (E(mr1->first, mr2->first)) {
        // Local complementation along the sum of their rows in C
        // mr2->first < mr1->first, and they are the first entries in their rows
        for (unsigned i = mr2->first; i < n_qbs; ++i) {
          if (C(mr1->second, i) ^ C(mr2->second, i)) {
            for (unsigned j = i + 1; j < n_qbs; ++j) {
              if (C(mr1->second, j) ^ C(mr2->second, j)) {
                E(i, j) ^= true;
                E(j, i) ^= true;
              }
            }
            // Add local phases around neighbourhood
            P(i) += 1;
          }
        }
      }
    }
  }
  // Remove mixed qubits from P
  for (const std::pair<const unsigned, unsigned>& mr : mixed_to_row) {
    unsigned Pm = (unsigned)P(mr.first) % 4;
    if (Pm % 2 == 1) {
      // Complementation within row of C (including the mixed qubit)
      for (unsigned q2 = mr.first; q2 < n_qbs; ++q2) {
        if (C(mr.second, q2)) {
          for (unsigned q3 = q2 + 1; q3 < n_qbs; ++q3) {
            if (C(mr.second, q3)) {
              E(q2, q3) ^= true;
              E(q3, q2) ^= true;
            }
          }
          // Local phases (we saved to Pm earlier so we can reset P(mr.first)
          // here)
          P(q2) += -Pm;
        }
      }
    } else if (Pm == 2) {
      // Local phases within row of C
      for (unsigned q2 = mr.first; q2 < n_qbs; ++q2) {
        if (C(mr.second, q2)) {
          P(q2) += 2;
        }
      }
    }
    // No global phase change
  }
}

bool APState::operator==(const APState& other) const {
  for (unsigned i = 0; i < P.size(); ++i) {
    if (((unsigned)P(i) % 4) != ((unsigned)other.P(i) % 4)) return false;
  }
  return (A == other.A) && (B == other.B) && (C == other.C) && (E == other.E) &&
         (phase == other.phase);
}

void to_json(nlohmann::json& j, const APState& aps) {
  j["nqubits"] = aps.A.cols();
  j["A"] = aps.A;
  j["B"] = aps.B;
  j["C"] = aps.C;
  j["E"] = aps.E;
  j["P"] = aps.P;
  j["phase"] = aps.phase;
}

void from_json(const nlohmann::json& j, APState& aps) {
  unsigned n_qbs = j.at("nqubits").get<unsigned>();
  MatrixXb A(n_qbs, n_qbs);
  VectorXb B(n_qbs);
  MatrixXb C(n_qbs, n_qbs);
  MatrixXb E(n_qbs, n_qbs);
  Eigen::VectorXi P(n_qbs);
  from_json(j.at("A"), A);
  from_json(j.at("B"), B);
  from_json(j.at("C"), C);
  from_json(j.at("E"), E);
  from_json(j.at("P"), P);
  Expr phase = j.at("phase").get<Expr>();
  aps = APState(A, B, C, E, P, phase);
}

/************************* ChoiAPState IMPLEMENTATION *************************/

static APState id_aps(unsigned n) {
  MatrixXb A(2 * n, 2 * n);
  A << MatrixXb::Identity(n, n), MatrixXb::Identity(n, n),
      MatrixXb::Zero(n, 2 * n);
  return APState(
      A, VectorXb::Zero(2 * n), MatrixXb::Zero(2 * n, 2 * n),
      MatrixXb::Zero(2 * n, 2 * n), Eigen::VectorXi::Zero(2 * n), 0);
}

ChoiAPState::ChoiAPState(unsigned n) : ap_(id_aps(n)), col_index_() {
  for (unsigned i = 0; i < n; ++i) {
    col_index_.insert({{Qubit(i), TableauSegment::Input}, i});
    col_index_.insert({{Qubit(i), TableauSegment::Output}, n + i});
  }
}

ChoiAPState::ChoiAPState(const qubit_vector_t& qbs)
    : ap_(id_aps(qbs.size())), col_index_() {
  unsigned n = qbs.size();
  unsigned i = 0;
  for (const Qubit& qb : qbs) {
    col_index_.insert({{qb, TableauSegment::Input}, i});
    col_index_.insert({{qb, TableauSegment::Output}, n + i});
    ++i;
  }
}

ChoiAPState::ChoiAPState(
    const MatrixXb& A, const VectorXb& B, const MatrixXb& C, const MatrixXb& E,
    const Eigen::VectorXi& P, const Expr& phase, unsigned n_ins)
    : ap_(A, B, C, E, P, phase), col_index_() {
  unsigned n_qbs = A.cols();
  if (n_ins > n_qbs)
    throw std::invalid_argument(
        "Number of inputs of a ChoiAPState cannot be larger than the number of "
        "qubits");
  for (unsigned i = 0; i < n_ins; ++i) {
    col_index_.insert({{Qubit(i), TableauSegment::Input}, i});
  }
  for (unsigned i = 0; i < n_qbs - n_ins; ++i) {
    col_index_.insert({{Qubit(i), TableauSegment::Output}, n_ins + i});
  }
}

unsigned ChoiAPState::get_n_boundaries() const { return ap_.A.cols(); }

unsigned ChoiAPState::get_n_inputs() const { return input_qubits().size(); }

unsigned ChoiAPState::get_n_outputs() const { return output_qubits().size(); }

qubit_vector_t ChoiAPState::input_qubits() const {
  qubit_vector_t ins;
  BOOST_FOREACH (
      tableau_col_index_t::left_const_reference entry, col_index_.left) {
    if (entry.first.second == TableauSegment::Input)
      ins.push_back(entry.first.first);
  }
  return ins;
}

qubit_vector_t ChoiAPState::output_qubits() const {
  qubit_vector_t outs;
  BOOST_FOREACH (
      tableau_col_index_t::left_const_reference entry, col_index_.left) {
    if (entry.first.second == TableauSegment::Output)
      outs.push_back(entry.first.first);
  }
  return outs;
}

void ChoiAPState::apply_gate(
    OpType type, const qubit_vector_t& qbs, TableauSegment seg) {
  std::vector<unsigned> u_qbs;
  for (const Qubit& q : qbs) {
    u_qbs.push_back(col_index_.left.at(col_key_t{q, seg}));
  }
  if (seg == TableauSegment::Output) {
    ap_.apply_gate(type, u_qbs);
  }
  if (seg == TableauSegment::Input) {
    switch (type) {
      case OpType::Z:
      case OpType::X:
      case OpType::S:
      case OpType::Sdg:
      case OpType::SX:
      case OpType::V:
      case OpType::SXdg:
      case OpType::Vdg:
      case OpType::H:
      case OpType::CX:
      case OpType::CZ:
      case OpType::ZZMax:
      case OpType::ISWAPMax:
      case OpType::SWAP:
      case OpType::BRIDGE:
      case OpType::noop: {
        ap_.apply_gate(type, u_qbs);
        break;
      }
      case OpType::Y: {
        ap_.apply_gate(OpType::Y, u_qbs);
        ap_.phase += 1.;
        break;
      }
      case OpType::CY: {
        ap_.apply_V(u_qbs.at(1));
        ap_.apply_X(u_qbs.at(1));
        ap_.apply_CZ(u_qbs.at(0), u_qbs.at(1));
        ap_.apply_V(u_qbs.at(1));
        break;
      }
      case OpType::ECR: {
        ap_.apply_S(u_qbs.at(1));
        ap_.apply_V(u_qbs.at(1));
        ap_.apply_S(u_qbs.at(1));
        ap_.apply_CZ(u_qbs.at(0), u_qbs.at(1));
        ap_.apply_V(u_qbs.at(1));
        ap_.apply_S(u_qbs.at(1));
        ap_.apply_X(u_qbs.at(0));
        ap_.apply_S(u_qbs.at(0));
        ap_.phase += .25;
        break;
      }
      case OpType::Reset: {
        unsigned q = u_qbs.at(0);
        unsigned removed_q = ap_.post_select(q);
        tableau_col_index_t::right_iterator it =
            col_index_.right.find(removed_q);
        col_key_t removed_key = it->second;
        col_index_.right.erase(it);
        it = col_index_.right.find(u_qbs.at(0));
        col_index_.right.erase(it);
        col_index_.insert({removed_key, q});
        unsigned new_q = ap_.init_qubit();
        ap_.apply_V(new_q);
        ap_.collapse_qubit(new_q);
        col_index_.insert({{qbs.at(0), seg}, new_q});
        break;
      }
      case OpType::Phase: {
        throw std::logic_error(
            "OpType::Phase cannot be applied via ChoiAPState::apply_gate");
      }
      default: {
        throw BadOpType(
            "Cannot be applied to a ChoiAPState: not a Clifford gate", type);
      }
    }
  }
}

void ChoiAPState::init_qubit(const Qubit& qb, TableauSegment seg) {
  unsigned u_qb = ap_.init_qubit();
  col_index_.insert({{qb, seg}, u_qb});
}

void ChoiAPState::post_select(const Qubit& qb, TableauSegment seg) {
  tableau_col_index_t::left_iterator lit =
      col_index_.left.find(col_key_t{qb, seg});
  unsigned u_qb = lit->second;
  unsigned removed_u_qb = ap_.post_select(u_qb);
  col_index_.left.erase(lit);
  tableau_col_index_t::right_iterator rit = col_index_.right.find(removed_u_qb);
  col_key_t removed_key = rit->second;
  col_index_.right.erase(rit);
  col_index_.insert({removed_key, u_qb});
}

void ChoiAPState::discard_qubit(const Qubit& qb, TableauSegment seg) {
  collapse_qubit(qb, seg);
  apply_gate(OpType::V, {qb}, seg);
  collapse_qubit(qb, seg);
  post_select(qb, seg);
}

void ChoiAPState::collapse_qubit(const Qubit& qb, TableauSegment seg) {
  unsigned u_qb = col_index_.left.at(col_key_t{qb, seg});
  ap_.collapse_qubit(u_qb);
}

void ChoiAPState::canonical_column_order(TableauSegment first) {
  std::set<Qubit> ins;
  std::set<Qubit> outs;
  BOOST_FOREACH (
      tableau_col_index_t::left_const_reference entry, col_index_.left) {
    if (entry.first.second == TableauSegment::Input)
      ins.insert(entry.first.first);
    else
      outs.insert(entry.first.first);
  }
  tableau_col_index_t new_index;
  unsigned i = 0;
  if (first == TableauSegment::Input) {
    for (const Qubit& q : ins) {
      new_index.insert({{q, TableauSegment::Input}, i});
      ++i;
    }
  }
  for (const Qubit& q : outs) {
    new_index.insert({{q, TableauSegment::Output}, i});
    ++i;
  }
  if (first == TableauSegment::Output) {
    for (const Qubit& q : ins) {
      new_index.insert({{q, TableauSegment::Input}, i});
      ++i;
    }
  }
  MatrixXb A = MatrixXb::Zero(i, i);
  MatrixXb C = MatrixXb::Zero(i, i);
  // E requires reordering both columns and rows
  // Reorder columns to get Etemp, then reorder rows for E
  MatrixXb Etemp = MatrixXb::Zero(i, i);
  MatrixXb E = MatrixXb::Zero(i, i);
  Eigen::VectorXi P = Eigen::VectorXi(i);
  for (unsigned j = 0; j < i; ++j) {
    col_key_t key = new_index.right.at(j);
    unsigned c = col_index_.left.at(key);
    A.col(j) = ap_.A.col(c);
    C.col(j) = ap_.C.col(c);
    Etemp.col(j) = ap_.E.col(c);
    P(j) = ap_.P(c);
  }
  for (unsigned j = 0; j < i; ++j) {
    col_key_t key = new_index.right.at(j);
    unsigned c = col_index_.left.at(key);
    E.row(j) = Etemp.row(c);
  }
  ap_ = APState(A, ap_.B, C, E, P, ap_.phase);
  col_index_ = new_index;
}

void ChoiAPState::normal_form() { ap_.normal_form(); }

void ChoiAPState::rename_qubits(const qubit_map_t& qmap, TableauSegment seg) {
  tableau_col_index_t new_index;
  BOOST_FOREACH (
      tableau_col_index_t::left_const_reference entry, col_index_.left) {
    auto found = qmap.find(entry.first.first);
    if (entry.first.second == seg && found != qmap.end())
      new_index.insert({{found->second, seg}, entry.second});
    else
      new_index.insert({entry.first, entry.second});
  }
  col_index_ = new_index;
}

bool ChoiAPState::operator==(const ChoiAPState& other) const {
  return (col_index_ == other.col_index_) && (ap_ == other.ap_);
}

void to_json(nlohmann::json& j, const ChoiAPState::TableauSegment& seg) {
  j = (seg == ChoiAPState::TableauSegment::Input) ? "In" : "Out";
}

void from_json(const nlohmann::json& j, ChoiAPState::TableauSegment& seg) {
  const std::string str_seg = j.get<std::string>();
  seg = (str_seg == "In") ? ChoiAPState::TableauSegment::Input
                          : ChoiAPState::TableauSegment::Output;
}

void to_json(nlohmann::json& j, const ChoiAPState& cap) {
  j["aps"] = cap.ap_;
  std::vector<ChoiAPState::col_key_t> qbs;
  for (unsigned i = 0; i < cap.get_n_boundaries(); ++i) {
    qbs.push_back(cap.col_index_.right.at(i));
  }
  j["qubits"] = qbs;
}

void from_json(const nlohmann::json& j, ChoiAPState& cap) {
  j.at("aps").get_to(cap.ap_);
  std::vector<ChoiAPState::col_key_t> qbs =
      j.at("qubits").get<std::vector<ChoiAPState::col_key_t>>();
  if ((unsigned)qbs.size() != (unsigned)cap.ap_.A.cols())
    throw std::invalid_argument(
        "Number of qubits in json ChoiAPState does not match APState "
        "size.");
  cap.col_index_.clear();
  for (unsigned i = 0; i < qbs.size(); ++i) {
    cap.col_index_.insert({qbs.at(i), i});
  }
}

}  // namespace tket