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

#include "tket/Clifford/APState.hpp"

#include <tkassert/Assert.hpp>

#include "tket/OpType/OpTypeInfo.hpp"

namespace tket {

APState::APState(
    const MatrixXb& A_, const VectorXb& B_, const MatrixXb& E_,
    const Eigen::VectorXi& P_, const Expr& phase_)
    : A(A_), B(B_), E(E_), P(P_), phase(phase_) {
  verify();
}

APState::APState(unsigned n_qubits)
    : A(MatrixXb::Identity(n_qubits, n_qubits)),
      B(VectorXb::Zero(n_qubits)),
      E(MatrixXb::Zero(n_qubits, n_qubits)),
      P(Eigen::VectorXi(n_qubits)),
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
      if (n_non_zero == (1 << offsets.size())) offsets.push_back(x);
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
          ((clifford_phase(complex_phase_cz) - P(qfree) - P(qfree2)) % 4) == 2;
    }
  }

  phase = Expr(std::arg(sv(neutral_z0)) / PI);
}

void APState::verify() {
  // Check dimensions all agree
  unsigned n_qubits = A.rows();
  TKET_ASSERT(A.cols() == n_qubits);
  TKET_ASSERT(B.size() == n_qubits);
  TKET_ASSERT(E.rows() == n_qubits);
  TKET_ASSERT(E.cols() == n_qubits);
  TKET_ASSERT(P.size() == n_qubits);

  for (unsigned r = 0; r < n_qubits; ++r) {
    for (unsigned c = 0; c < r; ++c) {
      // Check A is upper triangular
      TKET_ASSERT(!A(r, c));
      // Check E is symmetric
      TKET_ASSERT(E(r, c) == E(c, r));
    }
    if (A(r, r)) {
      // Check leaders are eliminated in A (row echelon)
      for (unsigned r2 = 0; r2 < r; ++r2) {
        TKET_ASSERT(!A(r2, r));
      }
    } else {
      // Check A only has entries in rows with a leader
      for (unsigned c = r + 1; c < n_qubits; ++c) {
        TKET_ASSERT(!A(r, c));
      }
    }
  }

  for (unsigned q = 0; q < n_qubits; ++q) {
    if (A(q, q)) {
      // Check P is zero on leaders
      TKET_ASSERT(P(q) % 4 == 0);
      // Check E is zero on leaders (simplified by symmetry)
      for (unsigned q2 = 0; q2 < n_qubits; ++q2) {
        TKET_ASSERT(!E(q, q2));
      }
    } else {
      // Check B is zero on free qubits
      TKET_ASSERT(!B(q));
      // Check diagonals of E are still zero
      TKET_ASSERT(!E(q, q));
    }
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

Eigen::VectorXcd APState::to_statevector() {
  unsigned n_qubits = A.cols();
  Eigen::VectorXcd sv = Eigen::VectorXcd::Zero(1 << n_qubits);
  unsigned n_terms = 0;
  Complex g_phase = std::exp(i_ * PI * eval_expr(phase).value());
  for (unsigned x = 0; x < 1 << n_qubits; ++x) {
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

void APState::apply_CZ(unsigned ctrl, unsigned trgt) {
  if (ctrl == trgt)
    throw std::logic_error(
        "APState::apply_CZ: cannot apply CZ when control and target qubits "
        "match");
  unsigned n_qubits = A.cols();
  if (ctrl >= n_qubits || trgt >= n_qubits)
    throw std::logic_error("APState::apply_CZ: qubit indices out of range");

  if (A(ctrl, ctrl)) {
    if (A(trgt, trgt)) {
      // ctrl and trgt are both leading

      std::set<unsigned> just_ctrl, just_trgt, both;
      // Add phases and determine connectivity
      for (unsigned q = 0; q < n_qubits; ++q) {
        if (A(ctrl, q)) {
          if (B(trgt)) P(q) += 2;
          if (A(trgt, q)) {
            if (!B(ctrl)) P(q) += 2;
            both.insert(q);
          } else {
            just_ctrl.insert(q);
          }
        } else if (A(trgt, q)) {
          if (B(ctrl)) P(q) += 2;
          just_trgt.insert(q);
        }
      }
      // ctrl and trgt were included in the previous loop, so reset them
      just_ctrl.erase(ctrl);
      just_trgt.erase(trgt);
      P(ctrl) = 0;
      P(trgt) = 0;
      // Pivoting-style update between those connected to ctrl, trgt, and both
      for (const unsigned q1 : just_ctrl) {
        for (const unsigned q2 : just_trgt) {
          E(q1, q2) ^= true;
          E(q2, q1) ^= true;
        }
        for (const unsigned q2 : both) {
          E(q1, q2) ^= true;
          E(q2, q1) ^= true;
        }
      }
      for (const unsigned q1 : just_trgt) {
        for (const unsigned q2 : both) {
          E(q1, q2) ^= true;
          E(q2, q1) ^= true;
        }
      }

      if (B(ctrl) && B(trgt)) phase += 1;
    } else {
      // ctrl is leading, trgt is free

      // Flip connectivity between trgt and every q on ctrl
      for (unsigned q = ctrl + 1; q < n_qubits; ++q) {
        if (A(ctrl, q)) {
          E(trgt, q) ^= true;
          E(q, trgt) ^= true;
          // Note that if q = trgt, this does nothing, preserving the zero
          // diagonal of E
        }
      }
      // Add phase to trgt
      if (A(ctrl, trgt) ^ B(ctrl)) P(trgt) += 2;

      // No global phase change
    }
  } else {
    if (A(trgt, trgt)) {
      // trgt is leading, ctrl is free

      // CZ is symmetric, so reuse the previous case
      apply_CZ(trgt, ctrl);
    } else {
      // ctrl and trgt are both free

      // Just add the CZ gate to E
      E(ctrl, trgt) ^= true;
      E(trgt, ctrl) ^= true;

      // No global phase change
    }
  }
}

void APState::apply_S(unsigned q) {
  unsigned n_qubits = A.cols();
  if (q >= n_qubits)
    throw std::logic_error("APState::apply_S: qubit index out of range");

  if (A(q, q)) {
    // q is leading

    // Local complementation update to E
    for (unsigned q1 = q + 1; q1 < n_qubits; ++q1) {
      if (A(q, q1)) {
        for (unsigned q2 = q + 1; q2 < n_qubits; ++q2) {
          E(q1, q2) ^= ((q1 != q2) && A(q, q2));
        }
      }
    }

    // Update global phase
    if (B(q)) phase += .5;

    // Update local phases
    unsigned local_phase_change = B(q) ? 3 : 1;
    for (unsigned q1 = q + 1; q1 < n_qubits; ++q1) {
      if (A(q, q1)) P(q1) += local_phase_change;
    }
  } else {
    // q is free

    // Add to local phase
    P(q) += 1;

    // No global phase change
  }
}

void APState::apply_V(unsigned q) {
  unsigned n_qubits = A.cols();
  if (q >= n_qubits)
    throw std::logic_error("APState::apply_V: qubit index out of range");

  if (A(q, q)) {
    // q is leading, but will become free

    // Local complementation for E
    for (unsigned q1 = q; q1 < n_qubits; ++q1) {
      if (A(q, q1)) {
        for (unsigned q2 = q1 + 1; q2 < n_qubits; ++q2) {
          E(q1, q2) ^= A(q, q2);
          E(q2, q1) ^= A(q, q2);
        }

        // Local phases
        P(q1) += B(q) ? 1 : 3;
      }
    }

    // Global phase update
    if (B(q)) phase += -.5;

    // q is no longer a leader
    for (unsigned q1 = q; q1 < n_qubits; ++q1) {
      A(q, q1) = false;
    }
    B(q) = false;
  } else {
    // q is free

    std::optional<unsigned> last_connected_leader;
    for (unsigned q1 = q; q1 > 0;) {
      --q1;
      if (A(q1, q)) {
        last_connected_leader = q1;
        break;
      }
    }

    if (last_connected_leader.has_value()) {
      // q is free and connected to a leader

      // For each other leader connected to q, add LCL's row in A
      std::list<unsigned> other_leaders;
      for (unsigned q1 = 0; q1 < *last_connected_leader; ++q1) {
        if (A(q1, q)) {
          other_leaders.push_back(q1);
          for (unsigned q2 = *last_connected_leader; q2 < n_qubits; ++q2) {
            A(q1, q2) ^= A(*last_connected_leader, q2);
          }
        }
      }

      // Sort the connected frees into just q, just LCL, and both
      std::list<unsigned> just_q, just_lcl, both;
      for (unsigned q1 = 0; q1 < n_qubits; ++q1) {
        if (q1 == q || q1 == *last_connected_leader) continue;
        if (E(q, q1)) {
          if (A(*last_connected_leader, q1))
            both.push_back(q1);
          else
            just_q.push_back(q1);
        } else if (A(*last_connected_leader, q1))
          just_lcl.push_back(q1);
      }

      // Split case by phase on q
      // In each case:
      // - Perform appropriate complementations between connected frees
      // - Set LCL's connectivity in E appropriately
      // - Modify q's connectivity in E
      // - Add all local and global phases by a case matrix
      if (P(q) % 2 == 0) {
        // Handle complementations between neighbours
        for (const unsigned q1 : just_lcl) {
          // Local complement within group
          for (const unsigned q2 : just_lcl) {
            E(q1, q2) ^= true;
          }
          E(q1, q1) ^= true;
          // Complement between groups
          for (const unsigned q2 : just_q) {
            E(q1, q2) ^= true;
            E(q2, q1) ^= true;
          }
          // Add connections with q
          E(q, q1) ^= true;
          E(q1, q) ^= true;
        }
        for (const unsigned q1 : both) {
          // Local complement within group
          for (const unsigned q2 : both) {
            E(q1, q2) ^= true;
          }
          E(q1, q1) ^= true;
          // Complement between groups
          for (const unsigned q2 : just_q) {
            E(q1, q2) ^= true;
            E(q2, q1) ^= true;
          }
        }
        // Remove connections with q
        for (const unsigned q1 : just_q) {
          E(q, q1) ^= true;
          E(q1, q) ^= true;
        }

        // Move LCL from A to E
        E(q, *last_connected_leader) = true;
        E(*last_connected_leader, q) = true;
        for (const unsigned q1 : just_q) {
          E(q1, *last_connected_leader) = true;
          E(*last_connected_leader, q1) = true;
        }
        for (const unsigned q1 : just_lcl) {
          E(q1, *last_connected_leader) = true;
          E(*last_connected_leader, q1) = true;
          A(*last_connected_leader, q1) = false;
        }
        for (const unsigned q1 : both) A(*last_connected_leader, q1) = false;
        A(*last_connected_leader, q) = false;
        A(*last_connected_leader, *last_connected_leader) = false;

        // Phases
        bool a = (P(q) % 4 == 2);
        bool b = B(*last_connected_leader);
        P(q) = b ? 1 : 3;
        P(*last_connected_leader) = (a ^ b) ? 1 : 3;
        B(*last_connected_leader) = false;
        for (const unsigned q1 : other_leaders) B(q1) ^= b;
        for (const unsigned q1 : both) P(q1) += a ? 3 : 1;
        for (const unsigned q1 : just_lcl) P(q1) += (b ^ a) ? 1 : 3;
        for (const unsigned q1 : just_q) P(q1) += b ? 2 : 0;
        if (b) phase += a ? .5 : -.5;
      } else {
        // Handle complementations between neighbours
        for (const unsigned q1 : just_lcl) {
          // Complement between groups
          for (const unsigned q2 : just_q) {
            E(q1, q2) ^= true;
            E(q2, q1) ^= true;
          }
          for (const unsigned q2 : both) {
            E(q1, q2) ^= true;
            E(q2, q1) ^= true;
          }
          // Add connections with q
          E(q, q1) = true;
          E(q1, q) = true;
        }
        for (const unsigned q1 : just_q) {
          // Complement between groups
          for (const unsigned q2 : both) {
            E(q1, q2) ^= true;
            E(q2, q1) ^= true;
          }
          // Remove connections with q
          E(q, q1) = false;
          E(q1, q) = false;
        }

        // Move LCL from A to E
        E(q, *last_connected_leader) = true;
        E(*last_connected_leader, q) = true;
        for (const unsigned q1 : just_q) {
          E(q1, *last_connected_leader) = true;
          E(*last_connected_leader, q1) = true;
        }
        for (const unsigned q1 : both) {
          E(q1, *last_connected_leader) = true;
          E(*last_connected_leader, q1) = true;
          A(*last_connected_leader, q1) = false;
        }
        for (const unsigned q1 : just_lcl)
          A(*last_connected_leader, q1) = false;
        A(*last_connected_leader, q) = false;
        A(*last_connected_leader, *last_connected_leader) = false;

        // Phases
        bool a = (P(q) % 4 == 3);
        bool b = B(*last_connected_leader);
        P(q) = b ? 1 : 3;
        P(*last_connected_leader) = a ? 2 : 0;
        B(*last_connected_leader) = false;
        for (const unsigned q1 : other_leaders) B(q1) ^= b;
        for (const unsigned q1 : just_q) P(q1) += b ? 2 : 0;
        for (const unsigned q1 : just_lcl) P(q1) += a ? 2 : 0;
        for (const unsigned q1 : both) P(q1) += (a ^ b) ? 0 : 2;
        if (a && b) phase += 1.;
      }
    } else if (P(q) % 2 == 0) {
      // q has phase 0/pi and no connected leader

      // Local complementation update to E
      std::list<unsigned> connected;
      for (unsigned q1 = 0; q1 < n_qubits; ++q1) {
        if (E(q, q1)) {
          for (const unsigned q2 : connected) {
            E(q1, q2) ^= true;
            E(q2, q1) ^= true;
          }
          connected.push_back(q1);

          // Add local phase
          P(q1) += (P(q) % 4 == 0) ? 1 : 3;
        }
      }

      // Local phase on q remains the same

      // Global phase
      phase += (P(q) % 4 == 0) ? -0.25 : 0.25;
    } else {
      // q has phase +-pi/2 and no connected leader

      std::optional<unsigned> first_connected_free = std::nullopt;
      for (unsigned q1 = 0; q1 < q; ++q1) {
        if (E(q, q1)) {
          first_connected_free = q1;
          break;
        }
      }

      if (first_connected_free.has_value()) {
        // q has phase +-pi/2, no connected leader, but a previous free which
        // will become the new leader

        // Make FCF a leader
        for (unsigned q1 = *first_connected_free; q1 < n_qubits; ++q1) {
          if (E(q, q1)) {
            A(*first_connected_free, q1) = true;
            // q ends up disconnected in E
            E(q, q1) = false;
            E(q1, q) = false;
          }
        }
        A(*first_connected_free, q) = true;
        B(*first_connected_free) = (P(q) % 4 == 3);

        // Set phase of q to -pi/2
        P(q) = 3;

        // No global phase change

        // Make A row reduced
        for (unsigned q1 = 0; q1 < *first_connected_free; ++q1) {
          if (A(q1, *first_connected_free)) {
            for (unsigned q2 = *first_connected_free; q2 < n_qubits; ++q2) {
              A(q1, q2) ^= A(*first_connected_free, q2);
            }
            B(q1) ^= B(*first_connected_free);
          }
        }

        // Need to apply rules for S and CZ gates on FCF to remove
        unsigned s_gates = P(*first_connected_free);
        P(*first_connected_free) = 0;
        std::list<unsigned> cz_targets;
        E(*first_connected_free, q) = false;
        E(q, *first_connected_free) = false;
        for (unsigned q1 = 0; q1 < n_qubits; ++q1) {
          if (E(*first_connected_free, q1)) {
            cz_targets.push_back(q1);
            E(*first_connected_free, q1) = false;
            E(q1, *first_connected_free) = false;
          }
        }
        for (unsigned i = 0; i < s_gates; ++i) apply_S(*first_connected_free);
        for (const unsigned trgt : cz_targets)
          apply_CZ(*first_connected_free, trgt);
      } else {
        // q has phase +-pi/2, no connected leader, and no previous free, so q
        // will become the new leader
        A(q, q) = true;

        std::list<unsigned> connected;
        for (unsigned q1 = q + 1; q1 < n_qubits; ++q1) {
          if (E(q, q1)) {
            // Local complementation
            for (const unsigned q2 : connected) {
              E(q1, q2) ^= true;
              E(q2, q1) ^= true;
            }
            connected.push_back(q1);

            // Add to A
            A(q, q1) = true;

            // Reset E rows
            E(q, q1) = false;
            E(q1, q) = false;
          }
        }

        // Update local and global phases
        if (P(q) % 4 == 1) {
          for (const unsigned q1 : connected) P(q1) += 3;
        } else {
          for (const unsigned q1 : connected) P(q1) += 1;
          B(q) = true;
          phase += -.5;
        }
        P(q) = 0;
      }
    }
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
      apply_V(qbs.at(0));
      apply_V(qbs.at(0));
      phase += .5;
      break;
    }
    case OpType::Y: {
      apply_S(qbs.at(0));
      apply_S(qbs.at(0));
      apply_V(qbs.at(0));
      apply_V(qbs.at(0));
      phase += 1.;
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
      apply_V(qbs.at(0));
      apply_V(qbs.at(0));
      phase += 1.;
      break;
    }
    case OpType::SXdg: {
      apply_V(qbs.at(0));
      apply_V(qbs.at(0));
      apply_V(qbs.at(0));
      phase += .75;
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
      apply_V(qbs.at(1));
      apply_V(qbs.at(1));
      phase += 1.;
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
      apply_V(qbs.at(0));
      apply_V(qbs.at(0));
      apply_S(qbs.at(1));
      apply_V(qbs.at(1));
      apply_CZ(qbs.at(0), qbs.at(1));
      apply_S(qbs.at(1));
      apply_V(qbs.at(1));
      apply_S(qbs.at(1));
      phase += .25;
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

}  // namespace tket