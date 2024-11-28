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
  unsigned n_qbs = ap.A.cols();
  Circuit circ(n_qbs);
  circ.qubit_create_all();

  MatrixXb ABgauss(n_qbs, n_qbs + 1);
  ABgauss.block(0, 0, n_qbs, n_qbs) = ap.A;
  ABgauss.col(n_qbs) = ap.B;
  std::vector<std::pair<unsigned, unsigned>> row_ops =
      gaussian_elimination_row_ops(ABgauss);
  for (const std::pair<unsigned, unsigned>& op : row_ops) {
    for (unsigned j = 0; j < n_qbs + 1; ++j) {
      ABgauss(op.second, j) ^= ABgauss(op.first, j);
    }
  }
  std::map<unsigned, unsigned> leader_to_row;
  for (unsigned r = 0; r < n_qbs; ++r) {
    bool empty_row = true;
    for (unsigned c = 0; c < n_qbs; ++c) {
      if (ABgauss(r, c)) {
        leader_to_row.insert({c, r});
        empty_row = false;
        break;
      }
    }
    if (empty_row) break;
  }

  for (unsigned q = 0; q < n_qbs; ++q) {
    auto lr = leader_to_row.find(q);
    if (lr == leader_to_row.end())
      circ.add_op<unsigned>(OpType::H, {q});
    else if (ABgauss(lr->second, n_qbs))
      circ.add_op<unsigned>(OpType::X, {q});
  }
  for (const std::pair<const unsigned, unsigned>& lr : leader_to_row) {
    for (unsigned ctrl = lr.first + 1; ctrl < n_qbs; ++ctrl) {
      if (ABgauss(lr.second, ctrl))
        circ.add_op<unsigned>(OpType::CX, {ctrl, lr.first});
    }
  }
  for (unsigned d = 0; d < n_qbs; ++d) {
    std::optional<unsigned> first = std::nullopt;
    for (unsigned c = 0; c < n_qbs; ++c) {
      if (ap.C(d, c)) {
        if (first) {
          circ.add_op<unsigned>(OpType::CX, {c, *first});
        } else {
          first = c;
        }
      }
    }
    if (first) {
      circ.add_op<unsigned>(OpType::Collapse, {*first});
      for (unsigned c = n_qbs; c > *first + 1;) {
        --c;
        if (ap.C(d, c)) circ.add_op<unsigned>(OpType::CX, {c, *first});
      }
    }
  }
  for (unsigned q1 = 0; q1 < n_qbs; ++q1) {
    for (unsigned q2 = q1 + 1; q2 < n_qbs; ++q2) {
      if (ap.E(q1, q2)) circ.add_op<unsigned>(OpType::CZ, {q1, q2});
    }
    switch ((unsigned)ap.P(q1) % 4) {
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

ChoiAPState circuit_to_choi_apstate(const Circuit& circ) {
  ChoiAPState ap(circ.all_qubits());
  for (const Qubit& q : circ.created_qubits()) {
    ap.post_select(q, ChoiAPState::TableauSegment::Input);
  }
  for (const Command& com : circ) {
    auto args = com.get_args();
    qubit_vector_t qbs = {args.begin(), args.end()};
    ap.apply_gate(com.get_op_ptr()->get_type(), qbs);
  }
  ap.rename_qubits(
      circ.implicit_qubit_permutation(), ChoiAPState::TableauSegment::Output);
  for (const Qubit& q : circ.discarded_qubits()) {
    ap.discard_qubit(q);
  }
  ap.canonical_column_order();
  ap.normal_form();
  return ap;
}

std::pair<Circuit, qubit_map_t> choi_apstate_to_exact_circuit(
    ChoiAPState ap, CXConfigType cx_config) {
  ChoiMixTableau tab = choi_apstate_to_cm_tableau(ap);
  std::pair<Circuit, qubit_map_t> res =
      cm_tableau_to_exact_circuit(tab, cx_config);
  ChoiAPState reconstructed = circuit_to_choi_apstate(res.first);
  reconstructed.rename_qubits(res.second, ChoiAPState::TableauSegment::Output);
  ap.canonical_column_order();
  ap.normal_form();
  reconstructed.canonical_column_order();
  reconstructed.normal_form();
  res.first.add_phase(ap.ap_.phase - reconstructed.ap_.phase);
  return res;
}

std::pair<Circuit, qubit_map_t> choi_apstate_to_unitary_extension_circuit(
    ChoiAPState ap, const std::vector<Qubit>& init_names,
    const std::vector<Qubit>& post_names, CXConfigType cx_config) {
  ChoiMixTableau tab = choi_apstate_to_cm_tableau(ap);
  std::pair<Circuit, qubit_map_t> res = cm_tableau_to_unitary_extension_circuit(
      tab, init_names, post_names, cx_config);
  ChoiAPState reconstructed = circuit_to_choi_apstate(res.first);
  reconstructed.rename_qubits(res.second, ChoiAPState::TableauSegment::Output);
  ap.canonical_column_order();
  ap.normal_form();
  reconstructed.canonical_column_order();
  reconstructed.normal_form();
  res.first.add_phase(ap.ap_.phase - reconstructed.ap_.phase);
  return res;
}

APState tableau_to_apstate(SymplecticTableau tab) {
  unsigned n_qbs = tab.get_n_qubits();
  unsigned n_rows = tab.get_n_rows();
  if (n_rows > n_qbs)
    throw std::logic_error(
        "tableau_to_apstate requires a tableau with up to n commuting rows for "
        "n qubits");
  MatrixXb fullmat = MatrixXb::Zero(n_rows, 2 * n_qbs);
  /**
   * Gaussian elimination by the x matrix first ensures the bottom rows are only
   * Zs, i.e. describing rows of A. Reversing the columns of the x matrix
   * guarantees that each row has an X on at most one free qubit, simplifying
   * the code for finding E and P
   */
  for (unsigned c = 0; c < n_qbs; ++c) {
    fullmat.col(c) = tab.xmat.col(n_qbs - 1 - c);
  }
  fullmat.block(0, n_qbs, n_rows, n_qbs) = tab.zmat;
  std::vector<std::pair<unsigned, unsigned>> row_ops =
      gaussian_elimination_row_ops(fullmat);
  for (const std::pair<unsigned, unsigned>& op : row_ops)
    tab.row_mult(op.first, op.second);

  MatrixXb A = MatrixXb::Zero(n_qbs, n_qbs);
  VectorXb B = VectorXb::Zero(n_qbs);
  MatrixXb C = MatrixXb::Zero(n_qbs, n_qbs);
  MatrixXb E = MatrixXb::Zero(n_qbs, n_qbs);
  Eigen::VectorXi P = Eigen::VectorXi::Zero(n_qbs);

  unsigned n_leading = 0;
  for (unsigned r = n_rows; r > 0;) {
    --r;
    bool is_a_row = true;
    for (unsigned c = 0; c < n_qbs; ++c) {
      if (tab.xmat(r, c)) {
        is_a_row = false;
        break;
      }
    }
    if (!is_a_row) break;
    ++n_leading;
  }
  A.block(0, 0, n_leading, n_qbs) =
      tab.zmat.block(n_rows - n_leading, 0, n_leading, n_qbs);
  B.head(n_leading) = tab.phase.tail(n_leading);
  std::map<unsigned, unsigned> leader_to_row;
  for (unsigned r = 0; r < n_leading; ++r) {
    for (unsigned c = r; c < n_qbs; ++c) {
      if (A(r, c)) {
        leader_to_row.insert({c, r});
        break;
      }
    }
  }

  /**
   * Each free qubit q is after all leaders connected to it, by reduced
   * row-echelon of A. Therefore Gaussian elimination of the x matrix in reverse
   * order gives rows corresponding to the columns of A of free qubits (plus the
   * free qubit itself).
   *
   * Then the corresponding row q of the z matrix is the free qubit's row of E,
   * plus the rows of any E for any mixed qubit connected to q in C, plus an
   * extra flip at q if P(q) is odd. We can first look at the column for each
   * mixed qubit to identify their rows of E and subtract from this matrix to
   * leave the interactions between free qubits and their local phases.
   */

  // Start by identifying the free and mixed qubits
  std::map<unsigned, unsigned> free_to_row;
  for (unsigned r = 0; r < n_qbs - n_leading; ++r) {
    for (unsigned c = n_qbs - r; c > 0;) {
      --c;
      if (tab.xmat(r, c)) {
        free_to_row.insert({c, r});
        break;
      }
    }
  }
  std::map<unsigned, unsigned> mixed_to_row;
  for (unsigned q = 0; q < n_qbs; ++q) {
    if (leader_to_row.find(q) == leader_to_row.end() &&
        free_to_row.find(q) == free_to_row.end()) {
      unsigned r = mixed_to_row.size();
      mixed_to_row.insert({q, r});
      C(r, q) = true;
    }
  }
  // Identify C and the mixed rows of E by looking at what mixed qubits appear
  // in the rows of each free qubit
  for (const std::pair<const unsigned, unsigned>& fr : free_to_row) {
    for (const std::pair<const unsigned, unsigned>& mr : mixed_to_row) {
      if (tab.xmat(fr.second, mr.first)) C(mr.second, fr.first) = true;
      if (tab.zmat(fr.second, mr.first)) {
        E(mr.first, fr.first) = true;
        E(fr.first, mr.first) = true;
      }
    }
  }
  // Identify connections in E between free qubits
  for (auto fr1 = free_to_row.begin(); fr1 != free_to_row.end(); ++fr1) {
    for (auto fr2 = free_to_row.begin(); fr2 != fr1; ++fr2) {
      unsigned n_shared_mixed = 0;
      for (const std::pair<const unsigned, unsigned>& mr : mixed_to_row) {
        if (C(mr.second, fr1->first) && C(mr.second, fr2->first))
          ++n_shared_mixed;
      }
      if (tab.zmat(fr1->second, fr2->first) ^ (n_shared_mixed % 2 == 1)) {
        E(fr1->first, fr2->first) = true;
        E(fr2->first, fr1->first) = true;
      }
    }
  }
  // Identify P
  for (const std::pair<const unsigned, unsigned>& fr : free_to_row) {
    unsigned n_mixed_in_c_and_e = 0;
    for (const std::pair<const unsigned, unsigned>& mr : mixed_to_row) {
      if (C(mr.second, fr.first) && E(mr.first, fr.first)) ++n_mixed_in_c_and_e;
    }
    if (tab.zmat(fr.second, fr.first) ^ (n_mixed_in_c_and_e % 2 == 1))
      P(fr.first) += 1;
    if (tab.phase(fr.second) ^ (n_mixed_in_c_and_e % 2 == 1)) P(fr.first) += 2;
  }

  return APState(A, B, C, E, P, 0);
}

SymplecticTableau apstate_to_tableau(APState ap) {
  unsigned n_qbs = ap.A.cols();
  // Want A and C in reduced row-echelon form to identify leaders and mixed
  // qubits, but don't need the rest in normal form
  std::vector<std::pair<unsigned, unsigned>> row_ops =
      gaussian_elimination_row_ops(ap.A);
  for (const std::pair<unsigned, unsigned>& op : row_ops) {
    for (unsigned j = 0; j < n_qbs; ++j) {
      ap.A(op.second, j) ^= ap.A(op.first, j);
    }
    ap.B(op.second) ^= ap.B(op.first);
  }
  std::map<unsigned, unsigned> leader_to_row;
  for (unsigned r = 0; r < n_qbs; ++r) {
    bool leader_found = false;
    for (unsigned c = 0; c < n_qbs; ++c) {
      if (ap.A(r, c)) {
        leader_found = true;
        leader_to_row.insert({c, r});
        break;
      }
    }
    if (!leader_found) break;
  }
  for (unsigned r = 0; r < n_qbs; ++r) {
    for (const std::pair<const unsigned, unsigned>& lr : leader_to_row) {
      if (ap.C(r, lr.first)) {
        for (unsigned c = 0; c < n_qbs; ++c) {
          ap.C(r, c) ^= ap.A(lr.second, c);
        }
      }
    }
  }
  row_ops = gaussian_elimination_row_ops(ap.C);
  for (const std::pair<unsigned, unsigned>& op : row_ops) {
    for (unsigned j = 0; j < n_qbs; ++j) {
      ap.C(op.second, j) ^= ap.C(op.first, j);
    }
  }
  std::map<unsigned, unsigned> mixed_to_row;
  for (unsigned r = 0; r < n_qbs; ++r) {
    bool mixed_found = false;
    for (unsigned c = 0; c < n_qbs; ++c) {
      if (ap.C(r, c)) {
        mixed_found = true;
        mixed_to_row.insert({c, r});
        break;
      }
    }
    if (!mixed_found) break;
  }
  unsigned n_rows = n_qbs - mixed_to_row.size();
  MatrixXb xmat = MatrixXb::Zero(n_rows, n_qbs);
  MatrixXb zmat = MatrixXb::Zero(n_rows, n_qbs);
  VectorXb phase = VectorXb::Zero(n_rows);

  // One stabiliser per leader, with Z on that qubit and Z on every neighbour in
  // A
  unsigned n_leading = leader_to_row.size();
  zmat.block(0, 0, n_leading, n_qbs) = ap.A.block(0, 0, n_leading, n_qbs);
  phase.head(n_leading) = ap.B.head(n_leading);

  // One stabiliser per free, with X on that qubit, every neighbour in C, and
  // the odd neighbourhood of this set in A; pushing this through the phase
  // polynomial adds Zs
  unsigned r = n_leading;
  for (unsigned q = 0; q < n_qbs; ++q) {
    if (leader_to_row.find(q) == leader_to_row.end() &&
        mixed_to_row.find(q) == mixed_to_row.end()) {
      // Calculate the Xs
      xmat(r, q) = true;
      for (const std::pair<const unsigned, unsigned>& mr : mixed_to_row) {
        xmat(r, mr.first) = ap.C(mr.second, q);
      }
      for (const std::pair<const unsigned, unsigned>& lr : leader_to_row) {
        unsigned n_mixed_in_between = 0;
        for (const std::pair<const unsigned, unsigned>& mr : mixed_to_row) {
          if (ap.A(lr.second, mr.first) && ap.C(mr.second, q))
            ++n_mixed_in_between;
        }
        xmat(r, lr.first) = ap.A(lr.second, q) ^ (n_mixed_in_between % 2 == 1);
      }
      // Push through the phase polynomial to calculate the Zs
      for (unsigned q2 = 0; q2 < n_qbs; ++q2) {
        if (xmat(r, q2)) {
          // Pushing an X through each CZ creates a Z
          for (unsigned q3 = 0; q3 < n_qbs; ++q3) {
            if (ap.E(q2, q3)) {
              zmat(r, q3) ^= true;
              // If both qubits of a CZ have Xs, we will need to reorder the X
              // and Z on one to match the other
              if (q2 < q3 && xmat(r, q3)) phase(r) ^= true;
            }
          }
          // Pushing an X through the local phase
          zmat(r, q2) ^= ((unsigned)ap.P(q2) % 2 == 1);
          phase(r) ^= ((unsigned)ap.P(q2) % 4 > 1);
        }
      }
      ++r;
    }
  }

  return SymplecticTableau(xmat, zmat, phase);
}

ChoiAPState cm_tableau_to_choi_apstate(const ChoiMixTableau& tab) {
  ChoiAPState ap(0);
  ap.ap_ = tableau_to_apstate(tab.tab_);
  BOOST_FOREACH (
      ChoiMixTableau::tableau_col_index_t::left_const_reference entry,
      tab.col_index_.left) {
    ChoiAPState::TableauSegment seg =
        (entry.first.second == ChoiMixTableau::TableauSegment::Input)
            ? ChoiAPState::TableauSegment::Input
            : ChoiAPState::TableauSegment::Output;
    ap.col_index_.insert({{entry.first.first, seg}, entry.second});
  }
  return ap;
}

ChoiMixTableau choi_apstate_to_cm_tableau(const ChoiAPState& ap) {
  ChoiMixTableau tab(0);
  tab.tab_ = apstate_to_tableau(ap.ap_);
  BOOST_FOREACH (
      ChoiAPState::tableau_col_index_t::left_const_reference entry,
      ap.col_index_.left) {
    ChoiMixTableau::TableauSegment seg =
        (entry.first.second == ChoiAPState::TableauSegment::Input)
            ? ChoiMixTableau::TableauSegment::Input
            : ChoiMixTableau::TableauSegment::Output;
    tab.col_index_.insert({{entry.first.first, seg}, entry.second});
  }
  return tab;
}

}  // namespace tket