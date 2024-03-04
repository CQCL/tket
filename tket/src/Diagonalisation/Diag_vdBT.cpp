#include "tket/Clifford/ChoiMixTableau.hpp"
#include "tket/Diagonalisation/Diagonalisation.hpp"

namespace tket {

unsigned diagonalise_X_block(ChoiMixTableau &tab, Circuit &cliff_circ, const MatrixXb &xmat, const MatrixXb &zmat) {
  // Algorithm 1: Diagonalise the X block
  unsigned k = 0;
  while (true) {
    bool found = false;
    for (unsigned i=k; i<xmat.rows(); ++i) {
      for (unsigned j=k; j<xmat.cols(); ++j) {
        if (xmat(i, j)) {
          if (i != k) {
            // Swap rows i and k
            tab.tab_.row_mult(i, k);
            tab.tab_.row_mult(k, i);
            tab.tab_.row_mult(i, k);
          }
          if (j != k) {
            // Swap columns j and k
            tab.tab_.apply_gate(OpType::SWAP, {j, k});
            // Update column index to make this into a wireswap
            Qubit j_qb = tab.col_index_.right.at(j).first;
            Qubit k_qb = tab.col_index_.right.at(k).first;
            tab.rename_qubits({{j_qb, k_qb}, {k_qb, j_qb}}, ChoiMixTableau::TableauSegment::Output);
          }
          for (unsigned i2=0; i2<xmat.rows(); ++i2) {
            if (i2 != k && xmat(i2, k)) {
              // Sweep row i2 with row k
              tab.tab_.row_mult(k, i2);
            }
          }
          found = true;
          break;
        }
      }
      if (found) break;
    }
    if (found) ++k;
    else break;
  }
  unsigned k2 = k;
  while (true) {
    bool found = false;
    for (unsigned i=k; i<xmat.rows(); ++i) {
      for (unsigned j=k; j<xmat.cols(); ++j) {
        if (zmat(i, j)) {
          if (i != k) {
            // Swap rows i and k
            tab.tab_.row_mult(i, k);
            tab.tab_.row_mult(k, i);
            tab.tab_.row_mult(i, k);
          }
          if (j != k) {
            // Swap columns j and k
            tab.tab_.apply_gate(OpType::SWAP, {j, k});
            // Update column index to make this into a wireswap
            Qubit j_qb = tab.col_index_.right.at(j).first;
            Qubit k_qb = tab.col_index_.right.at(k).first;
            tab.rename_qubits({{j_qb, k_qb}, {k_qb, j_qb}}, ChoiMixTableau::TableauSegment::Output);
          }
          for (unsigned i2=0; i2<xmat.rows(); ++i2) {
            if (i2 != k && xmat(i2, k)) {
              // Sweep row i2 with row k
              tab.tab_.row_mult(k, i2);
            }
          }
          found = true;
        }
      }
      if (found) break;
    }
    if (found) ++k;
    else break;
  }
  for (unsigned j=k2; j<k; ++j) {
    tab.tab_.apply_gate(OpType::H, {j});
    Qubit j_qb = tab.col_index_.right.at(j).first;
    cliff_circ.add_op<Qubit>(OpType::H, {j_qb});
  }
  for (unsigned i=0; i<k; ++i) {
    for (unsigned j=k; j<xmat.cols(); ++j) {
      if (xmat(i, j)) {
        tab.tab_.apply_gate(OpType::CX, {i,j});
        Qubit i_qb = tab.col_index_.right.at(i).first;
        Qubit j_qb = tab.col_index_.right.at(j).first;
        cliff_circ.add_op<Qubit>(OpType::CX, {i_qb,j_qb});
      }
    }
  }
  
  // Return rank of X matrix
  return k;
}

Circuit mutual_diagonalise_vdBT_PE(
    std::list<SpSymPauliTensor> &gadgets, std::set<Qubit> qubits,
    CXConfigType) {
  // See Section 4 of van den Berg & Temme, "Circuit optimization of Hamiltonian simulation by simultaneous diagonalization of Pauli clusters", https://quantum-journal.org/papers/q-2020-09-12-322
  Circuit cliff_circ;
  for (const Qubit &qb : qubits) {
    cliff_circ.add_qubit(qb);
  }
  ChoiMixTableau tab = tab_from_gadgets(gadgets);
  const MatrixXb &xmat = tab.tab_.xmat_;
  const MatrixXb &zmat = tab.tab_.zmat_;

  unsigned k = diagonalise_X_block(tab, cliff_circ, xmat, zmat);

  // Pairwise elimination method (Algorithm 2)
  for (unsigned i=1; i<k; ++i) {
    for (unsigned j=0; j<i; ++j) {
      if (zmat(i,j)) {
        tab.tab_.apply_gate(OpType::CZ, {i,j});
        Qubit i_qb = tab.col_index_.right.at(i).first;
        Qubit j_qb = tab.col_index_.right.at(j).first;
        cliff_circ.add_op<Qubit>(OpType::CZ, {i_qb, j_qb});
      }
    }
  }
  for (unsigned i=0; i<k; ++i) {
    Qubit i_qb = tab.col_index_.right.at(i).first;
    if (zmat(i,i)) {
      tab.tab_.apply_S(i);
      cliff_circ.add_op<Qubit>(OpType::S, {i_qb});
    }
    tab.tab_.apply_gate(OpType::H, {i});
    cliff_circ.add_op<Qubit>(OpType::H, {i_qb});
  }

  return cliff_circ;
}

Circuit mutual_diagonalise_vdBT_CX(
    std::list<SpSymPauliTensor> &gadgets, std::set<Qubit> qubits,
    CXConfigType) {
  // See Section 4 of van den Berg & Temme, "Circuit optimization of Hamiltonian simulation by simultaneous diagonalization of Pauli clusters", https://quantum-journal.org/papers/q-2020-09-12-322
  Circuit cliff_circ;
  for (const Qubit &qb : qubits) {
    cliff_circ.add_qubit(qb);
  }
  ChoiMixTableau tab = tab_from_gadgets(gadgets);
  const MatrixXb &xmat = tab.tab_.xmat_;
  const MatrixXb &zmat = tab.tab_.zmat_;

  unsigned k = diagonalise_X_block(tab, cliff_circ, xmat, zmat);

  // Elimination using CX operations (Algorithm 3)
  // Author's suggested using Patel-Markov-Hayes method for optimal CX synthesis, but we have omitted this for simplicity
  for (unsigned i=0; i<k; ++i) {
    unsigned n_ones = 0;
    for (unsigned j=0; j<=i; ++j) {
      if (zmat(i,j)) ++n_ones;
    }
    if (n_ones % 2 == 0) {
      tab.tab_.apply_S(i);
      Qubit i_qb = tab.col_index_.right.at(i).first;
      cliff_circ.add_op<Qubit>(OpType::S, {i_qb});
    }
    for (unsigned j=0; j<i; ++j) {
      if (zmat(i,j)) {
        tab.tab_.apply_CX(i,j);
        Qubit i_qb = tab.col_index_.right.at(i).first;
        Qubit j_qb = tab.col_index_.right.at(j).first;
        cliff_circ.add_op<Qubit>(OpType::CX, {i_qb, j_qb});
        tab.tab_.row_mult(j,i);
      }
    }
  }
  for (unsigned i=0; i<k; ++i) {
    tab.tab_.apply_S(i);
    tab.tab_.apply_gate(OpType::H, {i});
    Qubit i_qb = tab.col_index_.right.at(i).first;
    cliff_circ.add_op<Qubit>(OpType::S, {i_qb});
    cliff_circ.add_op<Qubit>(OpType::H, {i_qb});
  }

  return cliff_circ;
}

Circuit mutual_diagonalise_vdBT_greedy(
    std::list<SpSymPauliTensor> &gadgets, std::set<Qubit> qubits,
    CXConfigType, bool suborder_by_singles) {
  // See Section 4 of van den Berg & Temme, "Circuit optimization of Hamiltonian simulation by simultaneous diagonalization of Pauli clusters", https://quantum-journal.org/papers/q-2020-09-12-322
  Circuit cliff_circ;
  for (const Qubit &qb : qubits) {
    cliff_circ.add_qubit(qb);
  }
  ChoiMixTableau tab = tab_from_gadgets(gadgets);
  const MatrixXb &xmat = tab.tab_.xmat_;
  const MatrixXb &zmat = tab.tab_.zmat_;

  // std::cout << "diagonalise X" << std::endl;
  unsigned k = diagonalise_X_block(tab, cliff_circ, xmat, zmat);

  std::set<unsigned> active;
  for (unsigned i=0; i<k; ++i) active.insert(i);
  // unsigned iter = 0;
  while (!active.empty()) {
    // if (iter == 4) break;
    // ++iter;
    // std::cout << "Next search" << std::endl;
    unsigned min_cost = k;
    unsigned best_i = 0;
    auto best_it = active.end();
    unsigned best_j = 0;
    unsigned min_singles = 2*k;
    for (auto i_it = active.begin(); i_it != active.end(); ++i_it) {
      for (auto j_it = i_it; j_it != active.end(); ++j_it) {
        unsigned cost = 0;
        unsigned singles = 1; // Always involves a final H gate to remove X component
        if (*i_it == *j_it) {
          for (unsigned l = 0; l < k; ++l) {
            if (l != *i_it && zmat(l, *i_it)) ++cost;
          }
          if (suborder_by_singles) {
            // For the greedy-2 method, predict the number of single qubit gates required to solve column i by pairwise elimination
            if (zmat(*i_it, *i_it)) ++singles; // S gate
          }
        }
        else {
          ++cost;
          for (unsigned l = 0; l < k; ++l) {
            if (l != *i_it && l != *j_it && (zmat(l, *i_it) != zmat(l, *j_it))) ++cost;
          }
          if (suborder_by_singles) {
            // For the greedy-2 method, predict the number of single qubit gates required to solve column i using column j
            if (zmat(*i_it, *j_it) != zmat(*j_it, *j_it)) {
              // S required on j before sweeping
              ++singles;
            }
            if (zmat(*i_it, *i_it) != zmat(*j_it, *i_it)) {
              // S required on i after sweeping
              ++singles;
            }
          }
        }
        // std::cout << cost << std::endl;
        if (cost < min_cost || (suborder_by_singles && cost == min_cost && singles < min_singles)) {
          min_cost = cost;
          best_i = *i_it;
          best_it = i_it;
          best_j = *j_it;
          min_singles = singles;
        }
      }
    }
    // std::cout << "Best " << best_i << ", " << best_j << std::endl;

    Qubit i_qb = tab.col_index_.right.at(best_i).first;
    if (best_i != best_j) {
      // std::cout << "Sweep columns" << std::endl;
      Qubit j_qb = tab.col_index_.right.at(best_j).first;
      // Fix diagonal entry on j
      if (zmat(best_i, best_j) != zmat(best_j, best_j)) {
        tab.tab_.apply_S(best_j);
        cliff_circ.add_op<Qubit>(OpType::S, {j_qb});
      }
      // Sweep columns
      tab.tab_.apply_CX(best_i, best_j);
      cliff_circ.add_op<Qubit>(OpType::CX, {i_qb, j_qb});
      tab.tab_.row_mult(best_j, best_i);
      // Now minimal cost will be to solve column i outright, so we continue straight into that case
    }
    // std::cout << "Pairwise elimination" << std::endl;
    // Use simple case from pairwise elimination method
    for (unsigned l = 0; l < k; ++l) {
      if (l != best_i && zmat(best_i, l)) {
        tab.tab_.apply_gate(OpType::CZ, {best_i,l});
        Qubit l_qb = tab.col_index_.right.at(l).first;
        cliff_circ.add_op<Qubit>(OpType::CZ, {i_qb, l_qb});
      }
    }
    // std::cout << "Remove from X" << std::endl;
    if (zmat(best_i, best_i)) {
      tab.tab_.apply_S(best_i);
      cliff_circ.add_op<Qubit>(OpType::S, {i_qb});
    }
    tab.tab_.apply_gate(OpType::H, {best_i});
    cliff_circ.add_op<Qubit>(OpType::H, {i_qb});
    // std::cout << "Erase from active" << std::endl;
    active.erase(best_i);
  }

  return cliff_circ;
}

}  // namespace tket
