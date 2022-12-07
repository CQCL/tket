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

#include <boost/foreach.hpp>
#include <stdexcept>

#include "Converters.hpp"
#include "Diagonalisation/Diagonalisation.hpp"

namespace tket {

CoherentTableau circuit_to_coherent_tableau(const Circuit& circ) {
  CoherentTableau tab(circ.all_qubits());
  for (const Qubit& q : circ.created_qubits()) {
    tab.post_select(q, CoherentTableau::TableauSegment::Input);
  }
  for (const Command& com : circ) {
    auto args = com.get_args();
    qubit_vector_t qbs = {args.begin(), args.end()};
    tab.apply_gate(com.get_op_ptr()->get_type(), qbs);
  }
  for (const Qubit& q : circ.discarded_qubits()) {
    tab.discard_qubit(q);
  }
  return tab;
}

std::pair<Circuit, unit_map_t> coherent_tableau_to_circuit(
    const CoherentTableau& t) {
  /**
   * THE PLAN:
   * We first identify and solve all post-selections by Gaussian elimination and
   * diagonalisation of the input-only rows. This provides us with better
   * guarantees about the form of the stabilizers generated by the remaining
   * tableau - specifically that every possible stabilizer will involve at least
   * some output portion, so out_qb will always be set in later sections of the
   * synthesis. Diagonalisation of the input-only rows reduces them to just Z
   * strings but not necessarily over a minimal number of qubits. This can be
   * achieved by taking the Z matrix of these rows and performing row-wise
   * Gaussian elimination. The leading 1s indicate which qubits we are isolating
   * them onto and the rest gives the CXs required to reduce them to just the
   * leading qubits. After all post-selections have been identified, we enter
   * the main loop of handling all other inputs and reducing them one at a time
   * to either an identity wire or Z-decoherence (OpType::Collapse) connected to
   * an output, or to a discard. Let in_qb be this qubit to be solved. Pick a
   * row containing X_in_qb (if one exists). Since it must contain some
   * non-identity component in the output segment, we can pick one of these to
   * be out_qb and apply unitary gates at the input and output to reduce the row
   * to X_in_qb X_out_qb. We do the same with Z (if a row exists) to necessarily
   * give Z_in_qb Z_out_qb by commutativity properties, or we pick out_qb here
   * if no X row was found. If we had both an X row and a Z row, we have reduced
   * it to an identity wire just fine. If we have only one of them, e.g. X_in_qb
   * X_out_qb, we wish to make it the only row with X on either in_qb or out_qb
   * to leave a decoherence channel. We note that any other row containing
   * X_out_qb must have some other P_out2, since their combination cannot leave
   * an empty output segment. So we can apply a unitary gate to eliminate the
   * X_out_qb on the other row. By doing output-first gaussian elimination on
   * the output sub-tableau ignoring the target row and X_out_qb column, for
   * each such row with X_out_qb there is some other output column P_out2 for
   * which it is the unique entry, so we can be sure that applying the unitary
   * gate does not add X_out_qb onto other rows. Having this strategy available
   * to make it the unique row with X_out_qb still allows us to make it the
   * unique row with X_in_qb by row combinations. Once all rows with inputs have
   * been eliminated, any unsolved inputs must have been discarded and the
   * remaining tableau is an inverse diagonalisation tableau (only outputs).
   */

  // Operate on a copy of the tableau to track the remainder as we extract gates
  CoherentTableau tab(t);

  // Canonicalise tableau and perform output-first gaussian elimination to
  // isolate post-selected subspace in lower rows
  tab.canonical_column_order(CoherentTableau::TableauSegment::Output);
  tab.gaussian_form();

  // Set up circuits for extracting gates into (in_circ is set up after
  // diagonalising the post-selected subspace)
  qubit_vector_t input_qubits, output_qubits;
  for (unsigned i = 0; i < tab.get_n_boundaries(); ++i) {
    CoherentTableau::col_key_t key = tab.col_index_.right.at(i);
    if (key.second == CoherentTableau::TableauSegment::Input)
      input_qubits.push_back(key.first);
    else
      output_qubits.push_back(key.first);
  }
  Circuit out_circ_tp(output_qubits, {});
  unit_map_t join_permutation;
  std::set<Qubit> post_selected;
  std::map<Qubit, std::optional<unsigned>> in_x_row;
  std::map<Qubit, std::optional<unsigned>> in_z_row;
  boost::bimap<Qubit, Qubit> matched_qubits;

  // Call diagonalisation methods to diagonalise post-selected subspace
  std::list<std::pair<QubitPauliTensor, Expr>> to_diag;
  for (unsigned r = tab.get_n_rows(); r > 0;) {
    --r;
    CoherentTableau::row_tensor_t rten = tab.get_row(r);
    if (rten.second.string.map.size() != 0) {
      // Reached the rows with non-empty output segment
      break;
    }
    // Else, we add the row to the subspace
    to_diag.push_back({rten.first, 1.});
  }
  unsigned post_selected_size = to_diag.size();
  std::set<Qubit> diag_ins{input_qubits.begin(), input_qubits.end()};
  Circuit in_circ = mutual_diagonalise(to_diag, diag_ins, CXConfigType::Tree);
  // Extract the dagger of each gate in order from tab
  for (const Command& com : in_circ) {
    auto args = com.get_args();
    qubit_vector_t qbs = {args.begin(), args.end()};
    tab.apply_gate(
        com.get_op_ptr()->dagger()->get_type(), qbs,
        CoherentTableau::TableauSegment::Input);
  }

  // Diagonalised rows should still be at the bottom of the tableau - reduce
  // them to a minimal set of qubits for post-selection by first reducing to
  // upper echelon form
  std::vector<std::pair<unsigned, unsigned>> row_ops =
      gaussian_elimination_row_ops(
          tab.tab_.zmat_.bottomRows(post_selected_size));
  for (const std::pair<unsigned, unsigned>& op : row_ops) {
    tab.tab_.row_mult(op.first, op.second);
  }
  // Obtain CX instructions as column operations
  std::vector<std::pair<unsigned, unsigned>> col_ops =
      gaussian_elimination_col_ops(
          tab.tab_.zmat_.bottomRows(post_selected_size));
  // These gates will also swap qubits to isolate the post-selections on the
  // first few qubits - this is fine as can be cleaned up later with peephole
  // optimisations
  for (const std::pair<unsigned, unsigned>& op : col_ops) {
    tab.tab_.apply_CX(op.first, op.second);
    CoherentTableau::col_key_t ctrl = tab.col_index_.right.at(op.first);
    CoherentTableau::col_key_t trgt = tab.col_index_.right.at(op.second);
    in_circ.add_op<Qubit>(OpType::CX, {ctrl.first, trgt.first});
  }

  // Post-select rows
  // TODO Post-selection op is not yet available in tket - replace this once
  // implemented
  if (post_selected_size != 0)
    throw std::logic_error(
        "Not yet implemented: post-selection required during CoherentTableau "
        "synthesis");
  // for (unsigned r = 0; r < post_selected_size; ++r) {
  //   CoherentTableau::row_tensor_t row = tab.get_row(tab.get_n_rows() - 1);
  //   if (row.second.string.map.size() != 0 || row.first.string.map.size() != 1
  //   || row.first.string.map.begin()->second != Pauli::Z) throw
  //   std::logic_error("Unexpected error during post-selection identification
  //   in CoherentTableau synthesis"); Qubit post_selected_qb =
  //   row.first.string.map.begin()->first; if (row.second.coeff == -1.)
  //   in_circ.add_op<Qubit>(OpType::X, {post_selected_qb});
  //   tab.remove_row(tab.get_n_rows() - 1);
  //   post_selected.insert(post_selected_qb);
  //   throw std::logic_error("Not yet implemented: post-selection required
  //   during CoherentTableau synthesis");
  // }

  // Input-first gaussian elimination to solve input-sides of remaining rows
  tab.canonical_column_order(CoherentTableau::TableauSegment::Input);
  tab.gaussian_form();

  // Iterate through remaining inputs and reduce output portion to a single
  // qubit
  for (const Qubit& in_qb : input_qubits) {
    // Skip post-selected qubits
    if (post_selected.find(in_qb) != post_selected.end()) continue;

    unsigned col = tab.col_index_.left.at(CoherentTableau::col_key_t{
        in_qb, CoherentTableau::TableauSegment::Input});
    std::optional<Qubit> out_qb = std::nullopt;

    // Find the row with X_in_qb (if one exists)
    std::optional<unsigned> x_row = std::nullopt;
    for (unsigned r = 0; r < tab.get_n_rows(); ++r) {
      if (tab.tab_.xmat_(r, col)) {
        x_row = r;
        break;
      }
    }

    if (x_row) {
      // A possible optimisation could involve row multiplications to reduce the
      // Hamming weight of the row before applying gates, but minimising this
      // would solve minimum weight/distance of a binary linear code whose
      // decision problem is NP-complete (shown by Vardy, "The intractability of
      // computing the minimum distance of a code", 1997). Just settle on using
      // the first row for now, reducing the input and output to a single qubit
      CoherentTableau::row_tensor_t row_paulis = tab.get_row(*x_row);
      for (const std::pair<const Qubit, Pauli>& p :
           row_paulis.second.string.map) {
        if (matched_qubits.right.find(p.first) == matched_qubits.right.end()) {
          out_qb = p.first;
          matched_qubits.insert({in_qb, *out_qb});
          break;
        }
      }

      // Reduce input string to just X_in_qb
      if (row_paulis.first.string.get(in_qb) == Pauli::Y) {
        // If it is a Y, extract an Sdg gate so the Pauli is exactly X
        in_circ.add_op<Qubit>(OpType::Sdg, {in_qb});
        tab.apply_S(in_qb, CoherentTableau::TableauSegment::Input);
      }
      for (const std::pair<const Qubit, Pauli>& qbp :
           row_paulis.first.string.map) {
        if (qbp.first == in_qb) continue;
        // Extract an entangling gate to eliminate the qubit
        switch (qbp.second) {
          case Pauli::X: {
            in_circ.add_op<Qubit>(OpType::CX, {in_qb, qbp.first});
            tab.apply_CX(
                in_qb, qbp.first, CoherentTableau::TableauSegment::Input);
            break;
          }
          case Pauli::Y: {
            in_circ.add_op<Qubit>(OpType::CY, {in_qb, qbp.first});
            tab.apply_gate(
                OpType::CY, {in_qb, qbp.first},
                CoherentTableau::TableauSegment::Input);
            break;
          }
          case Pauli::Z: {
            in_circ.add_op<Qubit>(OpType::CZ, {in_qb, qbp.first});
            tab.apply_gate(
                OpType::CZ, {in_qb, qbp.first},
                CoherentTableau::TableauSegment::Input);
            break;
          }
          default: {
            break;
          }
        }
      }

      // And then the same for X_out_qb
      if (row_paulis.second.string.get(*out_qb) == Pauli::Y) {
        // If it is a Y, extract an Sdg gate so the Pauli is exactly X
        out_circ_tp.add_op<Qubit>(OpType::Sdg, {*out_qb});
        tab.apply_S(*out_qb, CoherentTableau::TableauSegment::Output);
      } else if (row_paulis.second.string.get(*out_qb) == Pauli::Z) {
        // If it is a Z, extract an Vdg and Sdg gate so the Pauli is exactly X
        out_circ_tp.add_op<Qubit>(OpType::Vdg, {*out_qb});
        out_circ_tp.add_op<Qubit>(OpType::Sdg, {*out_qb});
        tab.apply_V(*out_qb, CoherentTableau::TableauSegment::Output);
        tab.apply_S(*out_qb, CoherentTableau::TableauSegment::Output);
      }
      for (const std::pair<const Qubit, Pauli>& qbp :
           row_paulis.second.string.map) {
        if (qbp.first == *out_qb) continue;
        // Extract an entangling gate to eliminate the qubit
        switch (qbp.second) {
          case Pauli::X: {
            out_circ_tp.add_op<Qubit>(OpType::CX, {*out_qb, qbp.first});
            tab.apply_CX(
                *out_qb, qbp.first, CoherentTableau::TableauSegment::Output);
            break;
          }
          case Pauli::Y: {
            // CY does not have a transpose OpType defined so decompose
            out_circ_tp.add_op<Qubit>(OpType::S, {qbp.first});
            out_circ_tp.add_op<Qubit>(OpType::CX, {*out_qb, qbp.first});
            out_circ_tp.add_op<Qubit>(OpType::Sdg, {qbp.first});
            tab.apply_gate(
                OpType::CY, {*out_qb, qbp.first},
                CoherentTableau::TableauSegment::Output);
            break;
          }
          case Pauli::Z: {
            out_circ_tp.add_op<Qubit>(OpType::CZ, {*out_qb, qbp.first});
            tab.apply_gate(
                OpType::CZ, {*out_qb, qbp.first},
                CoherentTableau::TableauSegment::Output);
            break;
          }
          default: {
            break;
          }
        }
      }
    }

    // Find the row with Z_in_qb (if one exists)
    std::optional<unsigned> z_row = std::nullopt;
    for (unsigned r = 0; r < tab.get_n_rows(); ++r) {
      if (tab.tab_.zmat_(r, col)) {
        z_row = r;
        break;
      }
    }
    if (z_row) {
      // If both an X and Z row exist, then out_qb should have a
      // value and the rows should have anticommuting Paulis on out_qb to
      // preserve commutativity of rows
      CoherentTableau::row_tensor_t row_paulis = tab.get_row(*z_row);
      if (!x_row) {
        for (const std::pair<const Qubit, Pauli>& p :
             row_paulis.second.string.map) {
          if (matched_qubits.right.find(p.first) ==
              matched_qubits.right.end()) {
            out_qb = p.first;
            matched_qubits.insert({in_qb, *out_qb});
            break;
          }
        }
      }

      // Reduce input string to just Z_in_qb.
      // No need to consider different paulis on in_qb: if we had X or Y, this
      // row would have been identified as x_row instead of Z row, and if
      // another row was already chosen as x_row then canonical gaussian form
      // would imply all other rows do not contain X_in_qb or Y_in_qb
      for (const std::pair<const Qubit, Pauli>& qbp :
           row_paulis.first.string.map) {
        if (qbp.first == in_qb) continue;
        // Extract an entangling gate to eliminate the qubit
        switch (qbp.second) {
          case Pauli::X: {
            in_circ.add_op<Qubit>(OpType::H, {qbp.first});
            in_circ.add_op<Qubit>(OpType::CX, {qbp.first, in_qb});
            tab.apply_gate(
                OpType::H, {qbp.first}, CoherentTableau::TableauSegment::Input);
            tab.apply_CX(
                qbp.first, in_qb, CoherentTableau::TableauSegment::Input);
            break;
          }
          case Pauli::Y: {
            in_circ.add_op<Qubit>(OpType::Vdg, {qbp.first});
            in_circ.add_op<Qubit>(OpType::CX, {qbp.first, in_qb});
            tab.apply_V(qbp.first, CoherentTableau::TableauSegment::Input);
            tab.apply_CX(
                qbp.first, in_qb, CoherentTableau::TableauSegment::Input);
            break;
          }
          case Pauli::Z: {
            in_circ.add_op<Qubit>(OpType::CX, {qbp.first, in_qb});
            tab.apply_CX(
                qbp.first, in_qb, CoherentTableau::TableauSegment::Input);
            break;
          }
          default: {
            break;
          }
        }
      }

      // And then reduce output string to just Z_out_qb
      if (row_paulis.second.string.get(*out_qb) == Pauli::Y) {
        // If it is a Y, extract a Vdg gate so the Pauli is exactly Z
        out_circ_tp.add_op<Qubit>(OpType::Vdg, {*out_qb});
        tab.apply_V(*out_qb, CoherentTableau::TableauSegment::Output);
      } else if (row_paulis.second.string.get(*out_qb) == Pauli::X) {
        // If it is an X, extract an Sdg and Vdg gate so the Pauli is exactly Z
        // We do not need to care about messing up the X row here since if we
        // solved an X row then this row can't also have X by commutativity
        out_circ_tp.add_op<Qubit>(OpType::Sdg, {*out_qb});
        out_circ_tp.add_op<Qubit>(OpType::Vdg, {*out_qb});
        tab.apply_S(*out_qb, CoherentTableau::TableauSegment::Output);
        tab.apply_V(*out_qb, CoherentTableau::TableauSegment::Output);
      }
      for (const std::pair<const Qubit, Pauli>& qbp :
           row_paulis.second.string.map) {
        if (qbp.first == *out_qb) continue;
        // Extract an entangling gate to eliminate the qubit
        switch (qbp.second) {
          case Pauli::X: {
            out_circ_tp.add_op<Qubit>(OpType::H, {qbp.first});
            out_circ_tp.add_op<Qubit>(OpType::CX, {qbp.first, *out_qb});
            tab.apply_gate(
                OpType::H, {qbp.first},
                CoherentTableau::TableauSegment::Output);
            tab.apply_CX(
                qbp.first, *out_qb, CoherentTableau::TableauSegment::Output);
            break;
          }
          case Pauli::Y: {
            out_circ_tp.add_op<Qubit>(OpType::Vdg, {qbp.first});
            out_circ_tp.add_op<Qubit>(OpType::CX, {qbp.first, *out_qb});
            tab.apply_V(qbp.first, CoherentTableau::TableauSegment::Output);
            tab.apply_CX(
                qbp.first, *out_qb, CoherentTableau::TableauSegment::Output);
            break;
          }
          case Pauli::Z: {
            out_circ_tp.add_op<Qubit>(OpType::CX, {qbp.first, *out_qb});
            tab.apply_CX(
                qbp.first, *out_qb, CoherentTableau::TableauSegment::Output);
            break;
          }
          default: {
            break;
          }
        }
      }
    }

    in_x_row.insert({in_qb, x_row});
    in_z_row.insert({in_qb, z_row});
  }

  // X and Z row are guaranteed to be the unique rows with X_in_qb and Z_in_qb.
  // If the Z row exists, then all rows must commute with Z_in_qb Z_out_qb, so
  // if any row has X_out_qb it must also have X_in_qb. Hence if both exist then
  // the X row is also the unique row with X_out_qb, and by similar argument the
  // Z row is the unique row with Z_out_qb. However, if only one row exists,
  // this uniqueness may not hold but can be forced by applying some unitary
  // gates to eliminate e.g. Z_out_qb from all other rows. We use a gaussian
  // elimination subroutine to identify a combination of gates that won't add
  // Z_out_qb onto other qubits. Since we have already removed them from all
  // other rows containing input components, we only need to consider the
  // output-only rows, which all exist at the bottom of the tableau.
  unsigned out_stabs = 0;
  while (out_stabs < tab.get_n_rows()) {
    CoherentTableau::row_tensor_t rten =
        tab.get_row(tab.get_n_rows() - 1 - out_stabs);
    if (rten.first.string.map.size() != 0) {
      // Reached the rows with non-empty input segment
      break;
    }
    ++out_stabs;
  }
  MatrixXb out_rows = MatrixXb::Zero(out_stabs, 2 * output_qubits.size());
  out_rows(Eigen::all, Eigen::seq(0, Eigen::last, 2)) = tab.tab_.xmat_.block(
      tab.get_n_rows() - out_stabs, input_qubits.size(), out_stabs,
      output_qubits.size());
  out_rows(Eigen::all, Eigen::seq(1, Eigen::last, 2)) = tab.tab_.zmat_.block(
      tab.get_n_rows() - out_stabs, input_qubits.size(), out_stabs,
      output_qubits.size());
  // Identify other qubits to apply gates to for removing the extra matched
  // output terms
  MatrixXb unmatched_outs = MatrixXb::Zero(
      out_stabs, 2 * (output_qubits.size() - matched_qubits.size()));
  std::vector<std::pair<Qubit, Pauli>> col_lookup;
  for (unsigned out_col = 0; out_col < output_qubits.size(); ++out_col) {
    unsigned tab_col = input_qubits.size() + out_col;
    Qubit out_qb = tab.col_index_.right.at(tab_col).first;
    if (matched_qubits.right.find(out_qb) != matched_qubits.right.end())
      continue;
    unmatched_outs.col(col_lookup.size()) = out_rows.col(2 * out_col);
    col_lookup.push_back({out_qb, Pauli::X});
    unmatched_outs.col(col_lookup.size()) = out_rows.col(2 * out_col + 1);
    col_lookup.push_back({out_qb, Pauli::Z});
  }
  row_ops = gaussian_elimination_row_ops(unmatched_outs);
  for (const std::pair<unsigned, unsigned>& op : row_ops) {
    for (unsigned c = 0; c < unmatched_outs.cols(); ++c) {
      unmatched_outs(op.second, c) =
          unmatched_outs(op.second, c) ^ unmatched_outs(op.first, c);
    }
    tab.tab_.row_mult(
        tab.get_n_rows() - out_stabs + op.first,
        tab.get_n_rows() - out_stabs + op.second);
  }
  // Go through the output-only rows and remove terms from matched qubits
  for (unsigned r = 0; r < out_stabs; ++r) {
    unsigned leading_col = 0;
    for (unsigned c = 0; c < unmatched_outs.cols(); ++c) {
      if (unmatched_outs(r, c)) {
        leading_col = c;
        break;
      }
    }
    std::pair<Qubit, Pauli> alternate_qb = col_lookup.at(leading_col);

    if (alternate_qb.second != Pauli::X) {
      // Make the alternate point of contact X so we only need one set of rule
      // for eliminating Paulis
      tab.apply_gate(
          OpType::H, {alternate_qb.first},
          CoherentTableau::TableauSegment::Output);
      out_circ_tp.add_op<Qubit>(OpType::H, {alternate_qb.first});
    }

    CoherentTableau::row_tensor_t row_paulis =
        tab.get_row(tab.get_n_rows() - out_stabs + r);

    for (const std::pair<const Qubit, Pauli>& qbp :
         row_paulis.second.string.map) {
      if (matched_qubits.right.find(qbp.first) == matched_qubits.right.end())
        continue;
      // Alternate point is guaranteed to be unmatched, so always needs an
      // entangling gate
      switch (qbp.second) {
        case Pauli::X: {
          out_circ_tp.add_op<Qubit>(
              OpType::CX, {alternate_qb.first, qbp.first});
          tab.apply_CX(
              alternate_qb.first, qbp.first,
              CoherentTableau::TableauSegment::Output);
          break;
        }
        case Pauli::Z: {
          out_circ_tp.add_op<Qubit>(
              OpType::CZ, {alternate_qb.first, qbp.first});
          tab.apply_gate(
              OpType::CZ, {alternate_qb.first, qbp.first},
              CoherentTableau::TableauSegment::Output);
          break;
        }
        default: {
          // Don't have to care about Y since any matched qubit has a row that
          // is reduced to either X or Z and all other rows must commute with
          // that
          break;
        }
      }
    }
  }

  // Now that X_in_qb X_out_qb (or Zs) is the unique row for each of X_in_qb and
  // X_out_qb, we can actually link up the qubit wires and remove the rows
  for (const Qubit& in_qb : input_qubits) {
    if (post_selected.find(in_qb) != post_selected.end()) continue;

    std::optional<unsigned> x_row = in_x_row.at(in_qb);
    std::optional<unsigned> z_row = in_z_row.at(in_qb);
    auto found = matched_qubits.left.find(in_qb);
    std::optional<Qubit> out_qb = (found == matched_qubits.left.end())
                                      ? std::nullopt
                                      : std::optional<Qubit>{found->second};
    // Handle phases and resolve qubit connections
    if (x_row) {
      if (z_row) {
        // Hook up with an identity wire
        if (tab.tab_.phase_(*z_row)) in_circ.add_op<Qubit>(OpType::X, {in_qb});
        if (tab.tab_.phase_(*x_row)) in_circ.add_op<Qubit>(OpType::Z, {in_qb});
        join_permutation.insert({*out_qb, in_qb});
      } else {
        // Just an X row, so must be connected to out_qb via a decoherence
        if (tab.tab_.phase_(*x_row)) in_circ.add_op<Qubit>(OpType::Z, {in_qb});
        in_circ.add_op<Qubit>(OpType::H, {in_qb});
        in_circ.add_op<Qubit>(OpType::Collapse, {in_qb});
        in_circ.add_op<Qubit>(OpType::H, {in_qb});
        join_permutation.insert({*out_qb, in_qb});
      }
    } else {
      if (z_row) {
        // Just a Z row, so must be connected to out_qb via a decoherence
        if (tab.tab_.phase_(*z_row)) in_circ.add_op<Qubit>(OpType::X, {in_qb});
        in_circ.add_op<Qubit>(OpType::Collapse, {in_qb});
        join_permutation.insert({*out_qb, in_qb});
      } else {
        // No rows involving this input, so it is discarded
        in_circ.qubit_discard(in_qb);
      }
    }
  }

  // Remove rows with inputs from the tableau
  tab.tab_ = SymplecticTableau(
      tab.tab_.xmat_.bottomRows(out_stabs),
      tab.tab_.zmat_.bottomRows(out_stabs), tab.tab_.phase_.tail(out_stabs));

  // Can't use a template with multiple parameters within a macro since the
  // comma will register as an argument delimiter for the macro
  using match_entry = boost::bimap<Qubit, Qubit>::left_const_reference;
  BOOST_FOREACH (match_entry entry, matched_qubits.left) {
    tab.discard_qubit(entry.first, CoherentTableau::TableauSegment::Input);
    tab.discard_qubit(entry.second, CoherentTableau::TableauSegment::Output);
  }
  tab.canonical_column_order(CoherentTableau::TableauSegment::Output);

  // Only remaining rows must be completely over the outputs. Call
  // diagonalisation methods to diagonalise coherent subspace
  to_diag.clear();
  for (unsigned r = 0; r < tab.get_n_rows(); ++r) {
    CoherentTableau::row_tensor_t rten = tab.get_row(r);
    to_diag.push_back({rten.second, 1.});
  }
  std::set<Qubit> diag_outs;
  for (const Qubit& out : output_qubits) {
    if (matched_qubits.right.find(out) == matched_qubits.right.end())
      diag_outs.insert(out);
  }
  Circuit out_diag_circ =
      mutual_diagonalise(to_diag, diag_outs, CXConfigType::Tree);
  // Extract the dagger of each gate in order from tab
  for (const Command& com : out_diag_circ) {
    auto args = com.get_args();
    qubit_vector_t qbs = {args.begin(), args.end()};
    tab.apply_gate(com.get_op_ptr()->get_type(), qbs);
    out_circ_tp.add_op<Qubit>(com.get_op_ptr()->transpose()->dagger(), qbs);
  }

  // All rows are diagonalised so we can just focus on the Z matrix. Reduce them
  // to a minimal set of qubits for initialisation by first reducing to upper
  // echelon form
  row_ops = gaussian_elimination_row_ops(tab.tab_.zmat_);
  for (const std::pair<unsigned, unsigned>& op : row_ops) {
    tab.tab_.row_mult(op.first, op.second);
  }
  // Obtain CX instructions as column operations
  col_ops = gaussian_elimination_col_ops(tab.tab_.zmat_);
  for (const std::pair<unsigned, unsigned>& op : col_ops) {
    tab.tab_.apply_CX(op.second, op.first);
    CoherentTableau::col_key_t ctrl = tab.col_index_.right.at(op.second);
    CoherentTableau::col_key_t trgt = tab.col_index_.right.at(op.first);
    out_circ_tp.add_op<Qubit>(OpType::CX, {ctrl.first, trgt.first});
  }

  // Fix phases of zero_initialised qubits
  std::set<Qubit> zero_initialised;
  Circuit out_circ(output_qubits, {});
  for (unsigned r = 0; r < tab.get_n_rows(); ++r) {
    CoherentTableau::row_tensor_t row = tab.get_row(r);
    if (row.first.string.map.size() != 0 || row.second.string.map.size() != 1 ||
        row.second.string.map.begin()->second != Pauli::Z)
      throw std::logic_error(
          "Unexpected error during zero initialisation in CoherentTableau "
          "synthesis");
    Qubit initialised_qb = row.second.string.map.begin()->first;
    out_circ.qubit_create(initialised_qb);
    if (row.second.coeff == -1.) {
      out_circ.add_op<Qubit>(OpType::X, {initialised_qb});
    }
    zero_initialised.insert(initialised_qb);
  }

  // Remaining outputs that aren't zero initialised or matched need to be
  // initialised in the maximally-mixed state.
  // Also match up unmatched outputs to either unmatched inputs or reusable
  // output names (ones that are already matched up to other input names),
  // preferring the qubit of the same name
  std::list<Qubit> reusable_names;
  for (const Qubit& out_qb : output_qubits) {
    if (matched_qubits.right.find(out_qb) != matched_qubits.right.end() &&
        matched_qubits.left.find(out_qb) == matched_qubits.left.end())
      reusable_names.push_back(out_qb);
  }
  for (const Qubit& out_qb : output_qubits) {
    if (matched_qubits.right.find(out_qb) == matched_qubits.right.end()) {
      if (zero_initialised.find(out_qb) == zero_initialised.end()) {
        out_circ.qubit_create(out_qb);
        out_circ.add_op<Qubit>(OpType::H, {out_qb});
        out_circ.add_op<Qubit>(OpType::Collapse, {out_qb});
      }
      if (matched_qubits.left.find(out_qb) == matched_qubits.left.end()) {
        matched_qubits.insert({out_qb, out_qb});
        join_permutation.insert({out_qb, out_qb});
      } else {
        // Since the matching is a bijection, matched_qubits.left -
        // matched_qubits.right (set difference) is the same size as
        // matched_qubits.right - matched_qubits.left, so there will be exactly
        // the right number of reusable names to pull from here
        Qubit name = reusable_names.front();
        reusable_names.pop_front();
        matched_qubits.insert({name, out_qb});
        join_permutation.insert({out_qb, name});
      }
    }
  }

  // Initialise qubits with stabilizer rows and stitch subcircuits together
  out_circ.append(out_circ_tp.transpose());
  in_circ.append_with_map(out_circ, join_permutation);
  return {in_circ, join_permutation};
}

}  // namespace tket
