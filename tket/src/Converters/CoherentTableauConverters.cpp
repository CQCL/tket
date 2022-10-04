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

#include <stdexcept>

#include "Converters.hpp"

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
  CoherentTableau tab(t);
  qubit_vector_t input_qubits, output_qubits;
  for (unsigned i = 0; i < tab.get_n_boundaries(); ++i) {
    CoherentTableau::col_key_t key = tab.col_index_.right.at(i);
    if (key.second == CoherentTableau::TableauSegment::Input)
      input_qubits.push_back(key.first);
    else
      output_qubits.push_back(key.first);
  }
  Circuit in_circ(input_qubits, {});
  Circuit out_circ_tp(output_qubits, {});
  unit_map_t join_permutation;
  qubit_vector_t init_after_tp;
  std::set<CoherentTableau::col_key_t> solved;
  // Put tableau in Gaussian form to make synthesis identical irrespective of
  // current rows and to isolate any stabilisers over just inputs/just outputs
  // (i.e. corresponding to post-selection or pure initialisation) Column to
  // solve next
  unsigned col = 0;
  while (tab.get_n_rows() > 0) {
    // Identify qubit for column
    CoherentTableau::col_key_t key = tab.col_index_.right.at(col);
    if (solved.find(key) != solved.end()) {
      // Already handled this via a connection with a previous boundary
      ++col;
      continue;
    }
    // Isolate a single row with an X (if one exists)
    std::optional<unsigned> x_row = std::nullopt;
    for (unsigned r = 0; r < tab.get_n_rows(); ++r) {
      if (tab.tab_.xmat_(r, col)) {
        x_row = r;
        break;
      }
    }
    std::optional<Qubit> in_qb;
    std::optional<Qubit> out_qb;
    if (x_row) {
      // Make this row just X on that input/output.
      // A possible optimisation could involve row multiplications to reduce the
      // Hamming weight of the row before applying gates, but minimising this
      // would solve minimum weight/distance of a binary linear code whose
      // decision problem is NP-complete (shown by Vardy, "The intractability of
      // computing the minimum distance of a code", 1997). Just settle on using
      // the first row for now
      CoherentTableau::row_tensor_t row_paulis = tab.get_row(*x_row);
      // Start by fixing which qubits you are isolating as X
      if (key.second == CoherentTableau::TableauSegment::Input) {
        in_qb = key.first;
        if (row_paulis.second.string.map.size() != 0) {
          out_qb = row_paulis.second.string.map.begin()->first;
        }
      } else {
        out_qb = key.first;
        if (row_paulis.first.string.map.size() != 0) {
          in_qb = row_paulis.first.string.map.begin()->first;
        }
      }
      for (const std::pair<const Qubit, Pauli>& qbp :
           row_paulis.first.string.map) {
        if (qbp.first == *in_qb) {
          // If it is a Y, extract an Sdg gate so the Pauli is exactly X
          if (qbp.second == Pauli::Y) {
            in_circ.add_op<Qubit>(OpType::Sdg, {qbp.first});
            tab.apply_S(qbp.first, CoherentTableau::TableauSegment::Input);
          }
        } else {
          // Extract an entangling gate to eliminate the qubit
          switch (qbp.second) {
            case Pauli::X: {
              in_circ.add_op<Qubit>(OpType::CX, {*in_qb, qbp.first});
              tab.apply_CX(
                  *in_qb, qbp.first, CoherentTableau::TableauSegment::Input);
              break;
            }
            case Pauli::Y: {
              in_circ.add_op<Qubit>(OpType::CY, {*in_qb, qbp.first});
              tab.apply_gate(
                  OpType::CY, {*in_qb, qbp.first},
                  CoherentTableau::TableauSegment::Input);
              break;
            }
            case Pauli::Z: {
              in_circ.add_op<Qubit>(OpType::CZ, {*in_qb, qbp.first});
              tab.apply_gate(
                  OpType::CZ, {*in_qb, qbp.first},
                  CoherentTableau::TableauSegment::Input);
              break;
            }
            default: {
              break;
            }
          }
        }
      }
      // Same over the outputs
      for (const std::pair<const Qubit, Pauli>& qbp :
           row_paulis.second.string.map) {
        if (qbp.first == *out_qb) {
          // If it is a Y, extract an Sdg gate so the Pauli is exactly X
          if (qbp.second == Pauli::Y) {
            out_circ_tp.add_op<Qubit>(OpType::Sdg, {qbp.first});
            tab.apply_S(qbp.first, CoherentTableau::TableauSegment::Output);
          }
        } else {
          // Extract an entangling gate to eliminate the qubit
          switch (qbp.second) {
            case Pauli::X: {
              out_circ_tp.add_op<Qubit>(OpType::CX, {*out_qb, qbp.first});
              tab.apply_CX(
                  *in_qb, qbp.first, CoherentTableau::TableauSegment::Output);
              break;
            }
            case Pauli::Y: {
              // Transpose of CY is same as CY with a Z on control
              out_circ_tp.add_op<Qubit>(OpType::CY, {*out_qb, qbp.first});
              out_circ_tp.add_op<Qubit>(OpType::Z, {*out_qb});
              tab.apply_gate(
                  OpType::CY, {*out_qb, qbp.first},
                  CoherentTableau::TableauSegment::Output);
              tab.apply_gate(
                  OpType::Z, {*out_qb},
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
      // Multiply remaining rows with Xs to make it unique
      for (unsigned r = *x_row + 1; r < tab.get_n_rows(); ++r) {
        if (tab.tab_.xmat_(r, col)) {
          tab.tab_.row_mult(*x_row, r);
        }
      }
    }
    // And now do the same for Z
    std::optional<unsigned> z_row = std::nullopt;
    for (unsigned r = 0; r < tab.get_n_rows(); ++r) {
      if (tab.tab_.zmat_(r, col)) {
        z_row = r;
        break;
      }
    }
    if (z_row) {
      // If both an X and Z row exist, then both in_qb and out_qb should have
      // values and the rows should have anticommuting Paulis on each to
      // preserve commutativity of rows. If an X row existed over just the
      // inputs or outputs, no Z row could exist by commutativity of rows. So
      // only need to identify in_qb and out_qb if no X row existed
      CoherentTableau::row_tensor_t row_paulis = tab.get_row(*z_row);
      if (!x_row) {
        if (key.second == CoherentTableau::TableauSegment::Input) {
          in_qb = key.first;
          if (row_paulis.second.string.map.size() != 0) {
            out_qb = row_paulis.second.string.map.begin()->first;
          }
        } else {
          out_qb = key.first;
          if (row_paulis.first.string.map.size() != 0) {
            in_qb = row_paulis.first.string.map.begin()->first;
          }
        }
      }
      for (const std::pair<const Qubit, Pauli>& qbp :
           row_paulis.first.string.map) {
        if (qbp.first == *in_qb) {
          // If it is a Y, extract a Vdg gate so the Pauli is exactly Z
          if (qbp.second == Pauli::Y) {
            in_circ.add_op<Qubit>(OpType::Vdg, {qbp.first});
            tab.apply_V(qbp.first, CoherentTableau::TableauSegment::Input);
          }
        } else {
          // Extract an entangling gate to eliminate the qubit
          switch (qbp.second) {
            case Pauli::X: {
              in_circ.add_op<Qubit>(OpType::H, {qbp.first});
              in_circ.add_op<Qubit>(OpType::CX, {qbp.first, *in_qb});
              tab.apply_gate(
                  OpType::H, {qbp.first},
                  CoherentTableau::TableauSegment::Input);
              tab.apply_CX(
                  qbp.first, *in_qb, CoherentTableau::TableauSegment::Input);
              break;
            }
            case Pauli::Y: {
              in_circ.add_op<Qubit>(OpType::Vdg, {qbp.first});
              in_circ.add_op<Qubit>(OpType::CX, {qbp.first, *in_qb});
              tab.apply_V(qbp.first, CoherentTableau::TableauSegment::Input);
              tab.apply_CX(
                  qbp.first, *in_qb, CoherentTableau::TableauSegment::Input);
              break;
            }
            case Pauli::Z: {
              in_circ.add_op<Qubit>(OpType::CX, {qbp.first, *in_qb});
              tab.apply_CX(
                  qbp.first, *in_qb, CoherentTableau::TableauSegment::Input);
              break;
            }
            default: {
              break;
            }
          }
        }
      }
      // Same over the outputs
      for (const std::pair<const Qubit, Pauli>& qbp :
           row_paulis.second.string.map) {
        if (qbp.first == *out_qb) {
          // If it is a Y, extract an Vdg gate so the Pauli is exactly Z
          if (qbp.second == Pauli::Y) {
            out_circ_tp.add_op<Qubit>(OpType::Vdg, {qbp.first});
            tab.apply_V(qbp.first, CoherentTableau::TableauSegment::Output);
          }
        } else {
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
      // Multiply remaining rows with Zs to make it unique
      for (unsigned r = *z_row + 1; r < tab.get_n_rows(); ++r) {
        if (tab.tab_.zmat_(r, col)) {
          tab.tab_.row_mult(*z_row, r);
        }
      }
    }
    // Handle phases and resolve qubit connections
    if (x_row) {
      if (z_row) {
        // Hook up with an identity wire
        if (tab.tab_.phase_(*z_row)) in_circ.add_op<Qubit>(OpType::X, {*in_qb});
        if (tab.tab_.phase_(*x_row)) in_circ.add_op<Qubit>(OpType::Z, {*in_qb});
        join_permutation.insert({*out_qb, *in_qb});
        if (*x_row > *z_row) {
          // Remove highest index first to prevent changes to the indices
          // between removals
          tab.remove_row(*x_row);
          tab.remove_row(*z_row);
        } else {
          tab.remove_row(*z_row);
          tab.remove_row(*x_row);
        }
      } else {
        // Just an X row, but could be decoherence, post-selection, or
        // initialisation
        if (in_qb) {
          if (tab.tab_.phase_(*x_row))
            in_circ.add_op<Qubit>(OpType::Z, {*in_qb});
          if (out_qb) {
            // Decoherence
            in_circ.add_op<Qubit>(OpType::H, {*in_qb});
            in_circ.add_op<Qubit>(OpType::Collapse, {*in_qb});
            in_circ.add_op<Qubit>(OpType::H, {*in_qb});
            join_permutation.insert({*out_qb, *in_qb});
          } else {
            // Post-selection of input
            throw std::logic_error(
                "CoherentTableau synthesis requires post-selection which is "
                "not yet implemented");
          }
        } else {
          // At least one of in_qb and out_qb must be set, so out_qb is set here
          // Initialisation of output
          if (tab.tab_.phase_(*x_row))
            out_circ_tp.add_op<Qubit>(OpType::Z, {*out_qb});
          out_circ_tp.add_op<Qubit>(OpType::H, {*out_qb});
          init_after_tp.push_back(*out_qb);
        }
        tab.remove_row(*x_row);
      }
    } else {
      if (z_row) {
        // Just a Z row, but could be decoherence, post-selection, or
        // initialisation
        if (in_qb) {
          if (tab.tab_.phase_(*z_row))
            in_circ.add_op<Qubit>(OpType::X, {*in_qb});
          if (out_qb) {
            // Decoherence
            in_circ.add_op<Qubit>(OpType::Collapse, {*in_qb});
            join_permutation.insert({*out_qb, *in_qb});
          } else {
            // Post-selection of input
            throw std::logic_error(
                "CoherentTableau synthesis requires post-selection which is "
                "not yet implemented");
          }
        } else {
          // At least one of in_qb and out_qb must be set, so out_qb is set here
          // Initialisation of output
          if (tab.tab_.phase_(*z_row))
            out_circ_tp.add_op<Qubit>(OpType::X, {*out_qb});
          init_after_tp.push_back(*out_qb);
        }
        tab.remove_row(*z_row);
      } else {
        // No rows involving this qubit, so it is discarded or initialised in
        // the maximally mixed state
        if (in_qb) {
          in_circ.qubit_discard(*in_qb);
        } else {
          out_circ_tp.add_op<Qubit>(OpType::Collapse, {*out_qb});
          out_circ_tp.add_op<Qubit>(OpType::H, {*out_qb});
          init_after_tp.push_back(*out_qb);
        }
      }
    }
    ++col;
  }
  // Stitch subcircuits together
  Circuit out_circ = out_circ_tp.transpose();
  for (const Qubit& qb : init_after_tp) {
    out_circ.qubit_create(qb);
  }
  in_circ.append_with_map(out_circ, join_permutation);
  return {in_circ, join_permutation};
}

}  // namespace tket
