// Copyright 2019-2023 Cambridge Quantum Computing
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

#include "tket/Converters/Converters.hpp"
#include "tket/Diagonalisation/Diagonalisation.hpp"

namespace tket {

ChoiMixTableau circuit_to_cm_tableau(const Circuit& circ) {
  ChoiMixTableau tab(circ.all_qubits());
  for (const Qubit& q : circ.created_qubits()) {
    tab.post_select(q, ChoiMixTableau::TableauSegment::Input);
  }
  for (const Command& com : circ) {
    auto args = com.get_args();
    qubit_vector_t qbs = {args.begin(), args.end()};
    tab.apply_gate(com.get_op_ptr()->get_type(), qbs);
  }
  tab.rename_qubits(
      circ.implicit_qubit_permutation(),
      ChoiMixTableau::TableauSegment::Output);
  for (const Qubit& q : circ.discarded_qubits()) {
    tab.discard_qubit(q);
  }
  tab.canonical_column_order();
  tab.gaussian_form();
  return tab;
}

struct ChoiMixBuilder {
  /**
   * We will consider applying gates to either side of the tableau to reduce
   * qubits down to one of a few simple states (identity, collapse, zero
   * initialise, mixed initialised, post-selected, or discarded), allowing us to
   * remove that qubit and continue until the tableau contains no qubits left.
   * This gradually builds up a set of operations both before and after the
   * working tableau. We should ask for the following combination of temporary
   * states to compose to form a channel equivalent to that described by the
   * input tableau:
   * - in_circ: a unitary circuit (without implicit permutations)
   * - post_selected: a set of post-selection actions to apply to the outputs of
   * in_circ
   * - discarded: a set of discard actions to apply to the outputs of in_circ
   * - collapsed: a set of collapse actions to apply to the outputs of in_circ
   * (i.e. decoherence in the Z basis)
   * - tab: the remaining tableau still to be solved. Acting as an identity on
   * any qubits not contained within tab
   * - in_out_permutation: a permutation of the qubits, read as a map from the
   * input qubit name to the output qubit it is sent to
   * - zero_initialised: a set of initialisations of fresh output qubits
   * (in_out_permutation may join these qubits onto input qubits that have been
   * post-selected or discarded to reuse qubits)
   * - mix_initialised: a set of initialisations of output qubits into
   * maximally-mixed states (in_out_permutation may similarly join these onto
   * input qubits no longer in use)
   * - out_circ_tp: the transpose of a unitary circuit (without implicit
   * permuations); we store this as the transpose so we can build it up in
   * reverse
   */
  Circuit in_circ;
  std::set<Qubit> post_selected;
  std::set<Qubit> discarded;
  std::set<Qubit> collapsed;
  ChoiMixTableau tab;
  boost::bimap<Qubit, Qubit> in_out_permutation;
  std::set<Qubit> zero_initialised;
  std::set<Qubit> mix_initialised;
  Circuit out_circ_tp;

  // The CXConfigType preferred when invoking diagonalisation techniques
  CXConfigType cx_config;

  // Additional qubit names (distinct from qubits already on the respective
  // segment of the tableau) than can be used for zero initialisations and
  // post-selections when synthesising a unitary extension
  qubit_vector_t unitary_init_names;
  qubit_vector_t unitary_post_names;

  // Initialises the builder with a tableau
  explicit ChoiMixBuilder(const ChoiMixTableau& tab, CXConfigType cx_config);
  // For synthesis of a unitary extension, initialises the builder with a
  // tableau and some additional qubit names which the resulting circuit may
  // optionally use, representing zero-initialised or post-selected qubits.
  // These are the qubits on which we can freely add Zs to the rows of the given
  // tableau to guarantee a unitary extension; none of these names should appear
  // in tab
  explicit ChoiMixBuilder(
      const ChoiMixTableau& tab, CXConfigType cx_config,
      const qubit_vector_t& init_names, const qubit_vector_t& post_names);

  // Debug method: applies all staged operations back onto tab to provide the
  // tableau that the synthesis result is currently aiming towards. In exact
  // synthesis, this should remain invariant during synthesis. For synthesis of
  // a unitary extension, this should at least span the rows of the original
  // tableau up to additional Zs on spare input and output qubits.
  ChoiMixTableau realised_tableau() const;

  /**
   * STAGES OF SYNTHESIS
   */

  // Match up pairs of generators that anti-commute in the input segment but
  // commute with all others; such pairs of rows reduce to an identity wire
  // between a pair of qubits; we solve this by pairwise Pauli reduction methods
  void solve_id_subspace();
  // After removing the identity subspace, all remaining rows mutually commute
  // within each tableau segment; diagonalise each segment individually
  void diagonalise_segments();
  // Solve the post-selected subspace which has already been diagonalised
  void solve_postselected_subspace();
  // Solve the zero-initialised subspace which has already been diagonalised
  void solve_initialised_subspace();
  // All remaining rows are in the collapsed subspace (each row is the unique
  // stabilizer passing through some Collapse gate); solve it
  void solve_collapsed_subspace();

  // Simplifies the tableau by removing qubits on which all rows have I; such
  // qubits are either discarded inputs or mixed-initialised outputs
  void remove_unused_qubits();
  // For synthesis of a unitary extension, match up qubits from
  // post-selected/zero-initialised with unitary_post_names/unitary_init_names
  // and add them to in_out_permutation
  void assign_init_post_names();
  // Fill out in_out_permutation to map all qubits; this typically takes the
  // form of a standard qubit reuse pattern (e.g. discard and reinitialise)
  void assign_remaining_names();
  // Once tab has been completely reduced to no rows and no qubits, compose the
  // staged operations to build the output circuit and return the renaming map
  // from output names of original tableau to the qubits of the returned circuit
  // they are mapped to
  std::pair<Circuit, qubit_map_t> output_circuit();
  std::pair<Circuit, qubit_map_t> unitary_output_circuit();
};

std::pair<Circuit, qubit_map_t> cm_tableau_to_exact_circuit(
    const ChoiMixTableau& tab, CXConfigType cx_config) {
  ChoiMixBuilder builder(tab, cx_config);
  builder.remove_unused_qubits();
  builder.solve_id_subspace();
  builder.diagonalise_segments();
  builder.solve_postselected_subspace();
  builder.solve_initialised_subspace();
  builder.solve_collapsed_subspace();
  builder.remove_unused_qubits();
  builder.assign_remaining_names();
  return builder.output_circuit();
}

std::pair<Circuit, qubit_map_t> cm_tableau_to_unitary_extension_circuit(
    const ChoiMixTableau& tab, const std::vector<Qubit>& init_names,
    const std::vector<Qubit>& post_names, CXConfigType cx_config) {
  ChoiMixBuilder builder(tab, cx_config, init_names, post_names);
  builder.remove_unused_qubits();
  builder.solve_id_subspace();
  builder.diagonalise_segments();
  builder.solve_postselected_subspace();
  builder.solve_initialised_subspace();
  builder.solve_collapsed_subspace();
  builder.remove_unused_qubits();
  builder.assign_init_post_names();
  builder.assign_remaining_names();
  return builder.unitary_output_circuit();
}

ChoiMixBuilder::ChoiMixBuilder(const ChoiMixTableau& t, CXConfigType cx)
    : ChoiMixBuilder(t, cx, {}, {}) {}

ChoiMixBuilder::ChoiMixBuilder(
    const ChoiMixTableau& t, CXConfigType cx, const qubit_vector_t& inits,
    const qubit_vector_t& posts)
    : in_circ(),
      post_selected(),
      discarded(),
      collapsed(),
      tab(t),
      in_out_permutation(),
      zero_initialised(),
      mix_initialised(),
      out_circ_tp(),
      cx_config(cx),
      unitary_init_names(inits),
      unitary_post_names(posts) {
  for (unsigned i = 0; i < tab.get_n_boundaries(); ++i) {
    ChoiMixTableau::col_key_t key = tab.col_index_.right.at(i);
    if (key.second == ChoiMixTableau::TableauSegment::Input)
      in_circ.add_qubit(key.first);
    else
      out_circ_tp.add_qubit(key.first);
  }
  for (const Qubit& init_q : unitary_init_names) {
    if (tab.col_index_.left.find(ChoiMixTableau::col_key_t{
            init_q, ChoiMixTableau::TableauSegment::Input}) !=
        tab.col_index_.left.end())
      throw std::logic_error(
          "Free qubit name for initialisation conflicts with existing live "
          "input of ChoiMixTableau");
  }
  for (const Qubit& post_q : unitary_post_names) {
    if (tab.col_index_.left.find(ChoiMixTableau::col_key_t{
            post_q, ChoiMixTableau::TableauSegment::Output}) !=
        tab.col_index_.left.end())
      throw std::logic_error(
          "Free qubit name for post-selection conflicts with existing live "
          "output of ChoiMixTableau");
  }
}

ChoiMixTableau ChoiMixBuilder::realised_tableau() const {
  ChoiMixTableau in_tab = circuit_to_cm_tableau(in_circ);
  for (const Qubit& q : post_selected)
    in_tab.post_select(q, ChoiMixTableau::TableauSegment::Output);
  for (const Qubit& q : discarded)
    in_tab.discard_qubit(q, ChoiMixTableau::TableauSegment::Output);
  for (const Qubit& q : collapsed)
    in_tab.collapse_qubit(q, ChoiMixTableau::TableauSegment::Output);
  ChoiMixTableau out_tab = circuit_to_cm_tableau(out_circ_tp.transpose());
  for (const Qubit& q : zero_initialised)
    out_tab.post_select(q, ChoiMixTableau::TableauSegment::Input);
  for (const Qubit& q : mix_initialised)
    out_tab.discard_qubit(q, ChoiMixTableau::TableauSegment::Input);
  qubit_map_t out_in_permutation{};
  using perm_entry = boost::bimap<Qubit, Qubit>::left_const_reference;
  BOOST_FOREACH (perm_entry entry, in_out_permutation.left) {
    out_in_permutation.insert({entry.second, entry.first});
  }
  out_tab.rename_qubits(
      out_in_permutation, ChoiMixTableau::TableauSegment::Output);
  return ChoiMixTableau::compose(ChoiMixTableau::compose(in_tab, tab), out_tab);
}

void ChoiMixBuilder::solve_id_subspace() {
  // Input-first gaussian elimination to solve input-sides of remaining rows
  tab.canonical_column_order(ChoiMixTableau::TableauSegment::Input);
  tab.gaussian_form();

  std::set<unsigned> solved_rows;
  std::set<Qubit> solved_ins, solved_outs;
  for (unsigned r = 0; r < tab.get_n_rows(); ++r) {
    if (solved_rows.find(r) != solved_rows.end()) continue;

    // Look for a row which anticommutes with row r over the inputs
    std::list<unsigned> xcols, zcols;
    for (unsigned c = 0; c < tab.get_n_inputs(); ++c) {
      if (tab.tab_.xmat(r, c)) xcols.push_back(c);
      if (tab.tab_.zmat(r, c)) zcols.push_back(c);
    }
    for (unsigned r2 = r + 1; r2 < tab.get_n_rows(); ++r2) {
      if (solved_rows.find(r2) != solved_rows.end()) continue;
      bool anti = false;
      for (const unsigned c : xcols) anti ^= tab.tab_.zmat(r2, c);
      for (const unsigned c : zcols) anti ^= tab.tab_.xmat(r2, c);
      if (!anti) continue;

      // Found a candidate pair of rows. Because of the Gaussian elimination, it
      // is more likely that the first mismatching qubit is X for r and Z for
      // r2, so favour reducing r2 to Z and r to X
      ChoiMixTableau::row_tensor_t row_r = tab.get_row(r);
      ChoiMixTableau::row_tensor_t row_r2 = tab.get_row(r2);
      std::pair<Circuit, Qubit> in_diag_circ =
          reduce_anticommuting_paulis_to_z_x(
              row_r2.first, row_r.first, cx_config);
      in_circ.append(in_diag_circ.first);
      for (const Command& com : in_diag_circ.first) {
        auto args = com.get_args();
        qubit_vector_t qbs = {args.begin(), args.end()};
        tab.apply_gate(
            com.get_op_ptr()->dagger()->get_type(), qbs,
            ChoiMixTableau::TableauSegment::Input);
      }

      // Since the full rows must commute but they anticommute over the inputs,
      // they must also anticommute over the outputs; we similarly reduce these
      // down to Z and X
      std::pair<Circuit, Qubit> out_diag_circ_dag =
          reduce_anticommuting_paulis_to_z_x(
              row_r2.second, row_r.second, cx_config);
      out_circ_tp.append(out_diag_circ_dag.first.dagger().transpose());
      for (const Command& com : out_diag_circ_dag.first) {
        auto args = com.get_args();
        qubit_vector_t qbs = {args.begin(), args.end()};
        tab.apply_gate(
            com.get_op_ptr()->get_type(), qbs,
            ChoiMixTableau::TableauSegment::Output);
      }

      // Check that rows have been successfully reduced
      row_r = tab.get_row(r);
      row_r2 = tab.get_row(r2);
      if (row_r.first.string.map.size() != 1 ||
          row_r.first.string.map.begin()->second != Pauli::X ||
          row_r.second.string.map.size() != 1 ||
          row_r.second.string.map.begin()->second != Pauli::X ||
          row_r2.first.string.map.size() != 1 ||
          row_r2.first.string.map.begin()->second != Pauli::Z ||
          row_r2.second.string.map.size() != 1 ||
          row_r2.second.string.map.begin()->second != Pauli::Z)
        throw std::logic_error(
            "Unexpected error during identity reduction in ChoiMixTableau "
            "synthesis");
      // Solve phases
      if (row_r.second.coeff == -1.) {
        in_circ.add_op<Qubit>(OpType::Z, {in_diag_circ.second});
        tab.apply_gate(
            OpType::Z, {in_diag_circ.second},
            ChoiMixTableau::TableauSegment::Input);
      }
      if (row_r2.second.coeff == -1.) {
        in_circ.add_op<Qubit>(OpType::X, {in_diag_circ.second});
        tab.apply_gate(
            OpType::X, {in_diag_circ.second},
            ChoiMixTableau::TableauSegment::Input);
      }
      // Connect in permutation
      in_out_permutation.insert(
          {in_diag_circ.second, out_diag_circ_dag.second});
      solved_rows.insert(r);
      solved_rows.insert(r2);

      // Remove these solved qubits from other rows; by commutation of rows, a
      // row contains Z@in_diag_circ.second iff it contains
      // Z@out_diag_circ_dag.second and similarly for X
      unsigned in_c = tab.col_index_.left.at(ChoiMixTableau::col_key_t{
          in_diag_circ.second, ChoiMixTableau::TableauSegment::Input});
      for (unsigned r3 = 0; r3 < tab.get_n_rows(); ++r3) {
        if (r3 != r && tab.tab_.xmat(r3, in_c)) tab.tab_.row_mult(r, r3);
        if (r3 != r2 && tab.tab_.zmat(r3, in_c)) tab.tab_.row_mult(r2, r3);
      }
      solved_ins.insert({in_diag_circ.second});
      solved_outs.insert({out_diag_circ_dag.second});
      break;
    }
  }

  // Remove solved rows and qubits from tableau; since removing rows/columns
  // replaces them with the row/column from the end, remove in reverse order
  for (auto it = solved_rows.rbegin(); it != solved_rows.rend(); ++it)
    tab.remove_row(*it);
  for (auto it = solved_ins.rbegin(); it != solved_ins.rend(); ++it)
    tab.discard_qubit(*it, ChoiMixTableau::TableauSegment::Input);
  for (auto it = solved_outs.rbegin(); it != solved_outs.rend(); ++it)
    tab.discard_qubit(*it, ChoiMixTableau::TableauSegment::Output);
}

// Given a matrix that is already in upper echelon form, use the fact that the
// leading columns are already unique to give column operations that reduce it
// down to identity over the leading columns, eliminating extra swap gates to
// move to the first spaces
static std::vector<std::pair<unsigned, unsigned>>
leading_column_gaussian_col_ops(const MatrixXb& source) {
  std::vector<unsigned> col_list;
  std::set<unsigned> non_leads;
  for (unsigned r = 0; r < source.rows(); ++r) {
    bool leading_found = false;
    for (unsigned c = 0; c < source.cols(); ++c) {
      if (source(r, c)) {
        if (leading_found)
          non_leads.insert(c);
        else {
          leading_found = true;
          col_list.push_back(c);
        }
      }
    }
  }
  for (const unsigned& c : non_leads) col_list.push_back(c);
  MatrixXb reordered = MatrixXb::Zero(source.rows(), col_list.size());
  for (unsigned c = 0; c < col_list.size(); ++c)
    reordered.col(c) = source.col(col_list.at(c));
  std::vector<std::pair<unsigned, unsigned>> reordered_ops =
      gaussian_elimination_col_ops(reordered);
  std::vector<std::pair<unsigned, unsigned>> res;
  for (const std::pair<unsigned, unsigned>& op : reordered_ops)
    res.push_back({col_list.at(op.first), col_list.at(op.second)});
  return res;
}

void ChoiMixBuilder::diagonalise_segments() {
  // Canonicalise tableau
  tab.canonical_column_order(ChoiMixTableau::TableauSegment::Output);
  tab.gaussian_form();

  // Set up diagonalisation tasks
  std::list<std::pair<QubitPauliTensor, Expr>> to_diag_ins, to_diag_outs;
  for (unsigned r = 0; r < tab.get_n_rows(); ++r) {
    ChoiMixTableau::row_tensor_t rten = tab.get_row(r);
    if (!rten.first.string.map.empty()) to_diag_ins.push_back({rten.first, 1.});
    if (!rten.second.string.map.empty())
      to_diag_outs.push_back({rten.second, 1.});
  }
  qubit_vector_t input_qubits = tab.input_qubits();
  std::set<Qubit> diag_ins{input_qubits.begin(), input_qubits.end()};
  Circuit in_diag_circ = mutual_diagonalise(to_diag_ins, diag_ins, cx_config);
  for (const Command& com : in_diag_circ) {
    auto args = com.get_args();
    in_circ.add_op<UnitID>(com.get_op_ptr(), args);
    qubit_vector_t qbs = {args.begin(), args.end()};
    tab.apply_gate(
        com.get_op_ptr()->dagger()->get_type(), qbs,
        ChoiMixTableau::TableauSegment::Input);
  }
  qubit_vector_t output_qubits = tab.output_qubits();
  std::set<Qubit> diag_outs{output_qubits.begin(), output_qubits.end()};
  Circuit out_diag_circ =
      mutual_diagonalise(to_diag_outs, diag_outs, cx_config);
  for (const Command& com : out_diag_circ) {
    auto args = com.get_args();
    out_circ_tp.add_op<UnitID>(com.get_op_ptr()->dagger()->transpose(), args);
    qubit_vector_t qbs = {args.begin(), args.end()};
    tab.apply_gate(
        com.get_op_ptr()->get_type(), qbs,
        ChoiMixTableau::TableauSegment::Output);
  }

  // All rows are diagonalised, so we can just focus on the Z matrix
  if (tab.tab_.xmat != MatrixXb::Zero(tab.get_n_rows(), tab.get_n_boundaries()))
    throw std::logic_error(
        "Diagonalisation in ChoiMixTableau synthesis failed");
}

void ChoiMixBuilder::solve_postselected_subspace() {
  // As column order is currently output first, gaussian form will reveal the
  // post-selected space at the bottom of the tableau and the submatrix of those
  // rows will already be in upper echelon form
  tab.gaussian_form();
  // Reduce them to a minimal set of qubits using CX gates
  unsigned n_postselected = 0;
  for (; n_postselected < tab.get_n_rows(); ++n_postselected) {
    if (!tab.get_row(tab.get_n_rows() - 1 - n_postselected)
             .second.string.map.empty())
      break;
  }
  unsigned n_ins = tab.get_n_inputs();
  unsigned n_outs = tab.get_n_outputs();
  MatrixXb subtableau = tab.tab_.zmat.bottomRightCorner(n_postselected, n_ins);
  std::vector<std::pair<unsigned, unsigned>> col_ops =
      leading_column_gaussian_col_ops(subtableau);
  for (const std::pair<unsigned, unsigned>& op : col_ops) {
    unsigned tab_ctrl_col = n_outs + op.second;
    unsigned tab_trgt_col = n_outs + op.first;
    tab.tab_.apply_CX(tab_ctrl_col, tab_trgt_col);
    ChoiMixTableau::col_key_t ctrl = tab.col_index_.right.at(tab_ctrl_col);
    ChoiMixTableau::col_key_t trgt = tab.col_index_.right.at(tab_trgt_col);
    in_circ.add_op<Qubit>(OpType::CX, {ctrl.first, trgt.first});
  }
  // Postselect rows
  for (unsigned r = 0; r < n_postselected; ++r) {
    unsigned final_row = tab.get_n_rows() - 1;
    ChoiMixTableau::row_tensor_t row = tab.get_row(final_row);
    if (row.second.string.map.size() != 0 || row.first.string.map.size() != 1 ||
        row.first.string.map.begin()->second != Pauli::Z)
      throw std::logic_error(
          "Unexpected error during post-selection identification in "
          "ChoiMixTableau synthesis");
    Qubit post_selected_qb = row.first.string.map.begin()->first;
    // Multiply other rows to remove Z_qb components
    unsigned qb_col = tab.col_index_.left.at(ChoiMixTableau::col_key_t{
        post_selected_qb, ChoiMixTableau::TableauSegment::Input});
    for (unsigned s = 0; s < final_row; ++s)
      if (tab.tab_.zmat(s, qb_col)) tab.tab_.row_mult(final_row, s);
    // Post-select on correct phase
    if (row.second.coeff == -1.)
      in_circ.add_op<Qubit>(OpType::X, {post_selected_qb});
    tab.remove_row(final_row);
    post_selected.insert(post_selected_qb);
    tab.discard_qubit(post_selected_qb, ChoiMixTableau::TableauSegment::Input);
  }
}

void ChoiMixBuilder::solve_initialised_subspace() {
  // Input-first gaussian elimination now sorts the remaining rows into the
  // collapsed subspace followed by the zero-initialised subspace and the
  // collapsed subspace rows are in upper echelon form over the inputs, giving
  // unique leading columns and allowing us to solve them with CXs by
  // column-wise gaussian elimination; same for zero-initialised rows over the
  // outputs
  tab.canonical_column_order(ChoiMixTableau::TableauSegment::Input);
  tab.gaussian_form();

  // Reduce the zero-initialised space to a minimal set of qubits using CX gates
  unsigned n_collapsed = 0;
  for (; n_collapsed < tab.get_n_rows(); ++n_collapsed) {
    if (tab.get_row(n_collapsed).first.string.map.empty()) break;
  }
  unsigned n_ins = tab.get_n_inputs();
  unsigned n_outs = tab.get_n_outputs();
  MatrixXb subtableau =
      tab.tab_.zmat.bottomRightCorner(tab.get_n_rows() - n_collapsed, n_outs);
  std::vector<std::pair<unsigned, unsigned>> col_ops =
      leading_column_gaussian_col_ops(subtableau);
  for (const std::pair<unsigned, unsigned>& op : col_ops) {
    unsigned tab_ctrl_col = n_ins + op.second;
    unsigned tab_trgt_col = n_ins + op.first;
    tab.tab_.apply_CX(tab_ctrl_col, tab_trgt_col);
    ChoiMixTableau::col_key_t ctrl = tab.col_index_.right.at(tab_ctrl_col);
    ChoiMixTableau::col_key_t trgt = tab.col_index_.right.at(tab_trgt_col);
    out_circ_tp.add_op<Qubit>(OpType::CX, {ctrl.first, trgt.first});
  }
  // Initialise rows
  for (unsigned r = tab.get_n_rows(); r-- > n_collapsed;) {
    // r always refers to the final row in the tableau
    ChoiMixTableau::row_tensor_t row = tab.get_row(r);
    if (row.first.string.map.size() != 0 || row.second.string.map.size() != 1 ||
        row.second.string.map.begin()->second != Pauli::Z)
      throw std::logic_error(
          "Unexpected error during initialisation identification in "
          "ChoiMixTableau synthesis");
    Qubit initialised_qb = row.second.string.map.begin()->first;
    // Multiply other rows to remove Z_qb components
    unsigned qb_col = tab.col_index_.left.at(ChoiMixTableau::col_key_t{
        initialised_qb, ChoiMixTableau::TableauSegment::Output});
    for (unsigned s = 0; s < r; ++s)
      if (tab.tab_.zmat(s, qb_col)) tab.tab_.row_mult(r, s);
    // Initialise with correct phase
    if (row.second.coeff == -1.)
      out_circ_tp.add_op<Qubit>(OpType::X, {initialised_qb});
    tab.remove_row(r);
    zero_initialised.insert(initialised_qb);
    tab.discard_qubit(initialised_qb, ChoiMixTableau::TableauSegment::Output);
  }
}

void ChoiMixBuilder::solve_collapsed_subspace() {
  // Solving the initialised subspace will have preserved the upper echelon form
  // of the collapsed subspace; reduce the inputs of the collapsed space to a
  // minimal set of qubits using CX gates
  unsigned n_ins = tab.get_n_inputs();
  unsigned n_outs = tab.get_n_outputs();
  MatrixXb subtableau = tab.tab_.zmat.topLeftCorner(tab.get_n_rows(), n_ins);
  std::vector<std::pair<unsigned, unsigned>> col_ops =
      leading_column_gaussian_col_ops(subtableau);
  for (const std::pair<unsigned, unsigned>& op : col_ops) {
    tab.tab_.apply_CX(op.second, op.first);
    ChoiMixTableau::col_key_t ctrl = tab.col_index_.right.at(op.second);
    ChoiMixTableau::col_key_t trgt = tab.col_index_.right.at(op.first);
    in_circ.add_op<Qubit>(OpType::CX, {ctrl.first, trgt.first});
  }
  // Since row multiplications will unsolve the inputs, we cannot get the output
  // segment into upper echelon form for the same CX-saving trick; instead we
  // accept just removing any qubits that are now unused after solving the
  // initialised subspace
  remove_unused_qubits();
  tab.canonical_column_order(ChoiMixTableau::TableauSegment::Input);
  // Solve the output segment using CX gates
  n_ins = tab.get_n_inputs();
  n_outs = tab.get_n_outputs();
  col_ops = gaussian_elimination_col_ops(
      tab.tab_.zmat.topRightCorner(tab.get_n_rows(), n_outs));
  for (const std::pair<unsigned, unsigned>& op : col_ops) {
    tab.tab_.apply_CX(n_ins + op.second, n_ins + op.first);
    ChoiMixTableau::col_key_t ctrl = tab.col_index_.right.at(n_ins + op.second);
    ChoiMixTableau::col_key_t trgt = tab.col_index_.right.at(n_ins + op.first);
    out_circ_tp.add_op<Qubit>(OpType::CX, {ctrl.first, trgt.first});
  }
  // Connect up and remove rows and columns
  for (unsigned r = tab.get_n_rows(); r-- > 0;) {
    // r refers to the final row
    // Check that row r has been successfully reduced
    ChoiMixTableau::row_tensor_t row_r = tab.get_row(r);
    if (row_r.first.string.map.size() != 1 ||
        row_r.first.string.map.begin()->second != Pauli::Z ||
        row_r.second.string.map.size() != 1 ||
        row_r.second.string.map.begin()->second != Pauli::Z)
      throw std::logic_error(
          "Unexpected error during collapsed subspace reduction in "
          "ChoiMixTableau synthesis");
    Qubit in_q = row_r.first.string.map.begin()->first;
    Qubit out_q = row_r.second.string.map.begin()->first;
    // Solve phase
    if (row_r.second.coeff == -1.) {
      in_circ.add_op<Qubit>(OpType::X, {in_q});
      tab.apply_gate(OpType::X, {in_q}, ChoiMixTableau::TableauSegment::Input);
    }
    // Connect in permutation
    in_out_permutation.insert({in_q, out_q});
    collapsed.insert(in_q);
    tab.remove_row(r);
    tab.discard_qubit(in_q, ChoiMixTableau::TableauSegment::Input);
    tab.discard_qubit(out_q, ChoiMixTableau::TableauSegment::Output);
  }
}

void ChoiMixBuilder::remove_unused_qubits() {
  // Since removing a column replaces it with the last column, remove in reverse
  // order to examine each column exactly once
  for (unsigned c = tab.get_n_boundaries(); c-- > 0;) {
    bool used = false;
    for (unsigned r = 0; r < tab.get_n_rows(); ++r) {
      if (tab.tab_.zmat(r, c) || tab.tab_.xmat(r, c)) {
        used = true;
        break;
      }
    }
    if (used) continue;
    ChoiMixTableau::col_key_t col = tab.col_index_.right.at(c);
    if (col.second == ChoiMixTableau::TableauSegment::Input)
      discarded.insert(col.first);
    else
      mix_initialised.insert(col.first);
    tab.discard_qubit(col.first, col.second);
  }
}

void ChoiMixBuilder::assign_init_post_names() {
  auto it = unitary_post_names.begin();
  for (const Qubit& ps : post_selected) {
    if (it == unitary_post_names.end())
      throw std::logic_error(
          "Not enough additional qubit names for unitary extension of "
          "ChoiMixTableau to safely handle post-selected subspace");
    in_out_permutation.insert({ps, *it});
    ++it;
  }
  unitary_post_names = {it, unitary_post_names.end()};

  it = unitary_init_names.begin();
  for (const Qubit& zi : zero_initialised) {
    if (it == unitary_init_names.end())
      throw std::logic_error(
          "Not enough additional qubit names for unitary extension of "
          "ChoiMixTableau to safely handle initialised subspace");
    in_out_permutation.insert({*it, zi});
    ++it;
  }
  unitary_init_names = {it, unitary_init_names.end()};
}

void ChoiMixBuilder::assign_remaining_names() {
  // Some post-selected or initialised qubits might have already been matched up
  // for unitary synthesis, so we only need to match up the remainder
  std::set<Qubit> unsolved_ins = discarded;
  for (const Qubit& q : post_selected) {
    if (in_out_permutation.left.find(q) == in_out_permutation.left.end())
      unsolved_ins.insert(q);
  }
  std::set<Qubit> unsolved_outs = mix_initialised;
  for (const Qubit& q : zero_initialised) {
    if (in_out_permutation.right.find(q) == in_out_permutation.right.end())
      unsolved_outs.insert(q);
  }
  // If there are more unsolved_ins than unsolved_outs, we want to pad out
  // unsolved_outs with extra names that don't appear as output names of the
  // original tableau; between unsolved_ins and the inputs already in
  // in_out_permutation, there will be at least enough of these
  if (unsolved_ins.size() > unsolved_outs.size()) {
    for (const Qubit& q : unsolved_ins) {
      if (in_out_permutation.right.find(q) == in_out_permutation.right.end()) {
        unsolved_outs.insert(q);
        if (unsolved_ins.size() == unsolved_outs.size()) break;
      }
    }
    if (unsolved_ins.size() > unsolved_outs.size()) {
      using perm_entry = boost::bimap<Qubit, Qubit>::left_const_reference;
      BOOST_FOREACH (perm_entry entry, in_out_permutation.left) {
        if (in_out_permutation.right.find(entry.first) ==
            in_out_permutation.right.end()) {
          unsolved_outs.insert(entry.first);
          if (unsolved_ins.size() == unsolved_outs.size()) break;
        }
      }
    }
  } else if (unsolved_ins.size() < unsolved_outs.size()) {
    for (const Qubit& q : unsolved_outs) {
      if (in_out_permutation.left.find(q) == in_out_permutation.left.end()) {
        unsolved_ins.insert(q);
        if (unsolved_ins.size() == unsolved_outs.size()) break;
      }
    }
    if (unsolved_ins.size() < unsolved_outs.size()) {
      using perm_entry = boost::bimap<Qubit, Qubit>::left_const_reference;
      BOOST_FOREACH (perm_entry entry, in_out_permutation.left) {
        if (in_out_permutation.left.find(entry.second) ==
            in_out_permutation.left.end()) {
          unsolved_ins.insert(entry.second);
          if (unsolved_ins.size() == unsolved_outs.size()) break;
        }
      }
    }
  }
  // Prefer to connect qubits with the same names
  for (auto in_it = unsolved_ins.begin(); in_it != unsolved_ins.end();) {
    auto temp_it = in_it++;
    auto out_it = unsolved_outs.find(*temp_it);
    if (out_it != unsolved_outs.end()) {
      in_out_permutation.insert({*temp_it, *temp_it});
      unsolved_ins.erase(temp_it);
      unsolved_outs.erase(out_it);
    }
  }
  // Pair up remainders; by our earlier padding, they should have the exact same
  // number of elements, so pair them up exactly
  for (const Qubit& in : unsolved_ins) {
    auto it = unsolved_outs.begin();
    in_out_permutation.insert({in, *it});
    unsolved_outs.erase(it);
  }
}

std::pair<Circuit, qubit_map_t> ChoiMixBuilder::output_circuit() {
  if (tab.get_n_rows() != 0 || tab.get_n_boundaries() != 0)
    throw std::logic_error(
        "Unexpected error during ChoiMixTableau synthesis, reached the end "
        "with a non-empty tableau remaining");
  if (!post_selected.empty()) {
    throw std::logic_error(
        "Not yet implemented: post-selection required during ChoiMixTableau "
        "synthesis");
  }
  for (const Qubit& q : discarded) in_circ.qubit_discard(q);
  for (const Qubit& q : collapsed) in_circ.add_op<Qubit>(OpType::Collapse, {q});
  Circuit out_circ(out_circ_tp.all_qubits(), {});
  for (const Qubit& q : zero_initialised) out_circ.qubit_create(q);
  for (const Qubit& q : mix_initialised) {
    out_circ.qubit_create(q);
    out_circ.add_op<Qubit>(OpType::H, {q});
    out_circ.add_op<Qubit>(OpType::Collapse, {q});
  }
  out_circ.append(out_circ_tp.transpose());
  qubit_map_t return_perm;
  unit_map_t append_perm;
  using perm_entry = boost::bimap<Qubit, Qubit>::left_const_reference;
  BOOST_FOREACH (perm_entry entry, in_out_permutation.left) {
    return_perm.insert({entry.second, entry.first});
    append_perm.insert({entry.second, entry.first});
  }
  in_circ.append_with_map(out_circ, append_perm);
  return {in_circ, return_perm};
}

std::pair<Circuit, qubit_map_t> ChoiMixBuilder::unitary_output_circuit() {
  if (tab.get_n_rows() != 0 || tab.get_n_boundaries() != 0)
    throw std::logic_error(
        "Unexpected error during ChoiMixTableau synthesis, reached the end "
        "with a non-empty tableau remaining");
  qubit_map_t return_perm;
  unit_map_t append_perm;
  using perm_entry = boost::bimap<Qubit, Qubit>::left_const_reference;
  BOOST_FOREACH (perm_entry entry, in_out_permutation.left) {
    return_perm.insert({entry.second, entry.first});
    append_perm.insert({entry.second, entry.first});
  }
  in_circ.append_with_map(out_circ_tp.transpose(), append_perm);
  return {in_circ, return_perm};
}

}  // namespace tket
