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

#include "tket/Clifford/ChoiMixTableau.hpp"

#include <boost/foreach.hpp>

#include "tket/OpType/OpTypeInfo.hpp"

namespace tket {

static SymplecticTableau id_tab(unsigned n) {
  MatrixXb xmat(2 * n, 2 * n);
  xmat << MatrixXb::Identity(n, n), MatrixXb::Identity(n, n),
      MatrixXb::Zero(n, 2 * n);
  MatrixXb zmat(2 * n, 2 * n);
  zmat << MatrixXb::Zero(n, 2 * n), MatrixXb::Identity(n, n),
      MatrixXb::Identity(n, n);
  return SymplecticTableau(xmat, zmat, VectorXb::Zero(2 * n));
}

ChoiMixTableau::ChoiMixTableau(unsigned n) : tab_(id_tab(n)), col_index_() {
  for (unsigned i = 0; i < n; ++i) {
    col_index_.insert({{Qubit(i), TableauSegment::Input}, i});
    col_index_.insert({{Qubit(i), TableauSegment::Output}, n + i});
  }
}

ChoiMixTableau::ChoiMixTableau(const qubit_vector_t& qbs)
    : tab_(id_tab(qbs.size())), col_index_() {
  unsigned n = qbs.size();
  unsigned i = 0;
  for (const Qubit& qb : qbs) {
    col_index_.insert({{qb, TableauSegment::Input}, i});
    col_index_.insert({{qb, TableauSegment::Output}, n + i});
    ++i;
  }
}

ChoiMixTableau::ChoiMixTableau(
    const MatrixXb& xmat, const MatrixXb& zmat, const VectorXb& phase,
    unsigned n_ins)
    : tab_({}), col_index_() {
  unsigned n_rows = xmat.rows();
  unsigned n_bounds = xmat.cols();
  if (n_ins > n_bounds)
    throw std::invalid_argument(
        "Number of inputs of a Choi tableau cannot be larger than the "
        "number of qubits");
  if ((zmat.cols() != n_bounds) || (zmat.rows() != n_rows) ||
      (phase.size() != n_rows))
    throw std::invalid_argument(
        "Choi tableau requires equally-sized components");
  tab_ = SymplecticTableau(xmat, zmat, phase);
  if (tab_.anticommuting_rows() != MatrixXb::Zero(n_rows, n_rows))
    throw std::invalid_argument("Rows of Choi tableau do not commute");
  if (tab_.rank() != n_rows)
    throw std::invalid_argument("Rows of Choi tableau are not independent");
  for (unsigned i = 0; i < n_ins; ++i) {
    col_index_.insert({{Qubit(i), TableauSegment::Input}, i});
  }
  for (unsigned i = 0; i < n_bounds - n_ins; ++i) {
    col_index_.insert({{Qubit(i), TableauSegment::Output}, i});
  }
}

ChoiMixTableau::ChoiMixTableau(const std::list<row_tensor_t>& rows)
    : tab_({}), col_index_() {
  std::set<Qubit> in_qubits;
  std::set<Qubit> out_qubits;
  for (const row_tensor_t& row : rows) {
    for (const std::pair<const Qubit, Pauli>& qb : row.first.string) {
      in_qubits.insert(qb.first);
    }
    for (const std::pair<const Qubit, Pauli>& qb : row.second.string) {
      out_qubits.insert(qb.first);
    }
  }
  unsigned n_rows = rows.size();
  unsigned n_qbs = in_qubits.size() + out_qubits.size();
  unsigned i = 0;
  for (const Qubit& qb : in_qubits) {
    col_index_.insert({{qb, TableauSegment::Input}, i});
    ++i;
  }
  for (const Qubit& qb : out_qubits) {
    col_index_.insert({{qb, TableauSegment::Output}, i});
    ++i;
  }
  MatrixXb xmat = MatrixXb::Zero(n_rows, n_qbs);
  MatrixXb zmat = MatrixXb::Zero(n_rows, n_qbs);
  VectorXb phase = VectorXb::Zero(n_rows);
  unsigned r = 0;
  for (const row_tensor_t& row : rows) {
    unsigned n_ys = 0;
    for (const std::pair<const Qubit, Pauli>& qb : row.first.string) {
      unsigned c =
          col_index_.left.at(col_key_t{qb.first, TableauSegment::Input});
      if (qb.second == Pauli::X || qb.second == Pauli::Y) xmat(r, c) = true;
      if (qb.second == Pauli::Z || qb.second == Pauli::Y) zmat(r, c) = true;
      if (qb.second == Pauli::Y) ++n_ys;
    }
    for (const std::pair<const Qubit, Pauli>& qb : row.second.string) {
      unsigned c =
          col_index_.left.at(col_key_t{qb.first, TableauSegment::Output});
      if (qb.second == Pauli::X || qb.second == Pauli::Y) xmat(r, c) = true;
      if (qb.second == Pauli::Z || qb.second == Pauli::Y) zmat(r, c) = true;
    }
    phase(r) = row.first.is_real_negative() ^ row.second.is_real_negative() ^
               (n_ys % 2 == 1);
    ++r;
  }
  tab_ = SymplecticTableau(xmat, zmat, phase);
}

unsigned ChoiMixTableau::get_n_rows() const { return tab_.get_n_rows(); }

unsigned ChoiMixTableau::get_n_boundaries() const { return col_index_.size(); }

unsigned ChoiMixTableau::get_n_inputs() const { return input_qubits().size(); }

unsigned ChoiMixTableau::get_n_outputs() const {
  return output_qubits().size();
}

qubit_vector_t ChoiMixTableau::input_qubits() const {
  qubit_vector_t ins;
  BOOST_FOREACH (
      tableau_col_index_t::left_const_reference entry, col_index_.left) {
    if (entry.first.second == TableauSegment::Input)
      ins.push_back(entry.first.first);
  }
  return ins;
}

qubit_vector_t ChoiMixTableau::output_qubits() const {
  qubit_vector_t outs;
  BOOST_FOREACH (
      tableau_col_index_t::left_const_reference entry, col_index_.left) {
    if (entry.first.second == TableauSegment::Output)
      outs.push_back(entry.first.first);
  }
  return outs;
}

ChoiMixTableau::row_tensor_t ChoiMixTableau::stab_to_row_tensor(
    const PauliStabiliser& stab) const {
  QubitPauliMap in_qpm, out_qpm;
  for (unsigned i = 0; i < stab.string.size(); ++i) {
    Pauli p = stab.string.at(i);
    col_key_t col = col_index_.right.at(i);
    if (p != Pauli::I) {
      if (col.second == TableauSegment::Input)
        in_qpm.insert({col.first, p});
      else
        out_qpm.insert({col.first, p});
    }
  }
  return {SpPauliStabiliser(in_qpm), SpPauliStabiliser(out_qpm, stab.coeff)};
}

PauliStabiliser ChoiMixTableau::row_tensor_to_stab(
    const row_tensor_t& ten) const {
  std::vector<Pauli> ps;
  for (unsigned i = 0; i < col_index_.size(); ++i) {
    col_key_t qb = col_index_.right.at(i);
    if (qb.second == TableauSegment::Input)
      ps.push_back(ten.first.get(qb.first));
    else
      ps.push_back(ten.second.get(qb.first));
  }
  return PauliStabiliser(ps, (ten.first.coeff + ten.second.coeff) % 4);
}

ChoiMixTableau::row_tensor_t ChoiMixTableau::get_row(unsigned i) const {
  ChoiMixTableau::row_tensor_t res = stab_to_row_tensor(tab_.get_pauli(i));
  res.first.transpose();
  res.second.coeff = (res.first.coeff + res.second.coeff) % 4;
  res.first.coeff = 0;
  return res;
}

ChoiMixTableau::row_tensor_t ChoiMixTableau::get_row_product(
    const std::vector<unsigned>& rows) const {
  row_tensor_t result = {{}, {}};
  for (unsigned i : rows) {
    row_tensor_t row_i = stab_to_row_tensor(tab_.get_pauli(i));
    result.first = result.first * row_i.first;
    result.second = result.second * row_i.second;
  }
  result.first.transpose();
  result.second.coeff = (result.first.coeff + result.second.coeff) % 4;
  result.first.coeff = 0;
  return result;
}

void ChoiMixTableau::apply_S(const Qubit& qb, TableauSegment seg) {
  unsigned col = col_index_.left.at(col_key_t{qb, seg});
  tab_.apply_S(col);
}

void ChoiMixTableau::apply_Z(const Qubit& qb, TableauSegment seg) {
  unsigned col = col_index_.left.at(col_key_t{qb, seg});
  tab_.apply_Z(col);
}

void ChoiMixTableau::apply_V(const Qubit& qb, TableauSegment seg) {
  unsigned col = col_index_.left.at(col_key_t{qb, seg});
  tab_.apply_V(col);
}

void ChoiMixTableau::apply_X(const Qubit& qb, TableauSegment seg) {
  unsigned col = col_index_.left.at(col_key_t{qb, seg});
  tab_.apply_X(col);
}

void ChoiMixTableau::apply_H(const Qubit& qb, TableauSegment seg) {
  unsigned col = col_index_.left.at(col_key_t{qb, seg});
  tab_.apply_H(col);
}

void ChoiMixTableau::apply_CX(
    const Qubit& control, const Qubit& target, TableauSegment seg) {
  unsigned uc = col_index_.left.at(col_key_t{control, seg});
  unsigned ut = col_index_.left.at(col_key_t{target, seg});
  tab_.apply_CX(uc, ut);
}

void ChoiMixTableau::apply_gate(
    OpType type, const qubit_vector_t& qbs, TableauSegment seg) {
  switch (type) {
    case OpType::Z: {
      apply_Z(qbs.at(0), seg);
      break;
    }
    case OpType::X: {
      apply_X(qbs.at(0), seg);
      break;
    }
    case OpType::Y: {
      apply_Z(qbs.at(0), seg);
      apply_X(qbs.at(0), seg);
      break;
    }
    case OpType::S: {
      apply_S(qbs.at(0), seg);
      break;
    }
    case OpType::Sdg: {
      apply_S(qbs.at(0), seg);
      apply_Z(qbs.at(0), seg);
      break;
    }
    case OpType::SX:
    case OpType::V: {
      apply_V(qbs.at(0), seg);
      break;
    }
    case OpType::SXdg:
    case OpType::Vdg: {
      apply_V(qbs.at(0), seg);
      apply_X(qbs.at(0), seg);
      break;
    }
    case OpType::H: {
      apply_H(qbs.at(0), seg);
      break;
    }
    case OpType::CX: {
      apply_CX(qbs.at(0), qbs.at(1), seg);
      break;
    }
    case OpType::CY: {
      if (seg == TableauSegment::Input) {
        apply_S(qbs.at(1), seg);
        apply_CX(qbs.at(0), qbs.at(1), seg);
        apply_S(qbs.at(1), seg);
        apply_Z(qbs.at(1), seg);
      } else {
        apply_S(qbs.at(1), seg);
        apply_Z(qbs.at(1), seg);
        apply_CX(qbs.at(0), qbs.at(1), seg);
        apply_S(qbs.at(1), seg);
      }
      break;
    }
    case OpType::CZ: {
      apply_H(qbs.at(1), seg);
      apply_CX(qbs.at(0), qbs.at(1), seg);
      apply_H(qbs.at(1), seg);
      break;
    }
    case OpType::ZZMax: {
      apply_CX(qbs.at(0), qbs.at(1), seg);
      apply_S(qbs.at(1), seg);
      apply_CX(qbs.at(0), qbs.at(1), seg);
      break;
    }
    case OpType::ECR: {
      if (seg == TableauSegment::Input) {
        apply_X(qbs.at(0), seg);
        apply_S(qbs.at(0), seg);
        apply_V(qbs.at(1), seg);
        apply_CX(qbs.at(0), qbs.at(1), seg);
      } else {
        apply_CX(qbs.at(0), qbs.at(1), seg);
        apply_S(qbs.at(0), seg);
        apply_X(qbs.at(0), seg);
        apply_V(qbs.at(1), seg);
      }
      break;
    }
    case OpType::ISWAPMax: {
      apply_V(qbs.at(0), seg);
      apply_V(qbs.at(1), seg);
      apply_CX(qbs.at(0), qbs.at(1), seg);
      apply_V(qbs.at(0), seg);
      apply_S(qbs.at(1), seg);
      apply_Z(qbs.at(1), seg);
      apply_CX(qbs.at(0), qbs.at(1), seg);
      apply_V(qbs.at(0), seg);
      apply_V(qbs.at(1), seg);
      break;
    }
    case OpType::SWAP: {
      apply_CX(qbs.at(0), qbs.at(1), seg);
      apply_CX(qbs.at(1), qbs.at(0), seg);
      apply_CX(qbs.at(0), qbs.at(1), seg);
      break;
    }
    case OpType::BRIDGE: {
      apply_CX(qbs.at(0), qbs.at(2), seg);
      break;
    }
    case OpType::Phase:
    case OpType::noop: {
      break;
    }
    case OpType::Reset: {
      if (seg == TableauSegment::Input) {
        post_select(qbs.at(0), TableauSegment::Input);
        unsigned col = get_n_boundaries();
        unsigned rows = get_n_rows();
        // reinsert qubit initialised to maximally mixed state (no coherent
        // stabilizers)
        col_index_.insert({{qbs.at(0), TableauSegment::Input}, col});
        tab_.xmat.conservativeResize(rows, col + 1);
        tab_.xmat.col(col) = MatrixXb::Zero(rows, 1);
        tab_.zmat.conservativeResize(rows, col + 1);
        tab_.zmat.col(col) = MatrixXb::Zero(rows, 1);
      } else {
        discard_qubit(qbs.at(0), TableauSegment::Output);
        unsigned col = get_n_boundaries();
        unsigned rows = get_n_rows();
        // reinsert qubit initialised to |0> (add a Z stabilizer)
        col_index_.insert({{qbs.at(0), TableauSegment::Output}, col});
        tab_.xmat.conservativeResize(rows + 1, col + 1);
        tab_.xmat.col(col) = MatrixXb::Zero(rows + 1, 1);
        tab_.xmat.row(rows) = MatrixXb::Zero(1, col + 1);
        tab_.zmat.conservativeResize(rows + 1, col + 1);
        tab_.zmat.col(col) = MatrixXb::Zero(rows + 1, 1);
        tab_.zmat.row(rows) = MatrixXb::Zero(1, col + 1);
        tab_.zmat(rows, col) = true;
        tab_.phase.conservativeResize(rows + 1);
        tab_.phase(rows) = false;
      }
      break;
    }
    case OpType::Collapse: {
      collapse_qubit(qbs.at(0), seg);
      break;
    }
    default: {
      throw BadOpType(
          "Cannot be applied to a ChoiMixTableau: not a unitary Clifford gate",
          type);
    }
  }
}

void ChoiMixTableau::apply_pauli(
    const SpPauliStabiliser& pauli, unsigned half_pis, TableauSegment seg) {
  PauliStabiliser ps;
  if (seg == TableauSegment::Input) {
    SpPauliStabiliser tr = pauli;
    tr.transpose();
    ps = row_tensor_to_stab({tr, {}});
  } else {
    ps = row_tensor_to_stab({{}, pauli});
  }
  tab_.apply_pauli_gadget(ps, half_pis);
}

void ChoiMixTableau::post_select(const Qubit& qb, TableauSegment seg) {
  tab_.gaussian_form();
  // If +Z or -Z is a stabilizer, it will appear as the only row containing Z
  // after gaussian elimination Check for the deterministic cases
  unsigned n_rows = get_n_rows();
  unsigned n_cols = get_n_boundaries();
  unsigned col = col_index_.left.at(col_key_t{qb, seg});
  for (unsigned r = 0; r < n_rows; ++r) {
    if (tab_.zmat(r, col)) {
      bool only_z = true;
      for (unsigned c = 0; c < n_cols; ++c) {
        if ((tab_.xmat(r, c) || tab_.zmat(r, c)) && (c != col)) {
          only_z = false;
          break;
        }
      }
      if (!only_z) break;  // Not deterministic
      // From here, we know we are in a deterministic case
      // If deterministically fail, throw an exception
      if (tab_.phase(r))
        throw std::logic_error(
            "Post-selecting a tableau fails deterministically");
      // Otherwise, we succeed and remove the stabilizer
      remove_row(r);
      remove_col(col);
      return;
    }
  }
  // For the remainder of the method we handle the non-deterministic case
  // Isolate a single row with an X (if one exists)
  std::optional<unsigned> x_row = std::nullopt;
  for (unsigned r = 0; r < n_rows; ++r) {
    if (tab_.xmat(r, col)) {
      if (x_row) {
        // Already found another row with an X, so combine them
        tab_.row_mult(*x_row, r);
      } else {
        // This is the first row with an X
        // Continue searching the rest to make it unique
        x_row = r;
      }
    }
  }
  if (x_row) {
    // Remove the anti-commuting stabilizer
    remove_row(*x_row);
  }
  remove_col(col);
}

void ChoiMixTableau::discard_qubit(const Qubit& qb, TableauSegment seg) {
  unsigned col = col_index_.left.at(col_key_t{qb, seg});
  // Isolate a single row with an X (if one exists)
  std::optional<unsigned> x_row = std::nullopt;
  for (unsigned r = 0; r < get_n_rows(); ++r) {
    if (tab_.xmat(r, col)) {
      if (x_row) {
        // Already found another row with an X, so combine them
        tab_.row_mult(*x_row, r);
      } else {
        // This is the first row with an X
        // Continue searching the rest to make it unique
        x_row = r;
      }
    }
  }
  if (x_row) {
    // Remove the X stabilizer
    remove_row(*x_row);
  }
  // Isolate a single row with a Z (if one exists)
  std::optional<unsigned> z_row = std::nullopt;
  for (unsigned r = 0; r < get_n_rows(); ++r) {
    if (tab_.zmat(r, col)) {
      if (z_row) {
        // Already found another row with a Z, so combine them
        tab_.row_mult(*z_row, r);
      } else {
        // This is the first row with a Z
        // Continue searching the rest to make it unique
        z_row = r;
      }
    }
  }
  if (z_row) {
    // Remove the Z stabilizer
    remove_row(*z_row);
  }
  remove_col(col);
}

void ChoiMixTableau::collapse_qubit(const Qubit& qb, TableauSegment seg) {
  unsigned col = col_index_.left.at(col_key_t{qb, seg});
  // Isolate a single row with an X (if one exists)
  std::optional<unsigned> x_row = std::nullopt;
  for (unsigned r = 0; r < get_n_rows(); ++r) {
    if (tab_.xmat(r, col)) {
      if (x_row) {
        // Already found another row with an X, so combine them
        tab_.row_mult(*x_row, r);
      } else {
        // This is the first row with an X
        // Continue searching the rest to make it unique
        x_row = r;
      }
    }
  }
  if (x_row) {
    // Remove the X stabilizer
    remove_row(*x_row);
  }
}

void ChoiMixTableau::add_qubit(const Qubit& qb, TableauSegment seg) {
  unsigned n_rows = get_n_rows();
  unsigned n_cols = get_n_boundaries();
  auto inserted = col_index_.insert({{qb, seg}, n_cols});
  if (!inserted.second)
    throw std::logic_error(
        "Could not add qubit to ChoiMixTableau, it already exists");
  tab_.xmat.conservativeResize(n_rows, n_cols + 1);
  tab_.xmat.block(0, n_cols, n_rows, 1) = MatrixXb::Zero(n_rows, 1);
  tab_.zmat.conservativeResize(n_rows, n_cols + 1);
  tab_.zmat.block(0, n_cols, n_rows, 1) = MatrixXb::Zero(n_rows, 1);
}

void ChoiMixTableau::remove_row(unsigned row) {
  if (row >= get_n_rows())
    throw std::invalid_argument(
        "Cannot remove row " + std::to_string(row) + " from tableau with " +
        std::to_string(get_n_rows()) + " rows");
  unsigned n_rows = get_n_rows();
  unsigned n_cols = get_n_boundaries();
  if (row < n_rows - 1) {
    tab_.xmat.row(row) = tab_.xmat.row(n_rows - 1);
    tab_.zmat.row(row) = tab_.zmat.row(n_rows - 1);
    tab_.phase(row) = tab_.phase(n_rows - 1);
  }
  tab_.xmat.conservativeResize(n_rows - 1, n_cols);
  tab_.zmat.conservativeResize(n_rows - 1, n_cols);
  tab_.phase.conservativeResize(n_rows - 1);
}

void ChoiMixTableau::remove_col(unsigned col) {
  if (col >= get_n_boundaries())
    throw std::invalid_argument(
        "Cannot remove column " + std::to_string(col) + " from tableau with " +
        std::to_string(get_n_boundaries()) + " columns");
  unsigned n_rows = get_n_rows();
  unsigned n_cols = get_n_boundaries();
  if (col < n_cols - 1) {
    tab_.xmat.col(col) = tab_.xmat.col(n_cols - 1);
    tab_.zmat.col(col) = tab_.zmat.col(n_cols - 1);
  }
  tab_.xmat.conservativeResize(n_rows, n_cols - 1);
  tab_.zmat.conservativeResize(n_rows, n_cols - 1);
  col_index_.right.erase(col);
  if (col < n_cols - 1) {
    tableau_col_index_t::right_iterator it = col_index_.right.find(n_cols - 1);
    col_key_t last = it->second;
    col_index_.right.erase(it);
    col_index_.insert({last, col});
  }
}

void ChoiMixTableau::canonical_column_order(TableauSegment first) {
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
  unsigned n_rows = get_n_rows();
  MatrixXb xmat = MatrixXb::Zero(n_rows, i);
  MatrixXb zmat = MatrixXb::Zero(n_rows, i);
  for (unsigned j = 0; j < i; ++j) {
    col_key_t key = new_index.right.at(j);
    unsigned c = col_index_.left.at(key);
    xmat.col(j) = tab_.xmat.col(c);
    zmat.col(j) = tab_.zmat.col(c);
  }
  tab_ = SymplecticTableau(xmat, zmat, tab_.phase);
  col_index_ = new_index;
}

void ChoiMixTableau::gaussian_form() { tab_.gaussian_form(); }

void ChoiMixTableau::rename_qubits(
    const qubit_map_t& qmap, TableauSegment seg) {
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

ChoiMixTableau ChoiMixTableau::compose(
    const ChoiMixTableau& first, const ChoiMixTableau& second) {
  // Merge tableau into a single one with only output qubits with default
  // indexing
  tableau_col_index_t first_qubits_to_names;
  tableau_col_index_t second_qubits_to_names;
  unsigned f_rows = first.get_n_rows();
  unsigned f_cols = first.get_n_boundaries();
  unsigned s_rows = second.get_n_rows();
  unsigned s_cols = second.get_n_boundaries();
  for (unsigned i = 0; i < f_cols; ++i) {
    first_qubits_to_names.insert({first.col_index_.right.at(i), i});
  }
  for (unsigned i = 0; i < s_cols; ++i) {
    second_qubits_to_names.insert({second.col_index_.right.at(i), i + f_cols});
  }
  MatrixXb fullx(f_rows + s_rows, f_cols + s_cols),
      fullz(f_rows + s_rows, f_cols + s_cols);
  fullx << first.tab_.xmat, MatrixXb::Zero(f_rows, s_cols),
      MatrixXb::Zero(s_rows, f_cols), second.tab_.xmat;
  fullz << first.tab_.zmat, MatrixXb::Zero(f_rows, s_cols),
      MatrixXb::Zero(s_rows, f_cols), second.tab_.zmat;
  VectorXb fullph(f_rows + s_rows);
  fullph << first.tab_.phase, second.tab_.phase;
  ChoiMixTableau combined(fullx, fullz, fullph, 0);
  // For each connecting pair of qubits, compose via a Bell post-selection
  for (unsigned i = 0; i < f_cols; ++i) {
    col_key_t ind = first_qubits_to_names.right.at(i);
    if (ind.second == TableauSegment::Output) {
      auto found = second_qubits_to_names.left.find(
          col_key_t{ind.first, TableauSegment::Input});
      if (found != second_qubits_to_names.left.end()) {
        // Found a matching pair
        Qubit f_qb(i), s_qb(found->second);
        combined.apply_CX(f_qb, s_qb);
        combined.apply_gate(OpType::H, {f_qb});
        combined.post_select(f_qb);
        combined.post_select(s_qb);
      }
    }
  }
  // Rename qubits to original names
  tableau_col_index_t new_index;
  BOOST_FOREACH (
      tableau_col_index_t::left_const_reference col, combined.col_index_.left) {
    bool success = false;
    unsigned qb_num = col.first.first.index().at(0);
    auto found = first_qubits_to_names.right.find(qb_num);
    if (found != first_qubits_to_names.right.end()) {
      success = new_index.insert({found->second, col.second}).second;
    } else {
      success =
          new_index
              .insert({second_qubits_to_names.right.at(qb_num), col.second})
              .second;
    }
    if (!success)
      throw std::logic_error(
          "Qubits aliasing after composing two ChoiMixTableau objects");
  }
  combined.col_index_ = new_index;
  return combined;
}

std::ostream& operator<<(std::ostream& os, const ChoiMixTableau& tab) {
  for (unsigned i = 0; i < tab.get_n_rows(); ++i) {
    ChoiMixTableau::row_tensor_t row = tab.get_row(i);
    os << row.first.to_str() << "\t->\t" << row.second.to_str() << std::endl;
  }
  return os;
}

bool ChoiMixTableau::operator==(const ChoiMixTableau& other) const {
  return (col_index_ == other.col_index_) && (tab_ == other.tab_);
}

void to_json(nlohmann::json& j, const ChoiMixTableau::TableauSegment& seg) {
  j = (seg == ChoiMixTableau::TableauSegment::Input) ? "In" : "Out";
}

void from_json(const nlohmann::json& j, ChoiMixTableau::TableauSegment& seg) {
  const std::string str_seg = j.get<std::string>();
  seg = (str_seg == "In") ? ChoiMixTableau::TableauSegment::Input
                          : ChoiMixTableau::TableauSegment::Output;
}

void to_json(nlohmann::json& j, const ChoiMixTableau& tab) {
  j["tab"] = tab.tab_;
  std::vector<ChoiMixTableau::col_key_t> qbs;
  for (unsigned i = 0; i < tab.get_n_boundaries(); ++i) {
    qbs.push_back(tab.col_index_.right.at(i));
  }
  j["qubits"] = qbs;
}

void from_json(const nlohmann::json& j, ChoiMixTableau& tab) {
  j.at("tab").get_to(tab.tab_);
  std::vector<ChoiMixTableau::col_key_t> qbs =
      j.at("qubits").get<std::vector<ChoiMixTableau::col_key_t>>();
  if (qbs.size() != tab.tab_.get_n_qubits())
    throw std::invalid_argument(
        "Number of qubits in json ChoiMixTableau does not match tableau "
        "size.");
  tab.col_index_.clear();
  for (unsigned i = 0; i < qbs.size(); ++i) {
    tab.col_index_.insert({qbs.at(i), i});
  }
}

}  // namespace tket
