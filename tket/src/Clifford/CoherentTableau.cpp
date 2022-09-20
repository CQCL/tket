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

#include "CoherentTableau.hpp"

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

CoherentTableau::CoherentTableau(unsigned n) : tab_(id_tab(n)), col_index_() {
  for (unsigned i = 0; i < n; ++i) {
    col_index_.insert({{Qubit(i), TableauSegment::Input}, i});
    col_index_.insert({{Qubit(i), TableauSegment::Output}, n + i});
  }
}

CoherentTableau::CoherentTableau(const qubit_vector_t& qbs)
    : tab_(id_tab(qbs.size())), col_index_() {
  unsigned n = qbs.size();
  unsigned i = 0;
  for (const Qubit& qb : qbs) {
    col_index_.insert({{qb, TableauSegment::Input}, i});
    col_index_.insert({{qb, TableauSegment::Output}, n + i});
    ++i;
  }
}

CoherentTableau::CoherentTableau(
    const MatrixXb& xmat, const MatrixXb& zmat, const VectorXb& phase,
    unsigned n_ins)
    : tab_({}), col_index_() {
  unsigned n_rows = xmat.rows();
  unsigned n_bounds = xmat.cols();
  if ((zmat.cols() != n_bounds) || (zmat.rows() != n_rows) ||
      (phase.size() != n_rows))
    throw std::invalid_argument(
        "Partial tableau requires equally-sized components");
  tab_ = SymplecticTableau(xmat, zmat, phase);
  if (tab_.anticommuting_rows() != MatrixXb::Zero(n_rows, n_rows))
    throw std::invalid_argument("Rows of partial tableau do not commute");
  if (tab_.rank() != n_rows)
    throw std::invalid_argument("Rows of partial tableau are not independent");
  for (unsigned i = 0; i < n_ins; ++i) {
    col_index_.insert({{Qubit(i), TableauSegment::Input}, i});
  }
  for (unsigned i = 0; i < n_bounds - n_ins; ++i) {
    col_index_.insert({{Qubit(i), TableauSegment::Output}, i});
  }
}

CoherentTableau::CoherentTableau(const std::list<row_tensor_t>& rows)
    : tab_({}), col_index_() {
  std::set<Qubit> in_qubits;
  std::set<Qubit> out_qubits;
  for (const row_tensor_t& row : rows) {
    for (const std::pair<const Qubit, Pauli>& qb : row.first.string.map) {
      in_qubits.insert(qb.first);
    }
    for (const std::pair<const Qubit, Pauli>& qb : row.second.string.map) {
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
  for (unsigned r = 0; r < n_rows; ++r) {
    const row_tensor_t& row = rows.at(r);
    for (const std::pair<const Qubit, Pauli>& qb : row.first.string.map) {
      unsigned c = col_index_.at({qb.first, TableauSegment::Input});
      if (qb.second == Pauli::X || qb.second == Pauli::Y) xmat(r, c) = true;
      if (qb.second == Pauli::Z || qb.second == Pauli::Y) zmat(r, c) = true;
    }
    for (const std::pair<const Qubit, Pauli>& qb : row.second.string.map) {
      unsigned c = col_index_.at({qb.first, TableauSegment::Output});
      if (qb.second == Pauli::X || qb.second == Pauli::Y) xmat(r, c) = true;
      if (qb.second == Pauli::Z || qb.second == Pauli::Y) zmat(r, c) = true;
    }
    phase(r) = (row.first.coeff == -1.) ^ (row.second.coeff == -1.);
  }
  tab_ = SymplecticTableau(xmat, zmat, phase);
}

unsigned CoherentTableau::get_n_rows() const { return tab_.get_n_rows(); }

unsigned CoherentTableau::get_n_boundaries() const { return col_index_.size(); }

unsigned CoherentTableau::get_n_inputs() const {
  unsigned n = 0;
  for (const std::pair<const std::pair<Qubit, TableauSegment>, unsigned>&
           entry : col_index_.left) {
    if (entry.first.second == TableauSegment::Input) ++n;
  }
  return n;
}

unsigned CoherentTableau::get_n_outputs() const {
  unsigned n = 0;
  for (const std::pair<const std::pair<Qubit, TableauSegment>, unsigned>&
           entry : col_index_.left) {
    if (entry.first.second == TableauSegment::Output) ++n;
  }
  return n;
}

row_tensor_t CoherentTableau::stab_to_row_tensor(
    const PauliStabiliser& stab) const {
  QubitPauliMap in_qpm, out_qpm;
  double coeff = 1.;
  for (unsigned i = 0; i < stab.string.size(); ++i) {
    Pauli p = stab.string.at(i);
    std::pair<Qubit, TableauSegment> col = col_index_.right.at(i);
    if (p != Pauli::I) {
      if (col.second == TableauSegment::Input)
        in_qpm.insert({col.first, p});
      else
        out_qpm.insert({col.first, p});
    }
  }
  return {
      QubitPauliTensor(in_qpm),
      QubitPauliTensor(out_qpm, stab.coeff ? 1. : -1.)};
}

PauliStabiliser CoherentTableau::row_tensor_to_stab(
    const row_tensor_t& ten) const {
  std::vector<Pauli> ps;
  for (unsigned i = 0; i < col_index.size(); ++i) {
    std::pair<Qubit, TableauSegment> qb = col_index.right.at(i);
    if (qb.second == TableauSegment::Input)
      ps.push_back(ten.first.string.get(qb.first));
    else
      ps.push_back(ten.second.string.get(qb.first));
  }
  return PauliStabiliser(ps, (ten.first.coeff * ten.second.coeff == 1.));
}

row_tensor_t CoherentTableau::get_row(unsigned i) const {
  return stab_to_row_tensor(tab_.get_pauli(i));
}

row_tensor_t CoherentTableau::get_row_product(
    const std::vector<unsigned>& rows) const {
  row_tensor_t result = {{}, {}};
  for (unsigned i = 0; i < get_n_rows(); ++i) {
    row_tensor_t row_i = get_row(i);
    result.first = result.first * row_i.first;
    result.second = result.second * row_t.second;
  }
  return result;
}

void CoherentTableau::apply_S(const Qubit& qb, TableauSegment seg) {
  unsigned col = col_index_.left.at({qb, seg});
  tab_.apply_S(col);
}

void CoherentTableau::apply_V(const Qubit& qb, TableauSegment seg) {
  unsigned col = col_index_.left.at({qb, seg});
  tab_.apply_V(col);
}

void CoherentTableau::apply_CX(
    const Qubit& control, const Qubit& target, TableauSegment seg) {
  unsigned uc = col_index_.left.at({control, seg});
  unsigned ut = col_index_.left.at({target, seg});
  tab_.apply_CX(uc, ut);
}

void CoherentTableau::apply_gate(
    OpType type, const qubit_vector_t& qbs, TableauSegment seg) {
  switch (type) {
    case OpType::Z: {
      apply_S(qbs.at(0), seg);
      apply_S(qbs.at(0), seg);
      break;
    }
    case OpType::X: {
      apply_V(qbs.at(0), seg);
      apply_V(qbs.at(0), seg);
      break;
    }
    case OpType::Y: {
      apply_S(qbs.at(0), seg);
      apply_S(qbs.at(0), seg);
      apply_V(qbs.at(0), seg);
      apply_V(qbs.at(0), seg);
      break;
    }
    case OpType::S: {
      apply_S(qbs.at(0), seg);
      break;
    }
    case OpType::Sdg: {
      apply_S(qbs.at(0), seg);
      apply_S(qbs.at(0), seg);
      break;
    }
    case OpType::V: {
      apply_V(qbs.at(0), seg);
      break;
    }
    case OpType::Vdg: {
      apply_V(qbs.at(0), seg);
      apply_V(qbs.at(0), seg);
      break;
    }
    case OpType::H: {
      apply_S(qbs.at(0), seg);
      apply_V(qbs.at(0), seg);
      apply_S(qbs.at(0), seg);
      break;
    }
    case OpType::CX: {
      apply_CX(qbs.at(0), qbs.at(1), seg);
      break;
    }
    case OpType::CY: {
      apply_S(qbs.at(1), seg);
      apply_S(qbs.at(1), seg);
      apply_S(qbs.at(1), seg);
      apply_CX(qbs.at(0), qbs.at(1), seg);
      apply_S(qbs.at(1), seg);
      break;
    }
    case OpType::CZ: {
      apply_S(qbs.at(1), seg);
      apply_V(qbs.at(1), seg);
      apply_S(qbs.at(1), seg);
      apply_CX(qbs.at(0), qbs.at(1), seg);
      apply_S(qbs.at(1), seg);
      apply_V(qbs.at(1), seg);
      apply_S(qbs.at(1), seg);
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
    case OpType::noop: {
      break;
    }
    default: {
      throw BadOpType(
          "Cannot be applied to a CoherentTableau: not a unitary Clifford gate",
          type);
    }
  }
}

void CoherentTableau::apply_pauli(
    const QubitPauliTensor& pauli, unsigned half_pis, TableauSegment seg) {
  PauliStabiliser ps;
  if (seg == TableauSegment::Input) {
    tr = pauli;
    tr.transpose();
    ps = row_tableau_to_stab({tr, {}});
  } else {
    ps = row_tableau_to_stab({{}, pauli});
  }
  tab_.apply_pauli_gadget(ps, half_pis);
}

void CoherentTableau::post_select(const Qubit& qb, TableauSegment seg) {
  tab_.gaussian_form();
  // If +Z or -Z is a stabilizer, it will appear as the only row containing Z
  // after gaussian elimination Check for the deterministic cases
  unsigned n_rows = get_n_rows();
  unsigned n_cols = get_n_boundaries();
  unsigned col = col_index_.left.at({qb, seg});
  for (unsigned r = 0; r < n_rows; ++r) {
    if (tab_.zmat_(r, col)) {
      bool only_z = true;
      for (unsigned c = 0; c < n_cols; ++c) {
        if ((tab_.xmat_(r, c) || tab_.zmat_(r, c)) && (c != col)) {
          only_z = false;
          break;
        }
      }
      if (!only_z) break;  // Not deterministic
      // From here, we know we are in a deterministic case
      // If deterministically fail, throw an exception
      if (tab_.phase_(r))
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
    if (tab_.xmat_(r, col)) {
      if (x_row) {
        // Already found another row with an X, so combine them
        tab_.row_mult(
            tab_.xmat_.row(*x_row), tab_.zmat_.row(*x_row), tab_.phase_(*x_row),
            tab_.xmat_.row(r), tab_.zmat_.row(r), tab_.phase_(r),
            tab_.xmat_.row(r), tab_.zmat_.row(r), tab_.phase_(r));
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

void CoherentTableau::discard_qubit(const Qubit& qb, TableauSegment seg) {
  // Isolate a single row with an X (if one exists)
  std::optional<unsigned> x_row = std::nullopt;
  for (unsigned r = 0; r < n_rows; ++r) {
    if (tab_.xmat_(r, col)) {
      if (x_row) {
        // Already found another row with an X, so combine them
        tab_.row_mult(
            tab_.xmat_.row(*x_row), tab_.zmat_.row(*x_row), tab_.phase_(*x_row),
            tab_.xmat_.row(r), tab_.zmat_.row(r), tab_.phase_(r),
            tab_.xmat_.row(r), tab_.zmat_.row(r), tab_.phase_(r));
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
  for (unsigned r = 0; r < n_rows; ++r) {
    if (tab_.zmat_(r, col)) {
      if (z_row) {
        // Already found another row with a Z, so combine them
        tab_.row_mult(
            tab_.xmat_.row(*z_row), tab_.zmat_.row(*z_row), tab_.phase_(*z_row),
            tab_.xmat_.row(r), tab_.zmat_.row(r), tab_.phase_(r),
            tab_.xmat_.row(r), tab_.zmat_.row(r), tab_.phase_(r));
      } else {
        // This is the first row with a Z
        // Continue searching the rest to make it unique
        z_row = r;
      }
    }
  }
  if (x_row) {
    // Remove the Z stabilizer
    remove_row(*x_row);
  }
  remove_col(col);
}

void CoherentTableau::remove_row(unsigned row) {
  unsigned n_rows = get_n_rows();
  unsigned n_cols = get_n_boundaries();
  if (n_rows > 1) {
    tab_.xmat_.row(row) = tab_.xmat_.row(n_rows - 1);
    tab_.zmat_.row(row) = tab_.zmat_.row(n_rows - 1);
    tab_.phase_(row) = tab_.phase_(n_rows - 1);
  }
  tab_.xmat_.conservativeResize(n_rows - 1, n_cols);
  tab_.zmat_.conservativeResize(n_rows - 1, n_cols);
  tab_.phase_.conservativeResize(n_rows - 1);
}

void CoherentTableau::remove_col(unsigned col) {
  unsigned n_rows = get_n_rows();
  unsigned n_cols = get_n_boundaries();
  if (n_cols > 1) {
    tab_.xmat_.col(col) = tab_.xmat_.col(n_cols - 1);
    tab_.zmat_.col(col) = tab_.zmat_.col(n_cols - 1);
  }
  tab_.xmat_.conservativeResize(n_rows, n_cols - 1);
  tab_.zmat_.conservativeResize(n_rows, n_cols - 1);
  col_index_.right.erase(col);
  if (n_cols > 1)
    col_index_.right.replace_key(col_index_.right.find(n_cols - 1), col);
}

CoherentTableau CoherentTableau::compose(
    const CoherentTableau& first, const CoherentTableau& second) {
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
  fullx << first.tab_.xmat_, MatrixXb::Zero(f_rows, s_cols),
      MatrixXb::Zero(s_rows, f_cols), second.tab_.xmat_;
  fullz << first.tab_.zmat_, MatrixXb::Zero(f_rows, s_cols),
      MatrixXb::Zero(s_rows, f_cols), second.tab_.zmat_;
  VectorXb fullph(f_rows + s_rows);
  fullph << first.tab_.phase_, second.tab_.phase_;
  CoherentTableau combined(fullx, fullz, fullph, 0);
  // For each connecting pair of qubits, compose via a Bell post-selection
  for (unsigned i = 0; i < f_cols; ++i) {
    std::pair<Qubit, TableauSegment> ind = first_qubits_to_names.right.at(i);
    if (ind.second == TableauSegment::Output) {
      auto found =
          second_qubits_to_names.left.find({ind.first, TableauSegment::Input});
      if (found != second_qubits_to_names.left.end()) {
        // Found a matching pair
        Qubit f_qb(i), s_qb(found->second);
        combined.apply_CX(f_qb, s_qb);
        combined.apply_gate(OpType::H, f_qb);
        combined.post_select(f_qb);
        combined.post_select(s_qb);
      }
    }
  }
  // Rename qubits to original names
  col_index_t new_index;
  for (const std::pair<const std::pair<Qubit, TableauSegment>, unsigned>& col :
       combined.col_index_.left) {
    bool success = false;
    unsigned qb_num = col.first.first.index().at(0);
    auto found = first_qubits_to_names.right.find(qb_num);
    if (found != first_qubits_to_names.right.end()) {
      success = new_index.insert({found->first, col.second});
    } else {
      success = new_index.insert(
          {second_qubits_to_names.right.at(qb_num), col.second});
    }
    if (!success)
      throw std::logic_error(
          "Qubits aliasing after composing two CoherentTableau objects");
  }
  combined.col_index_ = new_index;
  return combined;
}

std::ostream& operator<<(std::ostream& os, const CoherentTableau& tab) {
  for (unsigned i = 0; i < get_n_rows(); ++i) {
    row_tensor_t row = get_row(i);
    os << row.first.repr() << "\t->\t" << row.second.repr();
  }
  return os;
}

void to_json(nlohmann::json& j, const CoherentTableau::TableauSegment& seg) {
  j = (seg == CoherentTableau::TableauSegment::Input) ? "In" : "Out";
}

void from_json(const nlohmann::json& j, CoherentTableau::TableauSegment& seg) {
  const std::string str_seg = j.get<std::string>();
  seg = (str_seg == "In") ? CoherentTableau::TableauSegment::Input
                          : CoherentTableau::TableauSegment::Output;
}

void to_json(nlohmann::json& j, const CoherentTableau& tab) {
  j["tab"] = tab.tab_;
  std::vector<std::pair<Qubit, TableauSegment>> qbs;
  for (unsigned i = 0; i < get_n_boundaries(); ++i) {
    qbs.push_back(tab.col_index_.right.at(i));
  }
  j["qubits"] = qbs;
}

void from_json(const nlohmann::json& j, CoherentTableau& tab) {
  j.at("tab").get_to(tab.tab_);
  std::vector<std::pair<Qubit, TableauSegment>> qbs =
      j.at("qubits").get<std::vector<std::pair<Qubit, TableauSegment>>>();
  if (qbs.size() != tab.tab_.get_n_qubits())
    throw std::invalid_argument(
        "Number of qubits in json CoherentTableau does not match tableau "
        "size.");
  tab.col_index_.clear();
  for (unsigned i = 0; i < qbs.size(); ++i) {
    tab_col_index_.insert({qbs.at(i), i});
  }
}

}  // namespace tket
