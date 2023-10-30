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

#include "tket/Clifford/UnitaryTableau.hpp"

#include <stdexcept>

#include "tkassert/Assert.hpp"
#include "tket/Gate/OpPtrFunctions.hpp"
#include "tket/OpType/OpTypeInfo.hpp"
#include "tket/Ops/Op.hpp"

namespace tket {

UnitaryTableau::UnitaryTableau(unsigned n) : tab_({}) {
  MatrixXb xmat(2 * n, n);
  xmat << MatrixXb::Identity(n, n), MatrixXb::Zero(n, n);
  MatrixXb zmat(2 * n, n);
  zmat << MatrixXb::Zero(n, n), MatrixXb::Identity(n, n);
  tab_ = SymplecticTableau(xmat, zmat, VectorXb::Zero(2 * n));
  qubits_ = boost::bimap<Qubit, unsigned>();
  for (unsigned i = 0; i < n; ++i) {
    qubits_.insert({Qubit(i), i});
  }
}

UnitaryTableau::UnitaryTableau(const qubit_vector_t& qbs)
    : UnitaryTableau(qbs.size()) {
  qubits_ = boost::bimap<Qubit, unsigned>();
  for (unsigned i = 0; i < qbs.size(); ++i) {
    qubits_.insert({qbs[i], i});
  }
}

UnitaryTableau::UnitaryTableau(
    const MatrixXb& xx, const MatrixXb& xz, const VectorXb& xph,
    const MatrixXb& zx, const MatrixXb& zz, const VectorXb& zph)
    : tab_({}) {
  unsigned n_qubits = xx.rows();
  if ((xx.cols() != n_qubits) || (xz.rows() != n_qubits) ||
      (xz.cols() != n_qubits) || (xph.size() != n_qubits) ||
      (zx.rows() != n_qubits) || (zx.cols() != n_qubits) ||
      (zz.rows() != n_qubits) || (zz.cols() != n_qubits) ||
      (zph.size() != n_qubits))
    throw std::invalid_argument(
        "Unitary tableau requires equally-sized square matrices and vectors");
  MatrixXb xmat(2 * n_qubits, n_qubits);
  xmat << xx, zx;
  MatrixXb zmat(2 * n_qubits, n_qubits);
  zmat << xz, zz;
  VectorXb phase(2 * n_qubits);
  phase << xph, zph;
  tab_ = SymplecticTableau(xmat, zmat, phase);
  MatrixXb expected_anticommutes(2 * n_qubits, 2 * n_qubits);
  expected_anticommutes << MatrixXb::Zero(n_qubits, n_qubits),
      MatrixXb::Identity(n_qubits, n_qubits),
      MatrixXb::Identity(n_qubits, n_qubits),
      MatrixXb::Zero(n_qubits, n_qubits);
  if (tab_.anticommuting_rows() != expected_anticommutes)
    throw std::invalid_argument(
        "Rows of tableau do not (anti-)commute as expected for UnitaryTableau");
  if (tab_.rank() != 2 * n_qubits)
    throw std::invalid_argument(
        "Rows of UnitaryTableau are not linearly independent");
  qubits_ = boost::bimap<Qubit, unsigned>();
  for (unsigned i = 0; i < n_qubits; ++i) {
    qubits_.insert({Qubit(i), i});
  }
}

UnitaryTableau::UnitaryTableau() : UnitaryTableau(0) {}

SpPauliStabiliser UnitaryTableau::get_xrow(const Qubit& qb) const {
  unsigned uq = qubits_.left.at(qb);
  PauliStabiliser stab = tab_.get_pauli(uq);
  QubitPauliMap qpm;
  for (unsigned i = 0; i < qubits_.size(); ++i) {
    qpm.insert({qubits_.right.at(i), stab.get(i)});
  }
  return SpPauliStabiliser(qpm, stab.coeff);
}

SpPauliStabiliser UnitaryTableau::get_zrow(const Qubit& qb) const {
  unsigned uqb = qubits_.left.at(qb);
  PauliStabiliser stab = tab_.get_pauli(uqb + qubits_.size());
  QubitPauliMap qpm;
  for (unsigned i = 0; i < qubits_.size(); ++i) {
    qpm.insert({qubits_.right.at(i), stab.get(i)});
  }
  return SpPauliStabiliser(qpm, stab.coeff);
}

SpPauliStabiliser UnitaryTableau::get_row_product(
    const SpPauliStabiliser& qpt) const {
  SpPauliStabiliser result({}, qpt.coeff);
  for (const std::pair<const Qubit, Pauli>& p : qpt.string) {
    auto qks_it = qubits_.left.find(p.first);
    if (qks_it == qubits_.left.end()) {
      // Acts as identity on p.first
      result = result * SpPauliStabiliser(p.first, p.second);
    } else {
      switch (p.second) {
        case Pauli::I: {
          break;
        }
        case Pauli::X: {
          result = result * get_xrow(p.first);
          break;
        }
        case Pauli::Y: {
          // Y = iXZ
          result = result * get_xrow(p.first);
          result = result * get_zrow(p.first);
          result.coeff += 1;
          break;
        }
        case Pauli::Z: {
          result = result * get_zrow(p.first);
          break;
        }
      }
    }
  }
  return result;
}

std::set<Qubit> UnitaryTableau::get_qubits() const {
  std::set<Qubit> result;
  for (boost::bimap<Qubit, unsigned>::const_iterator iter = qubits_.begin(),
                                                     iend = qubits_.end();
       iter != iend; ++iter) {
    result.insert(iter->left);
  }
  return result;
}

void UnitaryTableau::apply_S_at_front(const Qubit& qb) {
  unsigned uqb = qubits_.left.at(qb);
  tab_.row_mult(uqb + qubits_.size(), uqb, -i_);
}

void UnitaryTableau::apply_S_at_end(const Qubit& qb) {
  unsigned uqb = qubits_.left.at(qb);
  tab_.apply_S(uqb);
}

void UnitaryTableau::apply_Z_at_front(const Qubit& qb) {
  unsigned uqb = qubits_.left.at(qb);
  tab_.phase(uqb) = !tab_.phase(uqb);
}

void UnitaryTableau::apply_Z_at_end(const Qubit& qb) {
  unsigned uqb = qubits_.left.at(qb);
  tab_.apply_Z(uqb);
}

void UnitaryTableau::apply_V_at_front(const Qubit& qb) {
  unsigned uqb = qubits_.left.at(qb);
  tab_.row_mult(uqb, uqb + qubits_.size(), -i_);
}

void UnitaryTableau::apply_V_at_end(const Qubit& qb) {
  unsigned uqb = qubits_.left.at(qb);
  tab_.apply_V(uqb);
}

void UnitaryTableau::apply_X_at_front(const Qubit& qb) {
  unsigned uqb = qubits_.left.at(qb);
  tab_.phase(uqb + qubits_.size()) = !tab_.phase(uqb + qubits_.size());
}

void UnitaryTableau::apply_X_at_end(const Qubit& qb) {
  unsigned uqb = qubits_.left.at(qb);
  tab_.apply_X(uqb);
}

void UnitaryTableau::apply_H_at_front(const Qubit& qb) {
  unsigned uqb = qubits_.left.at(qb);
  unsigned n_qubits = qubits_.size();
  bool temp = tab_.phase(uqb);
  tab_.phase(uqb) = tab_.phase(uqb + n_qubits);
  tab_.phase(uqb + n_qubits) = temp;
  tab_.xmat.row(uqb).swap(tab_.xmat.row(uqb + n_qubits));
  tab_.zmat.row(uqb).swap(tab_.zmat.row(uqb + n_qubits));
}

void UnitaryTableau::apply_H_at_end(const Qubit& qb) {
  unsigned uqb = qubits_.left.at(qb);
  tab_.apply_H(uqb);
}

void UnitaryTableau::apply_CX_at_front(
    const Qubit& control, const Qubit& target) {
  unsigned uc = qubits_.left.at(control);
  unsigned ut = qubits_.left.at(target);
  tab_.row_mult(ut, uc, 1.);
  tab_.row_mult(uc + qubits_.size(), ut + qubits_.size());
}

void UnitaryTableau::apply_CX_at_end(
    const Qubit& control, const Qubit& target) {
  unsigned uc = qubits_.left.at(control);
  unsigned ut = qubits_.left.at(target);
  tab_.apply_CX(uc, ut);
}

void UnitaryTableau::apply_gate_at_front(
    OpType type, const qubit_vector_t& qbs) {
  switch (type) {
    case OpType::Z: {
      apply_Z_at_front(qbs.at(0));
      break;
    }
    case OpType::X: {
      apply_X_at_front(qbs.at(0));
      break;
    }
    case OpType::Y: {
      apply_Z_at_front(qbs.at(0));
      apply_X_at_front(qbs.at(0));
      break;
    }
    case OpType::S: {
      apply_S_at_front(qbs.at(0));
      break;
    }
    case OpType::Sdg: {
      apply_S_at_front(qbs.at(0));
      apply_Z_at_front(qbs.at(0));
      break;
    }
    case OpType::V:
    case OpType::SX: {
      apply_V_at_front(qbs.at(0));
      break;
    }
    case OpType::Vdg:
    case OpType::SXdg: {
      apply_V_at_front(qbs.at(0));
      apply_X_at_front(qbs.at(0));
      break;
    }
    case OpType::H: {
      apply_H_at_front(qbs.at(0));
      break;
    }
    case OpType::CX: {
      apply_CX_at_front(qbs.at(0), qbs.at(1));
      break;
    }
    case OpType::CY: {
      apply_S_at_front(qbs.at(1));
      apply_CX_at_front(qbs.at(0), qbs.at(1));
      apply_S_at_front(qbs.at(1));
      apply_Z_at_front(qbs.at(1));
      break;
    }
    case OpType::CZ: {
      apply_H_at_front(qbs.at(1));
      apply_CX_at_front(qbs.at(0), qbs.at(1));
      apply_H_at_front(qbs.at(1));
      break;
    }
    case OpType::SWAP: {
      apply_CX_at_front(qbs.at(0), qbs.at(1));
      apply_CX_at_front(qbs.at(1), qbs.at(0));
      apply_CX_at_front(qbs.at(0), qbs.at(1));
      break;
    }
    case OpType::BRIDGE: {
      apply_CX_at_front(qbs.at(0), qbs.at(2));
      break;
    }
    case OpType::ZZMax: {
      apply_H_at_front(qbs.at(1));
      apply_S_at_front(qbs.at(0));
      apply_V_at_front(qbs.at(1));
      apply_CX_at_front(qbs.at(0), qbs.at(1));
      apply_H_at_front(qbs.at(1));
      break;
    }
    case OpType::ECR: {
      apply_CX_at_front(qbs.at(0), qbs.at(1));
      apply_X_at_front(qbs.at(0));
      apply_S_at_front(qbs.at(0));
      apply_V_at_front(qbs.at(1));
      apply_X_at_front(qbs.at(1));
      break;
    }
    case OpType::ISWAPMax: {
      apply_V_at_front(qbs.at(0));
      apply_V_at_front(qbs.at(1));
      apply_CX_at_front(qbs.at(0), qbs.at(1));
      apply_V_at_front(qbs.at(0));
      apply_S_at_front(qbs.at(1));
      apply_Z_at_front(qbs.at(1));
      apply_CX_at_front(qbs.at(0), qbs.at(1));
      apply_V_at_front(qbs.at(0));
      apply_V_at_front(qbs.at(1));
      break;
    }
    case OpType::noop:
    case OpType::Phase: {
      break;
    }
    default: {
      throw BadOpType(
          "Cannot be applied to a UnitaryTableau: not a Clifford gate", type);
    }
  }
}

void UnitaryTableau::apply_gate_at_end(OpType type, const qubit_vector_t& qbs) {
  std::vector<unsigned> uqbs;
  for (const Qubit& q : qbs) {
    uqbs.push_back(qubits_.left.at(q));
  }
  tab_.apply_gate(type, uqbs);
}

void UnitaryTableau::apply_pauli_at_front(
    const SpPauliStabiliser& pauli, unsigned half_pis) {
  return apply_pauli_at_end(get_row_product(pauli), half_pis);
}

void UnitaryTableau::apply_pauli_at_end(
    const SpPauliStabiliser& pauli, unsigned half_pis) {
  std::vector<Pauli> string(qubits_.size(), Pauli::I);
  for (const std::pair<const Qubit, Pauli>& pair : pauli.string) {
    unsigned uqb = qubits_.left.at(pair.first);
    string.at(uqb) = pair.second;
  }
  tab_.apply_pauli_gadget(PauliStabiliser(string, pauli.coeff), half_pis);
}

UnitaryTableau UnitaryTableau::compose(
    const UnitaryTableau& first, const UnitaryTableau& second) {
  std::set<Qubit> qbs = first.get_qubits();
  for (const Qubit& q : second.get_qubits()) {
    qbs.insert(q);
  }
  UnitaryTableau result = UnitaryTableau({});
  const unsigned nqb = qbs.size();

  std::vector<SpPauliStabiliser> rows;

  unsigned qir = 0;
  for (const Qubit& qi : qbs) {
    auto qif_it = first.qubits_.left.find(qi);
    if (qif_it == first.qubits_.left.end()) {
      // First acts as identity on qi, so just take effect of second
      rows.push_back(second.get_xrow(qi));
    } else {
      // Sum rows of second according to entries of first
      SpPauliStabiliser fxrow = first.get_xrow(qi);
      SpPauliStabiliser rxrow = second.get_row_product(fxrow);
      rows.push_back(rxrow);
    }

    result.qubits_.insert({qi, qir});
    ++qir;
  }

  // Do the same for the Z rows
  for (const Qubit& qi : qbs) {
    auto qif_it = first.qubits_.left.find(qi);
    if (qif_it == first.qubits_.left.end()) {
      // First acts as identity on qi, so just take effect of second
      rows.push_back(second.get_zrow(qi));
    } else {
      // Sum rows of second according to entries of first
      SpPauliStabiliser fzrow = first.get_zrow(qi);
      SpPauliStabiliser rzrow = second.get_row_product(fzrow);
      rows.push_back(rzrow);
    }
  }

  // Combine row lists and convert to PauliStabilisers
  PauliStabiliserVec all_rows;
  for (const SpPauliStabiliser& row : rows) {
    std::vector<Pauli> ps(nqb, Pauli::I);
    for (const std::pair<const Qubit, Pauli>& p : row.string) {
      unsigned q = result.qubits_.left.at(p.first);
      ps[q] = p.second;
    }
    all_rows.push_back(PauliStabiliser(ps, row.coeff));
  }

  result.tab_ = SymplecticTableau(all_rows);

  return result;
}

static const std::map<
    std::pair<BoolPauli, BoolPauli>, std::pair<BoolPauli, BoolPauli>>&
invert_cell_map() {
  static const auto inv_map = []() {
    const BoolPauli I = {false, false};
    const BoolPauli X = {true, false};
    const BoolPauli Y = {true, true};
    const BoolPauli Z = {false, true};
    return std::map<
        std::pair<BoolPauli, BoolPauli>, std::pair<BoolPauli, BoolPauli>>{
        {{I, I}, {I, I}}, {{I, X}, {I, X}}, {{I, Y}, {X, X}}, {{I, Z}, {X, I}},
        {{X, I}, {I, Z}}, {{X, X}, {I, Y}}, {{X, Y}, {X, Y}}, {{X, Z}, {X, Z}},
        {{Y, I}, {Z, Z}}, {{Y, X}, {Z, Y}}, {{Y, Y}, {Y, Y}}, {{Y, Z}, {Y, Z}},
        {{Z, I}, {Z, I}}, {{Z, X}, {Z, X}}, {{Z, Y}, {Y, X}}, {{Z, Z}, {Y, I}},
    };
  }();
  return inv_map;
}

UnitaryTableau UnitaryTableau::dagger() const {
  // This method is following Craig Gidney's tableau inversion method
  // https://algassert.com/post/2002
  unsigned nqb = qubits_.size();
  MatrixXb dxx = MatrixXb::Zero(nqb, nqb);
  MatrixXb dxz = MatrixXb::Zero(nqb, nqb);
  VectorXb dxph = VectorXb::Zero(nqb);
  MatrixXb dzx = MatrixXb::Zero(nqb, nqb);
  MatrixXb dzz = MatrixXb::Zero(nqb, nqb);
  VectorXb dzph = VectorXb::Zero(nqb);
  for (unsigned i = 0; i < nqb; ++i) {
    for (unsigned j = 0; j < nqb; ++j) {
      // Take effect of some input on some output and invert
      auto inv_cell = invert_cell_map().at(
          {BoolPauli{tab_.xmat(i, j), tab_.zmat(i, j)},
           BoolPauli{tab_.xmat(i + nqb, j), tab_.zmat(i + nqb, j)}});
      // Transpose tableau and fill in cell
      dxx(j, i) = inv_cell.first.x;
      dxz(j, i) = inv_cell.first.z;
      dzx(j, i) = inv_cell.second.x;
      dzz(j, i) = inv_cell.second.z;
    }
  }

  UnitaryTableau dag(dxx, dxz, dxph, dzx, dzz, dzph);
  dag.qubits_ = qubits_;

  // Correct phases
  for (unsigned i = 0; i < nqb; ++i) {
    SpPauliStabiliser xr = dag.get_xrow(qubits_.right.at(i));
    dag.tab_.phase_(i) = get_row_product(xr).is_real_negative();
    SpPauliStabiliser zr = dag.get_zrow(qubits_.right.at(i));
    dag.tab_.phase_(i + nqb) = get_row_product(zr).is_real_negative();
  }

  return dag;
}

UnitaryTableau UnitaryTableau::transpose() const {
  return dagger().conjugate();
}

UnitaryTableau UnitaryTableau::conjugate() const {
  UnitaryTableau conj(0);
  conj.tab_ = tab_.conjugate();
  conj.qubits_ = qubits_;
  return conj;
}

std::ostream& operator<<(std::ostream& os, const UnitaryTableau& tab) {
  unsigned nqs = tab.qubits_.size();
  for (unsigned i = 0; i < nqs; ++i) {
    Qubit qi = tab.qubits_.right.at(i);
    os << "X@" << qi.repr() << "\t->\t" << tab.tab_.xmat.row(i) << "   "
       << tab.tab_.zmat.row(i) << "   " << tab.tab_.phase(i) << std::endl;
  }
  os << "--" << std::endl;
  for (unsigned i = 0; i < nqs; ++i) {
    Qubit qi = tab.qubits_.right.at(i);
    os << "Z@" << qi.repr() << "\t->\t" << tab.tab_.xmat.row(i + nqs) << "   "
       << tab.tab_.zmat.row(i + nqs) << "   " << tab.tab_.phase(i + nqs)
       << std::endl;
  }
  return os;
}

bool UnitaryTableau::operator==(const UnitaryTableau& other) const {
  if (get_qubits() != other.get_qubits()) return false;

  unsigned nq = qubits_.size();

  for (unsigned i = 0; i < nq; ++i) {
    Qubit qi = qubits_.right.at(i);
    unsigned oi = other.qubits_.left.at(qi);
    for (unsigned j = 0; j < nq; ++j) {
      Qubit qj = qubits_.right.at(j);
      unsigned oj = other.qubits_.left.at(qj);
      if (tab_.xmat(i, j) != other.tab_.xmat(oi, oj)) return false;
      if (tab_.zmat(i, j) != other.tab_.zmat(oi, oj)) return false;
      if (tab_.xmat(i + nq, j) != other.tab_.xmat(oi + nq, oj)) return false;
      if (tab_.zmat(i + nq, j) != other.tab_.zmat(oi + nq, oj)) return false;
    }
    if (tab_.phase(i) != other.tab_.phase(oi)) return false;
    if (tab_.phase(i + nq) != other.tab_.phase(oi + nq)) return false;
  }

  return true;
}

void to_json(nlohmann::json& j, const UnitaryTableau& tab) {
  j["tab"] = tab.tab_;
  qubit_vector_t qbs;
  for (unsigned i = 0; i < tab.qubits_.size(); ++i) {
    qbs.push_back(tab.qubits_.right.at(i));
  }
  j["qubits"] = qbs;
}

void from_json(const nlohmann::json& j, UnitaryTableau& tab) {
  j.at("tab").get_to(tab.tab_);
  if (tab.tab_.get_n_rows() != 2 * tab.tab_.get_n_qubits())
    throw std::invalid_argument(
        "Size of tableau does not match requirements for UnitaryTableau.");
  qubit_vector_t qbs = j.at("qubits").get<qubit_vector_t>();
  unsigned nqbs = qbs.size();
  if (nqbs != tab.tab_.get_n_qubits())
    throw std::invalid_argument(
        "Number of qubits in json UnitaryTableau does not match tableau size.");
  MatrixXb expected_anticommutes(2 * nqbs, 2 * nqbs);
  expected_anticommutes << MatrixXb::Zero(nqbs, nqbs),
      MatrixXb::Identity(nqbs, nqbs), MatrixXb::Identity(nqbs, nqbs),
      MatrixXb::Zero(nqbs, nqbs);
  if (tab.tab_.anticommuting_rows() != expected_anticommutes)
    throw std::invalid_argument(
        "Rows of tableau do not (anti-)commute as expected for UnitaryTableau");
  if (tab.tab_.rank() != 2 * nqbs)
    throw std::invalid_argument(
        "Rows of UnitaryTableau are not linearly independent");
  tab.qubits_.clear();
  for (unsigned i = 0; i < nqbs; ++i) {
    tab.qubits_.insert({qbs.at(i), i});
  }
}

UnitaryRevTableau::UnitaryRevTableau(unsigned n) : tab_(n) {}

UnitaryRevTableau::UnitaryRevTableau(const qubit_vector_t& qbs) : tab_(qbs) {}

UnitaryRevTableau::UnitaryRevTableau() : tab_() {}

SpPauliStabiliser UnitaryRevTableau::get_xrow(const Qubit& qb) const {
  return tab_.get_xrow(qb);
}

SpPauliStabiliser UnitaryRevTableau::get_zrow(const Qubit& qb) const {
  return tab_.get_zrow(qb);
}

SpPauliStabiliser UnitaryRevTableau::get_row_product(
    const SpPauliStabiliser& qpt) const {
  return tab_.get_row_product(qpt);
}

std::set<Qubit> UnitaryRevTableau::get_qubits() const {
  return tab_.get_qubits();
}

void UnitaryRevTableau::apply_S_at_front(const Qubit& qb) {
  tab_.apply_pauli_at_end(SpPauliStabiliser(qb, Pauli::Z), 3);
}

void UnitaryRevTableau::apply_S_at_end(const Qubit& qb) {
  tab_.apply_pauli_at_front(SpPauliStabiliser(qb, Pauli::Z), 3);
}

void UnitaryRevTableau::apply_Z_at_front(const Qubit& qb) {
  tab_.apply_Z_at_end(qb);
}

void UnitaryRevTableau::apply_Z_at_end(const Qubit& qb) {
  tab_.apply_Z_at_front(qb);
}

void UnitaryRevTableau::apply_V_at_front(const Qubit& qb) {
  tab_.apply_pauli_at_end(SpPauliStabiliser(qb, Pauli::X), 3);
}

void UnitaryRevTableau::apply_V_at_end(const Qubit& qb) {
  tab_.apply_pauli_at_front(SpPauliStabiliser(qb, Pauli::X), 3);
}

void UnitaryRevTableau::apply_X_at_front(const Qubit& qb) {
  tab_.apply_X_at_end(qb);
}

void UnitaryRevTableau::apply_X_at_end(const Qubit& qb) {
  tab_.apply_X_at_front(qb);
}

void UnitaryRevTableau::apply_H_at_front(const Qubit& qb) {
  tab_.apply_H_at_end(qb);
}

void UnitaryRevTableau::apply_H_at_end(const Qubit& qb) {
  tab_.apply_H_at_front(qb);
}

void UnitaryRevTableau::apply_CX_at_front(
    const Qubit& control, const Qubit& target) {
  tab_.apply_CX_at_end(control, target);
}

void UnitaryRevTableau::apply_CX_at_end(
    const Qubit& control, const Qubit& target) {
  tab_.apply_CX_at_front(control, target);
}

void UnitaryRevTableau::apply_gate_at_front(
    OpType type, const qubit_vector_t& qbs) {
  // Handle types whose dagger is not an optype
  switch (type) {
    case OpType::ZZMax: {
      tab_.apply_gate_at_end(OpType::ZZMax, qbs);
      tab_.apply_gate_at_end(OpType::Z, {qbs.at(0)});
      tab_.apply_gate_at_end(OpType::Z, {qbs.at(1)});
      break;
    }
    case OpType::ISWAPMax: {
      tab_.apply_gate_at_end(OpType::ISWAPMax, qbs);
      tab_.apply_gate_at_end(OpType::Z, {qbs.at(0)});
      tab_.apply_gate_at_end(OpType::Z, {qbs.at(1)});
      break;
    }
    case OpType::Phase: {
      break;
    }
    default: {
      tab_.apply_gate_at_end(get_op_ptr(type)->dagger()->get_type(), qbs);
      break;
    }
  }
}

void UnitaryRevTableau::apply_gate_at_end(
    OpType type, const qubit_vector_t& qbs) {
  // Handle types whose dagger is not an optype
  switch (type) {
    case OpType::ZZMax: {
      tab_.apply_gate_at_front(OpType::ZZMax, qbs);
      tab_.apply_gate_at_front(OpType::Z, {qbs.at(0)});
      tab_.apply_gate_at_front(OpType::Z, {qbs.at(1)});
      break;
    }
    case OpType::ISWAPMax: {
      tab_.apply_gate_at_front(OpType::ISWAPMax, qbs);
      tab_.apply_gate_at_front(OpType::Z, {qbs.at(0)});
      tab_.apply_gate_at_front(OpType::Z, {qbs.at(1)});
      break;
    }
    case OpType::Phase: {
      break;
    }
    default: {
      tab_.apply_gate_at_front(get_op_ptr(type)->dagger()->get_type(), qbs);
      break;
    }
  }
}

void UnitaryRevTableau::apply_pauli_at_front(
    const SpPauliStabiliser& pauli, unsigned half_pis) {
  tab_.apply_pauli_at_end(pauli, 4 - (half_pis % 4));
}

void UnitaryRevTableau::apply_pauli_at_end(
    const SpPauliStabiliser& pauli, unsigned half_pis) {
  tab_.apply_pauli_at_front(pauli, 4 - (half_pis % 4));
}

UnitaryRevTableau UnitaryRevTableau::compose(
    const UnitaryRevTableau& first, const UnitaryRevTableau& second) {
  UnitaryRevTableau result(0);
  result.tab_ = UnitaryTableau::compose(second.tab_, first.tab_);
  return result;
}

UnitaryRevTableau UnitaryRevTableau::dagger() const {
  UnitaryRevTableau result(0);
  result.tab_ = tab_.dagger();
  return result;
}

UnitaryRevTableau UnitaryRevTableau::transpose() const {
  UnitaryRevTableau result(0);
  result.tab_ = tab_.transpose();
  return result;
}

UnitaryRevTableau UnitaryRevTableau::conjugate() const {
  UnitaryRevTableau result(0);
  result.tab_ = tab_.conjugate();
  return result;
}

std::ostream& operator<<(std::ostream& os, const UnitaryRevTableau& tab) {
  unsigned nqs = tab.tab_.qubits_.size();
  for (unsigned i = 0; i < nqs; ++i) {
    Qubit qi = tab.tab_.qubits_.right.at(i);
    os << tab.tab_.tab_.xmat.row(i) << "   " << tab.tab_.tab_.zmat.row(i)
       << "   " << tab.tab_.tab_.phase(i) << "\t->\t"
       << "X@" << qi.repr() << std::endl;
  }
  os << "--" << std::endl;
  for (unsigned i = 0; i < nqs; ++i) {
    Qubit qi = tab.tab_.qubits_.right.at(i);
    os << tab.tab_.tab_.xmat.row(i + nqs) << "   "
       << tab.tab_.tab_.zmat.row(i + nqs) << "   "
       << tab.tab_.tab_.phase(i + nqs) << "\t->\t"
       << "Z@" << qi.repr() << std::endl;
  }
  return os;
}

bool UnitaryRevTableau::operator==(const UnitaryRevTableau& other) const {
  return (tab_ == other.tab_);
}

void to_json(nlohmann::json& j, const UnitaryRevTableau& tab) {
  j["tab"] = tab.tab_;
}

void from_json(const nlohmann::json& j, UnitaryRevTableau& tab) {
  j.at("tab").get_to(tab.tab_);
}

}  // namespace tket
