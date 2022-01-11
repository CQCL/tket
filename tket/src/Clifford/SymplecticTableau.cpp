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

#include "SymplecticTableau.hpp"

#include "OpType/OpTypeInfo.hpp"
#include "Utils/EigenConfig.hpp"
#include "Utils/Exceptions.hpp"

namespace tket {

bool BoolPauli::operator<(const BoolPauli &other) const {
  if (x == other.x) {
    return z < other.z;
  }
  return x < other.x;
}

Pauli BoolPauli::to_pauli() const {
  if (x) {
    if (z)
      return Pauli::Y;
    else
      return Pauli::X;
  } else if (z)
    return Pauli::Z;
  else
    return Pauli::I;
}

const std::map<std::pair<BoolPauli, BoolPauli>, std::pair<BoolPauli, Complex>>
    BoolPauli::mult_lut = {
        // {{{x1,z1},{x2,z2}}, {{x,z},ph}}
        {{{false, false}, {false, false}}, {{false, false}, 1.}},
        {{{false, false}, {false, true}}, {{false, true}, 1.}},
        {{{false, false}, {true, false}}, {{true, false}, 1.}},
        {{{false, false}, {true, true}}, {{true, true}, 1.}},
        {{{false, true}, {false, false}}, {{false, true}, 1.}},
        {{{false, true}, {false, true}}, {{false, false}, 1.}},
        {{{false, true}, {true, false}}, {{true, true}, i_}},
        {{{false, true}, {true, true}}, {{true, false}, -i_}},
        {{{true, false}, {false, false}}, {{true, false}, 1.}},
        {{{true, false}, {false, true}}, {{true, true}, -i_}},
        {{{true, false}, {true, false}}, {{false, false}, 1.}},
        {{{true, false}, {true, true}}, {{false, true}, i_}},
        {{{true, true}, {false, false}}, {{true, true}, 1.}},
        {{{true, true}, {false, true}}, {{true, false}, i_}},
        {{{true, true}, {true, false}}, {{false, true}, -i_}},
        {{{true, true}, {true, true}}, {{false, false}, 1.}}};

SymplecticTableau::SymplecticTableau(
    const MatrixXb &xmat, const MatrixXb &zmat, const VectorXb &phase)
    : n_rows_(xmat.rows()),
      n_qubits_(xmat.cols()),
      xmat_(xmat),
      zmat_(zmat),
      phase_(phase) {
  if (zmat.rows() != n_rows_ || phase_.size() != n_rows_)
    throw NotValid(
        "Tableau must have the same number of rows in each component.");
  if (zmat.cols() != n_qubits_)
    throw NotValid(
        "Tableau must have the same number of columns in x and z components.");
}

SymplecticTableau::SymplecticTableau(const PauliStabiliserList &rows) {
  n_rows_ = rows.size();
  if (n_rows_ == 0)
    n_qubits_ = 0;
  else
    n_qubits_ = rows[0].string.size();
  xmat_ = MatrixXb::Zero(n_rows_, n_qubits_);
  zmat_ = MatrixXb::Zero(n_rows_, n_qubits_);
  phase_ = VectorXb::Zero(n_rows_);
  for (unsigned i = 0; i < n_rows_; ++i) {
    const PauliStabiliser &stab = rows[i];
    if (stab.string.size() != n_qubits_)
      throw NotValid(
          "Tableau must have the same number of qubits in each row.");
    for (unsigned q = 0; q < n_qubits_; ++q) {
      const Pauli &p = stab.string[q];
      xmat_(i, q) = (p == Pauli::X) || (p == Pauli::Y);
      zmat_(i, q) = (p == Pauli::Z) || (p == Pauli::Y);
    }
    phase_(i) = !stab.coeff;
  }
}

unsigned SymplecticTableau::get_n_rows() const { return n_rows_; }
unsigned SymplecticTableau::get_n_qubits() const { return n_qubits_; }

PauliStabiliser SymplecticTableau::get_pauli(unsigned i) const {
  std::vector<Pauli> str(n_qubits_);
  for (unsigned q = 0; q < n_qubits_; ++q) {
    str[q] = BoolPauli{xmat_(i, q), zmat_(i, q)}.to_pauli();
  }
  return PauliStabiliser(str, !phase_(i));
}

std::ostream &operator<<(std::ostream &os, const SymplecticTableau &tab) {
  for (unsigned i = 0; i < tab.n_rows_; ++i) {
    os << tab.xmat_.row(i) << " " << tab.zmat_.row(i) << " " << tab.phase_(i)
       << std::endl;
  }
  return os;
}

bool SymplecticTableau::operator==(const SymplecticTableau &other) const {
  bool same = this->n_rows_ == other.n_rows_;
  same &= this->n_qubits_ == other.n_qubits_;
  same &= this->xmat_ == other.xmat_;
  same &= this->zmat_ == other.zmat_;
  same &= this->phase_ == other.phase_;
  return same;
}

void SymplecticTableau::row_mult(unsigned ra, unsigned rw, Complex coeff) {
  MatrixXb::RowXpr xa = xmat_.row(ra);
  MatrixXb::RowXpr za = zmat_.row(ra);
  MatrixXb::RowXpr xw = xmat_.row(rw);
  MatrixXb::RowXpr zw = zmat_.row(rw);
  row_mult(xa, za, phase_(ra), xw, zw, phase_(rw), coeff, xw, zw, phase_(rw));
}

void SymplecticTableau::apply_S(unsigned qb) {
  MatrixXb::ColXpr xcol = xmat_.col(qb);
  MatrixXb::ColXpr zcol = zmat_.col(qb);
  col_mult(xcol, zcol, true, zcol, phase_);
}

void SymplecticTableau::apply_V(unsigned qb) {
  MatrixXb::ColXpr xcol = xmat_.col(qb);
  MatrixXb::ColXpr zcol = zmat_.col(qb);
  col_mult(xcol, zcol, false, xcol, phase_);
}

void SymplecticTableau::apply_CX(unsigned qc, unsigned qt) {
  for (unsigned i = 0; i < n_rows_; ++i) {
    phase_(i) = phase_(i) ^ (xmat_(i, qc) && zmat_(i, qt) &&
                             !(xmat_(i, qt) ^ zmat_(i, qc)));
    xmat_(i, qt) = xmat_(i, qc) ^ xmat_(i, qt);
    zmat_(i, qc) = zmat_(i, qc) ^ zmat_(i, qt);
  }
}

void SymplecticTableau::apply_gate(
    OpType type, const std::vector<unsigned> &qbs) {
  switch (type) {
    case OpType::Z: {
      apply_S(qbs.at(0));
      apply_S(qbs.at(0));
      break;
    }
    case OpType::X: {
      apply_V(qbs.at(0));
      apply_V(qbs.at(0));
      break;
    }
    case OpType::Y: {
      apply_S(qbs.at(0));
      apply_S(qbs.at(0));
      apply_V(qbs.at(0));
      apply_V(qbs.at(0));
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
    case OpType::Vdg: {
      apply_V(qbs.at(0));
      apply_V(qbs.at(0));
      apply_V(qbs.at(0));
      break;
    }
    case OpType::H: {
      apply_S(qbs.at(0));
      apply_V(qbs.at(0));
      apply_S(qbs.at(0));
      break;
    }
    case OpType::CX: {
      apply_CX(qbs.at(0), qbs.at(1));
      break;
    }
    case OpType::CY: {
      apply_S(qbs.at(1));
      apply_S(qbs.at(1));
      apply_S(qbs.at(1));
      apply_CX(qbs.at(0), qbs.at(1));
      apply_S(qbs.at(1));
      break;
    }
    case OpType::CZ: {
      apply_S(qbs.at(1));
      apply_V(qbs.at(1));
      apply_S(qbs.at(1));
      apply_CX(qbs.at(0), qbs.at(1));
      apply_S(qbs.at(1));
      apply_V(qbs.at(1));
      apply_S(qbs.at(1));
      break;
    }
    case OpType::SWAP: {
      apply_CX(qbs.at(0), qbs.at(1));
      apply_CX(qbs.at(1), qbs.at(0));
      apply_CX(qbs.at(0), qbs.at(1));
      break;
    }
    case OpType::BRIDGE: {
      apply_CX(qbs.at(0), qbs.at(2));
      break;
    }
    case OpType::noop: {
      break;
    }
    default: {
      throw NotValid(
          optypeinfo().at(type).name +
          " cannot be applied to a SymplecticTableau; it is not a Clifford "
          "gate");
    }
  }
}

void SymplecticTableau::apply_pauli_gadget(
    const PauliStabiliser &pauli, unsigned half_pis) {
  if (pauli.string.size() != n_qubits_) {
    throw NotValid(
        "Cannot apply pauli gadget to SymplecticTableau; string and tableau "
        "have different numbers of qubits");
  }
  half_pis = half_pis % 4;
  if (half_pis == 0) return;  // Identity
  if (half_pis == 2) {        // Degenerates to product of PI rotations
    for (unsigned i = 0; i < pauli.string.size(); ++i) {
      switch (pauli.string.at(i)) {
        case Pauli::I: {
          break;
        }
        case Pauli::X: {
          apply_gate(OpType::X, {i});
          break;
        }
        case Pauli::Y: {
          apply_gate(OpType::Y, {i});
          break;
        }
        case Pauli::Z: {
          apply_gate(OpType::Z, {i});
          break;
        }
      }
    }
    return;
  }

  // From here, half_pis == 1 or 3
  // They act the same except for a phase flip on the product term
  MatrixXb pauli_xrow = MatrixXb::Zero(1, n_qubits_);
  MatrixXb pauli_zrow = MatrixXb::Zero(1, n_qubits_);
  for (unsigned i = 0; i < n_qubits_; ++i) {
    Pauli p = pauli.string.at(i);
    pauli_xrow(i) = (p == Pauli::X) || (p == Pauli::Y);
    pauli_zrow(i) = (p == Pauli::Z) || (p == Pauli::Y);
  }
  bool phase = (!pauli.coeff) ^ (half_pis == 3);

  for (unsigned i = 0; i < n_rows_; ++i) {
    bool anti = false;
    MatrixXb::RowXpr xr = xmat_.row(i);
    MatrixXb::RowXpr zr = zmat_.row(i);
    for (unsigned q = 0; q < n_qubits_; ++q) {
      anti ^= (xr(q) && pauli_zrow(q));
      anti ^= (zr(q) && pauli_xrow(q));
    }
    if (anti) {
      row_mult(
          xr, zr, phase_(i), pauli_xrow.row(0), pauli_zrow.row(0), phase, i_,
          xr, zr, phase_(i));
    }
  }
}

MatrixXb SymplecticTableau::anticommuting_rows() const {
  MatrixXb res = MatrixXb::Zero(n_rows_, n_rows_);
  for (unsigned i = 0; i < n_rows_; ++i) {
    for (unsigned j = 0; j < i; ++j) {
      bool anti = false;
      for (unsigned q = 0; q < n_qubits_; ++q) {
        anti ^= (xmat_(i, q) && zmat_(j, q));
        anti ^= (xmat_(j, q) && zmat_(i, q));
      }
      res(i, j) = anti;
      res(j, i) = anti;
    }
  }
  return res;
}

unsigned SymplecticTableau::rank() const {
  MatrixXb fullmat(n_rows_, 2 * n_qubits_);
  fullmat << xmat_, zmat_;
  Eigen::FullPivLU<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> lu(
      fullmat.cast<double>());
  return lu.rank();
}

SymplecticTableau SymplecticTableau::conjugate() const {
  SymplecticTableau conj(*this);
  for (unsigned i = 0; i < n_rows_; ++i) {
    unsigned sum = 0;
    for (unsigned j = 0; j < n_qubits_; ++j) {
      if (xmat_(i, j) && zmat_(i, j)) ++sum;
    }
    if (sum % 2 == 1) conj.phase_(i) ^= true;
  }
  return conj;
}

void SymplecticTableau::row_mult(
    const MatrixXb::RowXpr &xa, const MatrixXb::RowXpr &za, const bool &pa,
    const MatrixXb::RowXpr &xb, const MatrixXb::RowXpr &zb, const bool &pb,
    Complex phase, MatrixXb::RowXpr &xw, MatrixXb::RowXpr &zw, bool &pw) {
  if (pa) phase *= -1;
  if (pb) phase *= -1;
  for (unsigned i = 0; i < n_qubits_; i++) {
    std::pair<BoolPauli, Complex> res =
        BoolPauli::mult_lut.at({{xa(i), za(i)}, {xb(i), zb(i)}});
    xw(i) = res.first.x;
    zw(i) = res.first.z;
    phase *= res.second;
  }
  pw = (phase == -1.);
}

void SymplecticTableau::col_mult(
    const MatrixXb::ColXpr &a, const MatrixXb::ColXpr &b, bool flip,
    MatrixXb::ColXpr &w, VectorXb &pw) {
  for (unsigned i = 0; i < n_rows_; i++) {
    pw(i) = pw(i) ^ (a(i) && (b(i) ^ flip));
    w(i) = a(i) ^ b(i);
  }
}

void to_json(nlohmann::json &j, const SymplecticTableau &tab) {
  j["nrows"] = tab.n_rows_;
  j["nqubits"] = tab.n_qubits_;
  j["xmat"] = tab.xmat_;
  j["zmat"] = tab.zmat_;
  j["phase"] = tab.phase_;
}

void from_json(const nlohmann::json &j, SymplecticTableau &tab) {
  unsigned n_rows = j.at("nrows").get<unsigned>();
  unsigned n_qbs = j.at("nqubits").get<unsigned>();
  MatrixXb xmat(n_rows, n_qbs);
  MatrixXb zmat(n_rows, n_qbs);
  VectorXb phase(n_rows);
  from_json(j.at("xmat"), xmat);
  from_json(j.at("zmat"), zmat);
  from_json(j.at("phase"), phase);
  tab = SymplecticTableau(xmat, zmat, phase);
}

}  // namespace tket
