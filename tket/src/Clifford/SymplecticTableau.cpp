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

#include "tket/Clifford/SymplecticTableau.hpp"

#include <stdexcept>

#include "tket/OpType/OpTypeInfo.hpp"
#include "tket/Utils/EigenConfig.hpp"

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
    : xmat(xmat), zmat(zmat), phase(phase) {
  if (zmat.rows() != xmat.rows() || phase.size() != xmat.rows())
    throw std::invalid_argument(
        "Tableau must have the same number of rows in each component.");
  if (zmat.cols() != xmat.cols())
    throw std::invalid_argument(
        "Tableau must have the same number of columns in x and z components.");
}

SymplecticTableau::SymplecticTableau(const PauliStabiliserVec &rows) {
  unsigned n_rows = rows.size();
  unsigned n_qubits = 0;
  if (n_rows != 0) n_qubits = rows[0].string.size();
  xmat = MatrixXb::Zero(n_rows, n_qubits);
  zmat = MatrixXb::Zero(n_rows, n_qubits);
  phase = VectorXb::Zero(n_rows);
  for (unsigned i = 0; i < n_rows; ++i) {
    const PauliStabiliser &stab = rows[i];
    if (stab.string.size() != n_qubits)
      throw std::invalid_argument(
          "Tableau must have the same number of qubits in each row.");
    for (unsigned q = 0; q < n_qubits; ++q) {
      const Pauli &p = stab.get(q);
      xmat(i, q) = (p == Pauli::X) || (p == Pauli::Y);
      zmat(i, q) = (p == Pauli::Z) || (p == Pauli::Y);
    }
    phase(i) = stab.is_real_negative();
  }
}

unsigned SymplecticTableau::get_n_rows() const { return xmat.rows(); }
unsigned SymplecticTableau::get_n_qubits() const { return xmat.cols(); }

PauliStabiliser SymplecticTableau::get_pauli(unsigned i) const {
  unsigned n_qubits = get_n_qubits();
  std::vector<Pauli> str(n_qubits);
  for (unsigned q = 0; q < n_qubits; ++q) {
    str[q] = BoolPauli{xmat(i, q), zmat(i, q)}.to_pauli();
  }
  return PauliStabiliser(str, phase(i) ? 2 : 0);
}

std::ostream &operator<<(std::ostream &os, const SymplecticTableau &tab) {
  for (unsigned i = 0; i < tab.get_n_rows(); ++i) {
    os << tab.xmat.row(i) << " " << tab.zmat.row(i) << " " << tab.phase(i)
       << std::endl;
  }
  return os;
}

bool SymplecticTableau::operator==(const SymplecticTableau &other) const {
  // Need this to short-circuit before matrix checks as comparing matrices of
  // different sizes will throw an exception
  return (this->get_n_rows() == other.get_n_rows()) &&
         (this->get_n_qubits() == other.get_n_qubits()) &&
         (this->xmat == other.xmat) && (this->zmat == other.zmat) &&
         (this->phase == other.phase);
}

void SymplecticTableau::row_mult(unsigned ra, unsigned rw, Complex coeff) {
  MatrixXb::RowXpr xa = xmat.row(ra);
  MatrixXb::RowXpr za = zmat.row(ra);
  MatrixXb::RowXpr xw = xmat.row(rw);
  MatrixXb::RowXpr zw = zmat.row(rw);
  row_mult(xa, za, phase(ra), xw, zw, phase(rw), coeff, xw, zw, phase(rw));
}

void SymplecticTableau::apply_S(unsigned qb) {
  MatrixXb::ColXpr xcol = xmat.col(qb);
  MatrixXb::ColXpr zcol = zmat.col(qb);
  col_mult(xcol, zcol, false, zcol, phase);
}

void SymplecticTableau::apply_Z(unsigned qb) {
  for (unsigned i = 0; i < get_n_rows(); ++i) phase(i) = phase(i) ^ xmat(i, qb);
}

void SymplecticTableau::apply_V(unsigned qb) {
  MatrixXb::ColXpr xcol = xmat.col(qb);
  MatrixXb::ColXpr zcol = zmat.col(qb);
  col_mult(zcol, xcol, true, xcol, phase);
}

void SymplecticTableau::apply_X(unsigned qb) {
  for (unsigned i = 0; i < get_n_rows(); ++i) phase(i) = phase(i) ^ zmat(i, qb);
}

void SymplecticTableau::apply_H(unsigned qb) {
  for (unsigned i = 0; i < get_n_rows(); ++i) {
    phase(i) = phase(i) ^ (xmat(i, qb) && zmat(i, qb));
    bool temp = xmat(i, qb);
    xmat(i, qb) = zmat(i, qb);
    zmat(i, qb) = temp;
  }
}

void SymplecticTableau::apply_CX(unsigned qc, unsigned qt) {
  if (qc == qt)
    throw std::logic_error(
        "Attempting to apply a CX with equal control and target in a tableau");
  for (unsigned i = 0; i < get_n_rows(); ++i) {
    phase(i) =
        phase(i) ^ (xmat(i, qc) && zmat(i, qt) && !(xmat(i, qt) ^ zmat(i, qc)));
    xmat(i, qt) = xmat(i, qc) ^ xmat(i, qt);
    zmat(i, qc) = zmat(i, qc) ^ zmat(i, qt);
  }
}

void SymplecticTableau::apply_gate(
    OpType type, const std::vector<unsigned> &qbs) {
  switch (type) {
    case OpType::Z: {
      apply_Z(qbs.at(0));
      break;
    }
    case OpType::X: {
      apply_X(qbs.at(0));
      break;
    }
    case OpType::Y: {
      apply_Z(qbs.at(0));
      apply_X(qbs.at(0));
      break;
    }
    case OpType::S: {
      apply_S(qbs.at(0));
      break;
    }
    case OpType::Sdg: {
      apply_S(qbs.at(0));
      apply_Z(qbs.at(0));
      break;
    }
    case OpType::V:
    case OpType::SX: {
      apply_V(qbs.at(0));
      break;
    }
    case OpType::Vdg:
    case OpType::SXdg: {
      apply_V(qbs.at(0));
      apply_X(qbs.at(0));
      break;
    }
    case OpType::H: {
      apply_H(qbs.at(0));
      break;
    }
    case OpType::CX: {
      apply_CX(qbs.at(0), qbs.at(1));
      break;
    }
    case OpType::CY: {
      apply_S(qbs.at(1));
      apply_Z(qbs.at(1));
      apply_CX(qbs.at(0), qbs.at(1));
      apply_S(qbs.at(1));
      break;
    }
    case OpType::CZ: {
      apply_H(qbs.at(1));
      apply_CX(qbs.at(0), qbs.at(1));
      apply_H(qbs.at(1));
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
    case OpType::ZZMax: {
      apply_H(qbs.at(1));
      apply_S(qbs.at(0));
      apply_V(qbs.at(1));
      apply_CX(qbs.at(0), qbs.at(1));
      apply_H(qbs.at(1));
      break;
    }
    case OpType::ECR: {
      apply_S(qbs.at(0));
      apply_X(qbs.at(0));
      apply_V(qbs.at(1));
      apply_X(qbs.at(1));
      apply_CX(qbs.at(0), qbs.at(1));
      break;
    }
    case OpType::ISWAPMax: {
      apply_V(qbs.at(0));
      apply_V(qbs.at(1));
      apply_CX(qbs.at(0), qbs.at(1));
      apply_V(qbs.at(0));
      apply_S(qbs.at(1));
      apply_Z(qbs.at(1));
      apply_CX(qbs.at(0), qbs.at(1));
      apply_V(qbs.at(0));
      apply_V(qbs.at(1));
      break;
    }
    case OpType::noop:
    case OpType::Phase: {
      break;
    }
    default: {
      throw BadOpType(
          "Cannot be applied to a SymplecticTableau: not a Clifford gate",
          type);
    }
  }
}

void SymplecticTableau::apply_pauli_gadget(
    const PauliStabiliser &pauli, unsigned half_pis) {
  unsigned n_qubits = get_n_qubits();
  if (pauli.string.size() != n_qubits) {
    throw std::invalid_argument(
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
  MatrixXb pauli_xrow = MatrixXb::Zero(1, n_qubits);
  MatrixXb pauli_zrow = MatrixXb::Zero(1, n_qubits);
  for (unsigned i = 0; i < n_qubits; ++i) {
    Pauli p = pauli.string.at(i);
    pauli_xrow(i) = (p == Pauli::X) || (p == Pauli::Y);
    pauli_zrow(i) = (p == Pauli::Z) || (p == Pauli::Y);
  }
  bool phase_flip = pauli.is_real_negative() ^ (half_pis == 3);

  for (unsigned i = 0; i < get_n_rows(); ++i) {
    bool anti = false;
    MatrixXb::RowXpr xr = xmat.row(i);
    MatrixXb::RowXpr zr = zmat.row(i);
    for (unsigned q = 0; q < n_qubits; ++q) {
      anti ^= (xr(q) && pauli_zrow(q));
      anti ^= (zr(q) && pauli_xrow(q));
    }
    if (anti) {
      row_mult(
          xr, zr, phase(i), pauli_xrow.row(0), pauli_zrow.row(0), phase_flip,
          i_, xr, zr, phase(i));
    }
  }
}

MatrixXb SymplecticTableau::anticommuting_rows() const {
  unsigned n_rows = get_n_rows();
  MatrixXb res = MatrixXb::Zero(n_rows, n_rows);
  for (unsigned i = 0; i < n_rows; ++i) {
    for (unsigned j = 0; j < i; ++j) {
      bool anti = false;
      for (unsigned q = 0; q < get_n_qubits(); ++q) {
        anti ^= (xmat(i, q) && zmat(j, q));
        anti ^= (xmat(j, q) && zmat(i, q));
      }
      res(i, j) = anti;
      res(j, i) = anti;
    }
  }
  return res;
}

unsigned SymplecticTableau::rank() const {
  // Create a copy in gaussian form and count the empty rows
  SymplecticTableau copy(*this);
  copy.gaussian_form();
  unsigned empty_rows = 0;
  unsigned n_rows = get_n_rows();
  for (unsigned i = 0; i < n_rows; ++i) {
    if (copy.xmat.row(n_rows - 1 - i).isZero() &&
        copy.zmat.row(n_rows - 1 - i).isZero())
      ++empty_rows;
    else
      break;
  }
  return n_rows - empty_rows;
}

SymplecticTableau SymplecticTableau::conjugate() const {
  SymplecticTableau conj(*this);
  for (unsigned i = 0; i < get_n_rows(); ++i) {
    unsigned sum = 0;
    for (unsigned j = 0; j < get_n_qubits(); ++j) {
      if (xmat(i, j) && zmat(i, j)) ++sum;
    }
    if (sum % 2 == 1) conj.phase(i) ^= true;
  }
  return conj;
}

void SymplecticTableau::gaussian_form() {
  MatrixXb fullmat = MatrixXb::Zero(get_n_rows(), 2 * get_n_qubits());
  fullmat(Eigen::all, Eigen::seq(0, Eigen::last, 2)) = xmat;
  fullmat(Eigen::all, Eigen::seq(1, Eigen::last, 2)) = zmat;
  std::vector<std::pair<unsigned, unsigned>> row_ops =
      gaussian_elimination_row_ops(fullmat);
  for (const std::pair<unsigned, unsigned> &op : row_ops) {
    row_mult(op.first, op.second);
  }
}

void SymplecticTableau::row_mult(
    const MatrixXb::RowXpr &xa, const MatrixXb::RowXpr &za, const bool &pa,
    const MatrixXb::RowXpr &xb, const MatrixXb::RowXpr &zb, const bool &pb,
    Complex phase, MatrixXb::RowXpr &xw, MatrixXb::RowXpr &zw, bool &pw) {
  if (pa) phase *= -1;
  if (pb) phase *= -1;
  for (unsigned i = 0; i < get_n_qubits(); i++) {
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
  for (unsigned i = 0; i < get_n_rows(); i++) {
    pw(i) = pw(i) ^ (a(i) && (b(i) ^ flip));
    w(i) = a(i) ^ b(i);
  }
}

void to_json(nlohmann::json &j, const SymplecticTableau &tab) {
  j["nrows"] = tab.get_n_rows();
  j["nqubits"] = tab.get_n_qubits();
  j["xmat"] = tab.xmat;
  j["zmat"] = tab.zmat;
  j["phase"] = tab.phase;
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
