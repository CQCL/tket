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

#include "CliffTableau.hpp"

#include "OpType/OpType.hpp"
#include "OpType/OpTypeInfo.hpp"
#include "Utils/Exceptions.hpp"
#include "Utils/MatrixAnalysis.hpp"

namespace tket {

const std::map<
    std::pair<CliffTableau::BoolPauli, CliffTableau::BoolPauli>,
    std::pair<CliffTableau::BoolPauli, Complex>>
    CliffTableau::mult_lut = {
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

CliffTableau::CliffTableau(unsigned n) : size_(n) {
  xpauli_x = MatrixXb::Identity(n, n);
  xpauli_z = MatrixXb::Zero(n, n);
  xpauli_phase = VectorXb::Zero(n);
  zpauli_x = MatrixXb::Zero(n, n);
  zpauli_z = MatrixXb::Identity(n, n);
  zpauli_phase = VectorXb::Zero(n);
  for (unsigned i = 0; i < n; i++) {
    qubits_.insert({Qubit(q_default_reg(), i), i});
  }
}

CliffTableau::CliffTableau(const qubit_vector_t &qbs) : size_(qbs.size()) {
  const unsigned n = size_;
  xpauli_x = MatrixXb::Identity(n, n);
  xpauli_z = MatrixXb::Zero(n, n);
  xpauli_phase = VectorXb::Zero(n);
  zpauli_x = MatrixXb::Zero(n, n);
  zpauli_z = MatrixXb::Identity(n, n);
  zpauli_phase = VectorXb::Zero(n);

  unsigned i = 0;
  for (const Qubit &q : qbs) {
    qubits_.insert({q, i});
    i++;
  }
}

static QubitPauliTensor get_pauli(
    const Qubit &qb, const boost::bimap<Qubit, unsigned> &qubits_,
    const MatrixXb &pauli_x, const MatrixXb &pauli_z,
    const VectorXb &pauli_phase) {
  unsigned uqb = qubits_.left.at(qb);
  Complex phase = 1.;
  if (pauli_phase(uqb)) phase = -1.;
  QubitPauliTensor res(phase);
  for (boost::bimap<Qubit, unsigned>::const_iterator iter = qubits_.begin(),
                                                     iend = qubits_.end();
       iter != iend; ++iter) {
    unsigned origin = iter->right;
    if (pauli_x(uqb, origin)) {
      if (pauli_z(uqb, origin)) {
        res = res * QubitPauliTensor(iter->left, Pauli::Y);
      } else {
        res = res * QubitPauliTensor(iter->left, Pauli::X);
      }
    } else if (pauli_z(uqb, origin)) {
      res = res * QubitPauliTensor(iter->left, Pauli::Z);
    }
  }
  return res;
}

QubitPauliTensor CliffTableau::get_zpauli(const Qubit &qb) const {
  return get_pauli(qb, qubits_, zpauli_x, zpauli_z, zpauli_phase);
}

QubitPauliTensor CliffTableau::get_xpauli(const Qubit &qb) const {
  return get_pauli(qb, qubits_, xpauli_x, xpauli_z, xpauli_phase);
}

std::set<Qubit> CliffTableau::get_qubits() const {
  std::set<Qubit> result;
  for (boost::bimap<Qubit, unsigned>::const_iterator iter = qubits_.begin(),
                                                     iend = qubits_.end();
       iter != iend; ++iter) {
    result.insert(iter->left);
  }
  return result;
}

void CliffTableau::row_mult(
    const MatrixXb::RowXpr &xa, const MatrixXb::RowXpr &za, const bool &pa,
    const MatrixXb::RowXpr &xb, const MatrixXb::RowXpr &zb, const bool &pb,
    Complex phase, MatrixXb::RowXpr &xw, MatrixXb::RowXpr &zw, bool &pw) {
  if (pa) phase *= -1;
  if (pb) phase *= -1;
  for (unsigned i = 0; i < size_; i++) {
    std::pair<BoolPauli, Complex> res =
        mult_lut.at({{xa(i), za(i)}, {xb(i), zb(i)}});
    xw(i) = res.first.x;
    zw(i) = res.first.z;
    phase *= res.second;
  }
  pw = (phase == -1.);
}

void CliffTableau::col_mult(
    const MatrixXb::ColXpr &a, const MatrixXb::ColXpr &b, bool flip,
    MatrixXb::ColXpr &w, VectorXb &pw) {
  for (unsigned i = 0; i < size_; i++) {
    pw(i) = pw(i) ^ (a(i) && (b(i) ^ flip));
    w(i) = a(i) ^ b(i);
  }
}

void CliffTableau::apply_S_at_front(unsigned qb) {
  MatrixXb::ColXpr xx = xpauli_x.col(qb);
  MatrixXb::ColXpr xz = xpauli_z.col(qb);
  MatrixXb::ColXpr zx = zpauli_x.col(qb);
  MatrixXb::ColXpr zz = zpauli_z.col(qb);
  col_mult(xx, xz, true, xz, xpauli_phase);
  col_mult(zx, zz, true, zz, zpauli_phase);
}

void CliffTableau::apply_S_at_end(unsigned qb) {
  MatrixXb::RowXpr xx = xpauli_x.row(qb);
  MatrixXb::RowXpr xz = xpauli_z.row(qb);
  MatrixXb::RowXpr zx = zpauli_x.row(qb);
  MatrixXb::RowXpr zz = zpauli_z.row(qb);
  row_mult(
      zx, zz, zpauli_phase(qb), xx, xz, xpauli_phase(qb), i_, xx, xz,
      xpauli_phase(qb));
}

void CliffTableau::apply_V_at_front(unsigned qb) {
  MatrixXb::ColXpr xx = xpauli_x.col(qb);
  MatrixXb::ColXpr xz = xpauli_z.col(qb);
  MatrixXb::ColXpr zx = zpauli_x.col(qb);
  MatrixXb::ColXpr zz = zpauli_z.col(qb);
  col_mult(xx, xz, false, xx, xpauli_phase);
  col_mult(zx, zz, false, zx, zpauli_phase);
}

void CliffTableau::apply_V_at_end(unsigned qb) {
  MatrixXb::RowXpr xx = xpauli_x.row(qb);
  MatrixXb::RowXpr xz = xpauli_z.row(qb);
  MatrixXb::RowXpr zx = zpauli_x.row(qb);
  MatrixXb::RowXpr zz = zpauli_z.row(qb);
  row_mult(
      xx, xz, xpauli_phase(qb), zx, zz, zpauli_phase(qb), i_, zx, zz,
      zpauli_phase(qb));
}

void CliffTableau::apply_CX_at_front(unsigned control, unsigned target) {
  for (unsigned i = 0; i < size_; i++) {
    xpauli_phase(i) =
        xpauli_phase(i) ^ (xpauli_x(i, control) && xpauli_z(i, target) &&
                           !(xpauli_x(i, target) ^ xpauli_z(i, control)));
    xpauli_x(i, target) = xpauli_x(i, control) ^ xpauli_x(i, target);
    xpauli_z(i, control) = xpauli_z(i, control) ^ xpauli_z(i, target);
    zpauli_phase(i) =
        zpauli_phase(i) ^ (zpauli_x(i, control) && zpauli_z(i, target) &&
                           !(zpauli_x(i, target) ^ zpauli_z(i, control)));
    zpauli_x(i, target) = zpauli_x(i, control) ^ zpauli_x(i, target);
    zpauli_z(i, control) = zpauli_z(i, control) ^ zpauli_z(i, target);
  }
}

void CliffTableau::apply_CX_at_end(unsigned control, unsigned target) {
  MatrixXb::RowXpr xxc = xpauli_x.row(control);
  MatrixXb::RowXpr xzc = xpauli_z.row(control);
  MatrixXb::RowXpr zxt = zpauli_x.row(target);
  MatrixXb::RowXpr zzt = zpauli_z.row(target);
  row_mult(
      xxc, xzc, xpauli_phase(control), xpauli_x.row(target),
      xpauli_z.row(target), xpauli_phase(target), 1., xxc, xzc,
      xpauli_phase(control));
  row_mult(
      zpauli_x.row(control), zpauli_z.row(control), zpauli_phase(control), zxt,
      zzt, zpauli_phase(target), 1., zxt, zzt, zpauli_phase(target));
}

void CliffTableau::apply_gate_at_front(
    OpType type, const std::vector<unsigned> &qbs) {
  switch (type) {
    case OpType::Z: {
      apply_S_at_front(qbs.at(0));
      apply_S_at_front(qbs.at(0));
      break;
    }
    case OpType::X: {
      apply_V_at_front(qbs.at(0));
      apply_V_at_front(qbs.at(0));
      break;
    }
    case OpType::Y: {
      apply_S_at_front(qbs.at(0));
      apply_S_at_front(qbs.at(0));
      apply_V_at_front(qbs.at(0));
      apply_V_at_front(qbs.at(0));
      break;
    }
    case OpType::S: {
      apply_S_at_front(qbs.at(0));
      break;
    }
    case OpType::Sdg: {
      apply_S_at_front(qbs.at(0));
      apply_S_at_front(qbs.at(0));
      apply_S_at_front(qbs.at(0));
      break;
    }
    case OpType::V: {
      apply_V_at_front(qbs.at(0));
      break;
    }
    case OpType::Vdg: {
      apply_V_at_front(qbs.at(0));
      apply_V_at_front(qbs.at(0));
      apply_V_at_front(qbs.at(0));
      break;
    }
    case OpType::H: {
      apply_S_at_front(qbs.at(0));
      apply_V_at_front(qbs.at(0));
      apply_S_at_front(qbs.at(0));
      break;
    }
    case OpType::CX: {
      apply_CX_at_front(qbs.at(0), qbs.at(1));
      break;
    }
    case OpType::CY: {
      apply_V_at_front(qbs.at(1));
      apply_V_at_front(qbs.at(1));
      apply_V_at_front(qbs.at(1));
      apply_CX_at_front(qbs.at(0), qbs.at(1));
      apply_V_at_front(qbs.at(1));
      break;
    }
    case OpType::CZ: {
      apply_S_at_front(qbs.at(1));
      apply_V_at_front(qbs.at(1));
      apply_S_at_front(qbs.at(1));
      apply_CX_at_front(qbs.at(0), qbs.at(1));
      apply_S_at_front(qbs.at(1));
      apply_V_at_front(qbs.at(1));
      apply_S_at_front(qbs.at(1));
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
    case OpType::noop: {
      break;
    }
    default: {
      throw NotValid(optypeinfo().at(type).name + " is not a Clifford gate");
    }
  }
}

void CliffTableau::apply_gate_at_front(OpType type, const qubit_vector_t &qbs) {
  std::vector<unsigned> uqbs;
  for (const Qubit &qb : qbs) {
    uqbs.push_back(qubits_.left.at(qb));
  }
  apply_gate_at_front(type, uqbs);
}

void CliffTableau::apply_gate_at_end(
    OpType type, const std::vector<unsigned> &qbs) {
  switch (type) {
    case OpType::Z: {
      apply_S_at_end(qbs.at(0));
      apply_S_at_end(qbs.at(0));
      break;
    }
    case OpType::X: {
      apply_V_at_end(qbs.at(0));
      apply_V_at_end(qbs.at(0));
      break;
    }
    case OpType::Y: {
      apply_S_at_end(qbs.at(0));
      apply_S_at_end(qbs.at(0));
      apply_V_at_end(qbs.at(0));
      apply_V_at_end(qbs.at(0));
      break;
    }
    case OpType::S: {
      apply_S_at_end(qbs.at(0));
      break;
    }
    case OpType::Sdg: {
      apply_S_at_end(qbs.at(0));
      apply_S_at_end(qbs.at(0));
      apply_S_at_end(qbs.at(0));
      break;
    }
    case OpType::V: {
      apply_V_at_end(qbs.at(0));
      break;
    }
    case OpType::Vdg: {
      apply_V_at_end(qbs.at(0));
      apply_V_at_end(qbs.at(0));
      apply_V_at_end(qbs.at(0));
      break;
    }
    case OpType::H: {
      apply_S_at_end(qbs.at(0));
      apply_V_at_end(qbs.at(0));
      apply_S_at_end(qbs.at(0));
      break;
    }
    case OpType::CX: {
      apply_CX_at_end(qbs.at(0), qbs.at(1));
      break;
    }
    case OpType::CY: {
      apply_V_at_end(qbs.at(1));
      apply_V_at_end(qbs.at(1));
      apply_V_at_end(qbs.at(1));
      apply_CX_at_end(qbs.at(0), qbs.at(1));
      apply_V_at_end(qbs.at(1));
      break;
    }
    case OpType::CZ: {
      apply_S_at_end(qbs.at(1));
      apply_V_at_end(qbs.at(1));
      apply_S_at_end(qbs.at(1));
      apply_CX_at_end(qbs.at(0), qbs.at(1));
      apply_S_at_end(qbs.at(1));
      apply_V_at_end(qbs.at(1));
      apply_S_at_end(qbs.at(1));
      break;
    }
    case OpType::SWAP: {
      apply_CX_at_end(qbs.at(0), qbs.at(1));
      apply_CX_at_end(qbs.at(1), qbs.at(0));
      apply_CX_at_end(qbs.at(0), qbs.at(1));
      break;
    }
    case OpType::BRIDGE: {
      apply_CX_at_end(qbs.at(0), qbs.at(2));
      break;
    }
    case OpType::noop: {
      break;
    }
    default: {
      throw NotValid(optypeinfo().at(type).name + " is not a Clifford gate");
    }
  }
}

void CliffTableau::apply_gate_at_end(OpType type, const qubit_vector_t &qbs) {
  std::vector<unsigned> uqbs;
  for (const Qubit &qb : qbs) {
    uqbs.push_back(qubits_.left.at(qb));
  }
  apply_gate_at_end(type, uqbs);
}

void CliffTableau::apply_pauli_at_front(
    const QubitPauliTensor &pauli, unsigned half_pis) {
  CliffTableau pauli_tab(size_);
  pauli_tab.qubits_ = this->qubits_;
  pauli_tab.apply_pauli_at_end(pauli, half_pis);
  CliffTableau product = CliffTableau::compose(pauli_tab, *this);
  this->xpauli_x = product.xpauli_x;
  this->xpauli_z = product.xpauli_z;
  this->xpauli_phase = product.xpauli_phase;
  this->zpauli_x = product.zpauli_x;
  this->zpauli_z = product.zpauli_z;
  this->zpauli_phase = product.zpauli_phase;
}

void CliffTableau::apply_pauli_at_end(
    const QubitPauliTensor &pauli, unsigned half_pis) {
  half_pis = half_pis % 4;
  if (half_pis == 0) return;  // Identity
  if (half_pis == 2) {        // Degenerates to product of PI rotations
    for (const std::pair<const Qubit, Pauli> &term : pauli.string.map) {
      switch (term.second) {
        case Pauli::I: {
          break;
        }
        case Pauli::X: {
          apply_gate_at_end(OpType::X, {term.first});
          break;
        }
        case Pauli::Y: {
          apply_gate_at_end(OpType::Y, {term.first});
          break;
        }
        case Pauli::Z: {
          apply_gate_at_end(OpType::Z, {term.first});
          break;
        }
      }
    }
    return;
  }

  // From here, half_pis == 1 or 3
  // They act the same except for a phase flip on the product term
  MatrixXb product_x = MatrixXb::Zero(1, size_);
  MatrixXb product_z = MatrixXb::Zero(1, size_);
  MatrixXb::RowXpr px = product_x.row(0);
  MatrixXb::RowXpr pz = product_z.row(0);
  if (pauli.coeff != 1. && pauli.coeff != -1.)
    throw NotValid(
        "Can only apply Paulis with real unit coefficients to "
        "CliffTableaus");
  bool phase = (pauli.coeff == -1.) ^ (half_pis == 3);

  // Collect the product term
  for (const std::pair<const Qubit, Pauli> &term : pauli.string.map) {
    unsigned uqb = qubits_.left.at(term.first);
    switch (term.second) {
      case Pauli::I: {
        break;
      }
      case Pauli::X: {
        row_mult(
            xpauli_x.row(uqb), xpauli_z.row(uqb), xpauli_phase(uqb), px, pz,
            phase, 1., px, pz, phase);
        break;
      }
      case Pauli::Y: {
        row_mult(
            zpauli_x.row(uqb), zpauli_z.row(uqb), zpauli_phase(uqb), px, pz,
            phase, 1., px, pz, phase);
        row_mult(
            xpauli_x.row(uqb), xpauli_z.row(uqb), xpauli_phase(uqb), px, pz,
            phase, i_, px, pz, phase);
        break;
      }
      case Pauli::Z: {
        row_mult(
            zpauli_x.row(uqb), zpauli_z.row(uqb), zpauli_phase(uqb), px, pz,
            phase, 1., px, pz, phase);
        break;
      }
    }
  }

  // Apply the product term on the anti-commuting rows
  for (const std::pair<const Qubit, Pauli> &term : pauli.string.map) {
    unsigned uqb = qubits_.left.at(term.first);
    MatrixXb::RowXpr xx = xpauli_x.row(uqb);
    MatrixXb::RowXpr xz = xpauli_z.row(uqb);
    MatrixXb::RowXpr zx = zpauli_x.row(uqb);
    MatrixXb::RowXpr zz = zpauli_z.row(uqb);
    switch (term.second) {
      case Pauli::I: {
        break;
      }
      case Pauli::X: {
        row_mult(
            px, pz, phase, zx, zz, zpauli_phase(uqb), i_, zx, zz,
            zpauli_phase(uqb));
        break;
      }
      case Pauli::Y: {
        row_mult(
            px, pz, phase, zx, zz, zpauli_phase(uqb), i_, zx, zz,
            zpauli_phase(uqb));
        row_mult(
            px, pz, phase, xx, xz, xpauli_phase(uqb), i_, xx, xz,
            xpauli_phase(uqb));
        break;
      }
      case Pauli::Z: {
        row_mult(
            px, pz, phase, xx, xz, xpauli_phase(uqb), i_, xx, xz,
            xpauli_phase(uqb));
        break;
      }
    }
  }
}

std::ostream &operator<<(std::ostream &os, CliffTableau const &tab) {
  MatrixXb full(2 * tab.size_, 2 * tab.size_ + 1);
  full << tab.xpauli_x, tab.xpauli_z, tab.xpauli_phase, tab.zpauli_x,
      tab.zpauli_z, tab.zpauli_phase;
  os << full;
  return os;
}

bool CliffTableau::operator==(const CliffTableau &other) const {
  bool same = this->size_ == other.size_;
  same &= this->qubits_ == other.qubits_;
  same &= this->xpauli_x == other.xpauli_x;
  same &= this->xpauli_z == other.xpauli_z;
  same &= this->xpauli_phase == other.xpauli_phase;
  same &= this->zpauli_x == other.zpauli_x;
  same &= this->zpauli_z == other.zpauli_z;
  same &= this->zpauli_phase == other.zpauli_phase;
  return same;
}

static std::vector<Complex> cphase_from_tableau(
    const MatrixXb &pauli_x, const MatrixXb &pauli_z,
    const VectorXb &pauli_phase) {
  std::vector<Complex> cphase;
  const unsigned n = pauli_phase.size();
  for (unsigned i = 0; i < n; i++) {
    Complex cp = pauli_phase(i) ? -1. : 1.;
    for (unsigned j = 0; j < n; j++) {
      if (pauli_x(i, j) && pauli_z(i, j)) cp *= i_;
    }
    cphase.push_back(cp);
  }
  return cphase;
}

CliffTableau CliffTableau::compose(
    const CliffTableau &first, const CliffTableau &second) {
  if (first.qubits_ != second.qubits_)
    throw NotImplemented(
        "Cannot compose Clifford Tableaus with different qubit maps");
  const unsigned n = first.size_;
  std::vector<Complex> first_x_cphase =
      cphase_from_tableau(first.xpauli_x, first.xpauli_z, first.xpauli_phase);
  std::vector<Complex> first_z_cphase =
      cphase_from_tableau(first.zpauli_x, first.zpauli_z, first.zpauli_phase);
  std::vector<Complex> second_x_cphase = cphase_from_tableau(
      second.xpauli_x, second.xpauli_z, second.xpauli_phase);
  std::vector<Complex> second_z_cphase = cphase_from_tableau(
      second.zpauli_x, second.zpauli_z, second.zpauli_phase);

  CliffTableau result(n);
  result.xpauli_x.setZero();
  result.xpauli_z.setZero();
  result.zpauli_x.setZero();
  result.zpauli_z.setZero();
  result.qubits_ = first.qubits_;

  // Target cphases
  std::vector<Complex> xpauli_cphase(n);
  std::vector<Complex> zpauli_cphase(n);

  // Calculate xpauli rows
  for (unsigned i = 0; i < n; i++) {  // row in second
    // Sum rows from first according to second.xpauli
    xpauli_cphase[i] = second_x_cphase[i];
    bool phase_flip = false;
    for (unsigned j = 0; j < n; j++) {  // col in second
      if (second.xpauli_x(i, j)) {
        // add first.xpauli to xpauli
        for (unsigned k = 0; k < n; k++) {
          phase_flip ^= result.xpauli_z(i, k) && first.xpauli_x(j, k);
          result.xpauli_x(i, k) = result.xpauli_x(i, k) != first.xpauli_x(j, k);
          result.xpauli_z(i, k) = result.xpauli_z(i, k) != first.xpauli_z(j, k);
        }
        xpauli_cphase[i] *= first_x_cphase[j];
      }
      if (second.xpauli_z(i, j)) {
        // add first.zpauli to xpauli
        for (unsigned k = 0; k < n; k++) {
          phase_flip ^= result.xpauli_z(i, k) && first.zpauli_x(j, k);
          result.xpauli_x(i, k) = result.xpauli_x(i, k) != first.zpauli_x(j, k);
          result.xpauli_z(i, k) = result.xpauli_z(i, k) != first.zpauli_z(j, k);
        }
        xpauli_cphase[i] *= first_z_cphase[j];
      }
    }
    if (phase_flip) xpauli_cphase[i] *= -1.;
  }

  // Calculate zpauli rows
  for (unsigned i = 0; i < n; i++) {  // row in second
    // Sum rows from first according to second.zpauli
    zpauli_cphase[i] = second_z_cphase[i];
    bool phase_flip = false;
    for (unsigned j = 0; j < n; j++) {  // col in second
      if (second.zpauli_x(i, j)) {
        // add first.xpauli to zpauli
        for (unsigned k = 0; k < n; k++) {
          phase_flip ^= result.zpauli_z(i, k) && first.xpauli_x(j, k);
          result.zpauli_x(i, k) = result.zpauli_x(i, k) != first.xpauli_x(j, k);
          result.zpauli_z(i, k) = result.zpauli_z(i, k) != first.xpauli_z(j, k);
        }
        zpauli_cphase[i] *= first_x_cphase[j];
      }
      if (second.zpauli_z(i, j)) {
        // add first.zpauli to zpauli
        for (unsigned k = 0; k < n; k++) {
          phase_flip ^= result.zpauli_z(i, k) && first.zpauli_x(j, k);
          result.zpauli_x(i, k) = result.zpauli_x(i, k) != first.zpauli_x(j, k);
          result.zpauli_z(i, k) = result.zpauli_z(i, k) != first.zpauli_z(j, k);
        }
        zpauli_cphase[i] *= first_z_cphase[j];
      }
    }
    if (phase_flip) zpauli_cphase[i] *= -1.;
  }

  std::vector<Complex> current_x_cphase = cphase_from_tableau(
      result.xpauli_x, result.xpauli_z, result.xpauli_phase);
  std::vector<Complex> current_z_cphase = cphase_from_tableau(
      result.zpauli_x, result.zpauli_z, result.zpauli_phase);

  for (unsigned i = 0; i < n; i++) {
    if ((current_x_cphase[i] * xpauli_cphase[i]).imag() != 0)
      throw NotValid("Error in Tableau phase calculations");
    result.xpauli_phase(i) = (current_x_cphase[i] == -xpauli_cphase[i]);
    if ((current_z_cphase[i] * zpauli_cphase[i]).imag() != 0)
      throw NotValid("Error in Tableau phase calculations");
    result.zpauli_phase(i) = (current_z_cphase[i] == -zpauli_cphase[i]);
  }

  return result;
}

}  // namespace tket
