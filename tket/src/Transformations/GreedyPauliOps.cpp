// Copyright 2019-2024 Cambridge Quantum Computing
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

#include <algorithm>

#include "tket/Circuit/PauliExpBoxes.hpp"
#include "tket/OpType/OpType.hpp"
#include "tket/Transformations/CliffordOptimisation.hpp"
#include "tket/Transformations/GreedyPauliOptimisation.hpp"
#include "tket/Transformations/GreedyPauliOptimisationLookupTables.hpp"

namespace tket {

namespace Transforms {

namespace GreedyPauliSimp {

static std::string pauli_str(const Pauli& p) {
  switch (p) {
    case Pauli::I:
      return "I";
    case Pauli::Z:
      return "Z";
    case Pauli::X:
      return "X";
    case Pauli::Y:
      return "Y";
  }
}

static std::string pauli_vec_str(const std::vector<Pauli>& pvec) {
  std::string s;
  for (unsigned i = 0; i < pvec.size(); i++) {
    s += pauli_str(pvec[i]);
    if (i != pvec.size() - 1) s += "-";
  }
  return s;
}

static std::string pauli_pair_str(
    const std::vector<Pauli>& zvec, const std::vector<Pauli>& xvec) {
  std::string s;
  for (unsigned i = 0; i < zvec.size(); i++) {
    s += pauli_str(zvec[i]);
    s += "|";
    s += pauli_str(xvec[i]);
    if (i != zvec.size() - 1) s += "-";
  }
  return s;
}

std::string TQE::repr() const {
  std::stringstream ss;
  ss << "TQE [";
  switch (type) {
    case TQEType::XX:
      ss << "XX";
      break;
    case TQEType::XY:
      ss << "XY";
      break;
    case TQEType::XZ:
      ss << "XZ";
      break;
    case TQEType::YX:
      ss << "YX";
      break;
    case TQEType::YY:
      ss << "YY";
      break;
    case TQEType::YZ:
      ss << "YZ";
      break;
    case TQEType::ZX:
      ss << "ZX";
      break;
    case TQEType::ZY:
      ss << "ZY";
      break;
    case TQEType::ZZ:
      ss << "ZZ";
  }
  ss << " " << a << " " << b << "]";
  return ss.str();
}

std::string Rotation2Q::repr() const {
  std::stringstream ss;
  ss << "Rotation2Q [" << pauli_str(p_a) << pauli_str(p_b) << angle << " " << a
     << " " << b << "]";
  return ss.str();
}

static CommuteType get_pauli_pair_commute_type(
    const Pauli& p0, const Pauli& p1) {
  if (p0 == Pauli::I && p1 == Pauli::I) {
    return CommuteType::I;
  }
  if (p0 == p1 || p0 == Pauli::I || p1 == Pauli::I) {
    return CommuteType::C;
  }
  return CommuteType::A;
}

// PauliNode abstract class

PauliNode::~PauliNode() {}

void PauliNode::update(const OpType& /*sq_cliff*/, const unsigned& /*a*/) {
  throw GreedyPauliSimpError("Single qubit Clifford update not implemented.");
}

void PauliNode::swap(const unsigned& /*a*/, const unsigned& /*b*/) {
  throw GreedyPauliSimpError("SWAP update not implemented.");
}

// SingleNode

SingleNode::SingleNode(std::vector<Pauli> string, bool sign)
    : string_(string), sign_(sign) {
  if (string.empty()) {
    throw GreedyPauliSimpError("SingleNode cannot have empty strings.");
  }
  weight_ =
      string_.size() - std::count(string_.begin(), string_.end(), Pauli::I);
  if (weight_ == 0) {
    throw GreedyPauliSimpError(
        "SingleNode cannot be constructed with identity strings.");
  }
}

unsigned SingleNode::tqe_cost() const { return weight_ - 1; }

int SingleNode::tqe_cost_increase(const TQE& tqe) const {
  const TQEType& g = tqe.type;
  const unsigned& a = tqe.a;
  const unsigned& b = tqe.b;
  Pauli p0 = string_[a];
  Pauli p1 = string_[b];
  auto [new_p0, new_p1, sign] = TQE_PAULI_MAP.at({g, p0, p1});
  return (p0 == Pauli::I) + (p1 == Pauli::I) - (new_p0 == Pauli::I) -
         (new_p1 == Pauli::I);
}

void SingleNode::update(const TQE& tqe) {
  const TQEType& g = tqe.type;
  const unsigned& a = tqe.a;
  const unsigned& b = tqe.b;
  Pauli p0 = string_[a];
  Pauli p1 = string_[b];
  auto [new_p0, new_p1, sign] = TQE_PAULI_MAP.at({g, p0, p1});
  string_[a] = new_p0;
  string_[b] = new_p1;
  weight_ += (p0 == Pauli::I) + (p1 == Pauli::I) - (new_p0 == Pauli::I) -
             (new_p1 == Pauli::I);
  if (!sign) {
    sign_ = !sign_;
  }
}

std::vector<TQE> SingleNode::reduction_tqes() const {
  std::vector<TQE> tqes;
  // qubits with support
  std::vector<unsigned> sqs;
  for (unsigned i = 0; i < string_.size(); i++) {
    if (string_[i] != Pauli::I) sqs.push_back(i);
  }
  TKET_ASSERT(!sqs.empty());
  for (unsigned i = 0; i < sqs.size() - 1; i++) {
    for (unsigned j = i + 1; j < sqs.size(); j++) {
      std::vector<TQEType> tqe_types =
          TQE_REDUCTION_MAP.at({string_[sqs[i]], string_[sqs[j]]});
      for (const TQEType& tt : tqe_types) {
        tqes.push_back({tt, sqs[i], sqs[j]});
      }
    }
  }
  return tqes;
}

std::pair<unsigned, Pauli> SingleNode::first_support() const {
  for (unsigned i = 0; i < string_.size(); i++) {
    if (string_[i] != Pauli::I) {
      return {i, string_[i]};
    }
  }
  // Should be impossible to reach here
  TKET_ASSERT(false);
}

// ACPairNode

ACPairNode::ACPairNode(
    std::vector<Pauli> z_string, std::vector<Pauli> x_string, bool z_sign,
    bool x_sign)
    : z_string_(z_string),
      x_string_(x_string),
      z_sign_(z_sign),
      x_sign_(x_sign) {
  if (z_string.empty() || x_string.empty()) {
    throw GreedyPauliSimpError("ACPairNode cannot have empty strings.");
  }
  n_commute_entries_ = 0;
  n_anti_commute_entries_ = 0;
  for (unsigned i = 0; i < z_string_.size(); i++) {
    CommuteType commute_type =
        get_pauli_pair_commute_type(z_string_[i], x_string_[i]);
    commute_type_vec_.push_back(commute_type);
    if (commute_type == CommuteType::C) {
      n_commute_entries_ += 1;
    }
    if (commute_type == CommuteType::A) {
      n_anti_commute_entries_ += 1;
    }
  }
}

unsigned ACPairNode::tqe_cost() const {
  // for a node with n A pairs and m C pairs
  // it takes (n-1)/2 TQE gates to convert n-1 A
  // pairs to C pairs. It then takes n-1+m TQE gates
  // to convert all the C pairs to I pairs.
  // total TQEs required is 1.5n - 1.5 + m
  return static_cast<unsigned>(
      1.5 * (n_anti_commute_entries_ - 1) + n_commute_entries_);
}

int ACPairNode::tqe_cost_increase(const TQE& tqe) const {
  const TQEType& g = tqe.type;
  const unsigned& a = tqe.a;
  const unsigned& b = tqe.b;
  Pauli z_p0 = z_string_[a];
  Pauli z_p1 = z_string_[b];
  Pauli x_p0 = x_string_[a];
  Pauli x_p1 = x_string_[b];
  auto [new_z_p0, new_z_p1, z_sign] = TQE_PAULI_MAP.at({g, z_p0, z_p1});
  auto [new_x_p0, new_x_p1, x_sign] = TQE_PAULI_MAP.at({g, x_p0, x_p1});
  CommuteType new_a_type = get_pauli_pair_commute_type(new_z_p0, new_x_p0);
  CommuteType new_b_type = get_pauli_pair_commute_type(new_z_p1, new_x_p1);
  unsigned old_anti_commutes = (commute_type_vec_[a] == CommuteType::A) +
                               (commute_type_vec_[b] == CommuteType::A);
  unsigned old_commutes = (commute_type_vec_[a] == CommuteType::C) +
                          (commute_type_vec_[b] == CommuteType::C);
  unsigned new_anti_commutes =
      (new_a_type == CommuteType::A) + (new_b_type == CommuteType::A);
  unsigned new_commutes =
      (new_a_type == CommuteType::C) + (new_b_type == CommuteType::C);
  int anti_commute_increase = new_anti_commutes - old_anti_commutes;
  int commute_increase = new_commutes - old_commutes;
  return static_cast<int>(1.5 * anti_commute_increase + commute_increase);
}

void ACPairNode::update(const TQE& tqe) {
  const TQEType& g = tqe.type;
  const unsigned& a = tqe.a;
  const unsigned& b = tqe.b;
  Pauli z_p0 = z_string_[a];
  Pauli z_p1 = z_string_[b];
  Pauli x_p0 = x_string_[a];
  Pauli x_p1 = x_string_[b];
  auto [new_z_p0, new_z_p1, z_sign] = TQE_PAULI_MAP.at({g, z_p0, z_p1});
  auto [new_x_p0, new_x_p1, x_sign] = TQE_PAULI_MAP.at({g, x_p0, x_p1});
  CommuteType new_a_type = get_pauli_pair_commute_type(new_z_p0, new_x_p0);
  CommuteType new_b_type = get_pauli_pair_commute_type(new_z_p1, new_x_p1);
  unsigned old_anti_commutes = (commute_type_vec_[a] == CommuteType::A) +
                               (commute_type_vec_[b] == CommuteType::A);
  unsigned old_commutes = (commute_type_vec_[a] == CommuteType::C) +
                          (commute_type_vec_[b] == CommuteType::C);
  unsigned new_anti_commutes =
      (new_a_type == CommuteType::A) + (new_b_type == CommuteType::A);
  unsigned new_commutes =
      (new_a_type == CommuteType::C) + (new_b_type == CommuteType::C);
  int anti_commute_increase = new_anti_commutes - old_anti_commutes;
  int commute_increase = new_commutes - old_commutes;
  n_anti_commute_entries_ += anti_commute_increase;
  n_commute_entries_ += commute_increase;
  commute_type_vec_[a] = new_a_type;
  commute_type_vec_[b] = new_b_type;
  z_string_[a] = new_z_p0;
  z_string_[b] = new_z_p1;
  x_string_[a] = new_x_p0;
  x_string_[b] = new_x_p1;
  if (!z_sign) {
    z_sign_ = !z_sign_;
  }
  if (!x_sign) {
    x_sign_ = !x_sign_;
  }
}

void ACPairNode::update(const OpType& sq_cliff, const unsigned& a) {
  auto [new_z_p, z_sign] = SQ_CLIFF_MAP.at({sq_cliff, z_string_[a]});
  auto [new_x_p, x_sign] = SQ_CLIFF_MAP.at({sq_cliff, x_string_[a]});
  z_string_[a] = new_z_p;
  x_string_[a] = new_x_p;
  if (!z_sign) {
    z_sign_ = !z_sign_;
  }
  if (!x_sign) {
    x_sign_ = !x_sign_;
  }
}

void ACPairNode::swap(const unsigned& a, const unsigned& b) {
  std::swap(z_string_[a], z_string_[b]);
  std::swap(x_string_[a], x_string_[b]);
  std::swap(commute_type_vec_[a], commute_type_vec_[b]);
}

std::vector<TQE> ACPairNode::reduction_tqes() const {
  std::vector<TQE> tqes;
  // qubits with support
  std::vector<unsigned> sqs;
  for (unsigned i = 0; i < commute_type_vec_.size(); i++) {
    if (commute_type_vec_[i] != CommuteType::I) sqs.push_back(i);
  }
  TKET_ASSERT(!sqs.empty());
  for (unsigned i = 0; i < sqs.size() - 1; i++) {
    for (unsigned j = i + 1; j < sqs.size(); j++) {
      std::vector<TQEType> tqe_types;
      unsigned a = sqs[i];
      unsigned b = sqs[j];
      CommuteType ctype0 = commute_type_vec_[a];
      CommuteType ctype1 = commute_type_vec_[b];
      if (ctype0 == CommuteType::A) {
        if (ctype1 == CommuteType::A) {
          // TQEs that transform a AA pair to CC
          tqe_types = AA_TO_CC_MAP.at(
              {z_string_[a], z_string_[b], x_string_[a], x_string_[b]});
        } else {
          // TQEs that transform a AC pair to AI
          tqe_types = AC_TO_AI_MAP.at(
              {z_string_[a], z_string_[b], x_string_[a], x_string_[b]});
        }
      } else {
        if (ctype1 == CommuteType::A) {
          // TQEs that transform a CA pair to a IA
          tqe_types = AC_TO_AI_MAP.at(
              {z_string_[b], z_string_[a], x_string_[b], x_string_[a]});
          // flip qubits
          a = sqs[j];
          b = sqs[i];
        } else {
          // TQEs that transform a CC pair to CI or IC, not always
          // possible
          tqe_types = CC_TO_IC_OR_CI_MAP.at(
              {z_string_[a], z_string_[b], x_string_[a], x_string_[b]});
        }
      }
      for (const TQEType& tt : tqe_types) {
        tqes.push_back({tt, a, b});
      }
    }
  }
  return tqes;
}

std::tuple<unsigned, Pauli, Pauli> ACPairNode::first_support() const {
  for (unsigned i = 0; i < commute_type_vec_.size(); i++) {
    if (commute_type_vec_[i] != CommuteType::I) {
      return {i, z_string_[i], x_string_[i]};
    }
  }
  // Should be impossible to reach here
  TKET_ASSERT(false);
}

// PauliRotation
PauliRotation::PauliRotation(std::vector<Pauli> string, bool sign, Expr theta)
    : SingleNode(string, sign), theta_(theta) {}

CommuteInfo PauliRotation::get_commute_info() const { return {{string_}, {}}; }

std::string PauliRotation::repr() const {
  std::stringstream ss;
  ss << "PauliRotation [" << pauli_vec_str(string_) << "," << this->angle()
     << "]";
  return ss.str();
}

ConditionalBlock::ConditionalBlock(
    std::vector<std::tuple<std::vector<Pauli>, bool, Expr>> rotations,
    std::vector<unsigned> cond_bits, unsigned cond_value)
    : rotations_(rotations), cond_bits_(cond_bits), cond_value_(cond_value) {
  total_weight_ = 0;
  for (const std::tuple<std::vector<Pauli>, bool, Expr>& rot : rotations_) {
    total_weight_ +=
        std::get<0>(rot).size() -
        std::count(std::get<0>(rot).begin(), std::get<0>(rot).end(), Pauli::I);
  }
}

unsigned ConditionalBlock::tqe_cost() const { return total_weight_ - 1; }

int ConditionalBlock::tqe_cost_increase(const TQE& tqe) const {
  int total_increase = 0;
  for (const std::tuple<std::vector<Pauli>, bool, Expr>& rot : rotations_) {
    const TQEType& g = tqe.type;
    const unsigned& a = tqe.a;
    const unsigned& b = tqe.b;
    Pauli p0 = std::get<0>(rot)[a];
    Pauli p1 = std::get<0>(rot)[b];
    auto [new_p0, new_p1, sign] = TQE_PAULI_MAP.at({g, p0, p1});
    total_increase += (p0 == Pauli::I) + (p1 == Pauli::I) -
                      (new_p0 == Pauli::I) - (new_p1 == Pauli::I);
  }
  return total_increase;
}

void ConditionalBlock::update(const TQE& tqe) {
  for (std::tuple<std::vector<Pauli>, bool, Expr>& rot : rotations_) {
    const TQEType& g = tqe.type;
    const unsigned& a = tqe.a;
    const unsigned& b = tqe.b;
    Pauli p0 = std::get<0>(rot)[a];
    Pauli p1 = std::get<0>(rot)[b];
    auto [new_p0, new_p1, sign] = TQE_PAULI_MAP.at({g, p0, p1});
    std::get<0>(rot)[a] = new_p0;
    std::get<0>(rot)[b] = new_p1;
    total_weight_ += (p0 == Pauli::I) + (p1 == Pauli::I) -
                     (new_p0 == Pauli::I) - (new_p1 == Pauli::I);
    if (!sign) {
      std::get<1>(rot) = !std::get<1>(rot);
    }
  }
}

CommuteInfo ConditionalBlock::get_commute_info() const {
  std::vector<std::pair<UnitID, BitType>> bits_info;
  for (unsigned b : cond_bits_) {
    bits_info.push_back({Bit(b), BitType::READ});
  }
  std::vector<std::vector<Pauli>> strings;
  for (const std::tuple<std::vector<Pauli>, bool, Expr>& rot : rotations_) {
    strings.push_back(std::get<0>(rot));
  }
  return {strings, bits_info};
}

void ConditionalBlock::append(const ConditionalBlock& other) {
  for (const auto& rot : other.rotations()) {
    rotations_.push_back(rot);
    total_weight_ +=
        std::get<0>(rot).size() -
        std::count(std::get<0>(rot).begin(), std::get<0>(rot).end(), Pauli::I);
  }
}

std::string ConditionalBlock::repr() const { return "ConditionalBlock []"; }

// PauliPropagation
PauliPropagation::PauliPropagation(
    std::vector<Pauli> z_string, std::vector<Pauli> x_string, bool z_sign,
    bool x_sign, unsigned qubit_index)
    : ACPairNode(z_string, x_string, z_sign, x_sign),
      qubit_index_(qubit_index) {}

CommuteInfo PauliPropagation::get_commute_info() const {
  return {{z_string_, x_string_}, {}};
}

std::string PauliPropagation::repr() const {
  return "PauliPropagation [" + pauli_pair_str(z_string_, x_string_) + "]";
}

// ClassicalNode
ClassicalNode::ClassicalNode(std::vector<UnitID> args, Op_ptr op)
    : args_(args), op_(op) {}

CommuteInfo ClassicalNode::get_commute_info() const {
  std::vector<std::pair<UnitID, BitType>> bits_info;
  for (const UnitID& b : args_) {
    bits_info.push_back({b, BitType::WRITE});
  }
  return {{}, bits_info};
}

std::string ClassicalNode::repr() const { return "Classical []"; }

// MidMeasure
MidMeasure::MidMeasure(std::vector<Pauli> string, bool sign, unsigned bit)
    : SingleNode(string, sign), bit_(bit) {}

CommuteInfo MidMeasure::get_commute_info() const {
  return {{string_}, {{Bit(bit_), BitType::WRITE}}};
}

std::string MidMeasure::repr() const {
  return "MidMeasure [" + pauli_vec_str(string_) + "]";
}

// Reset
Reset::Reset(
    std::vector<Pauli> z_string, std::vector<Pauli> x_string, bool z_sign,
    bool x_sign)
    : ACPairNode(z_string, x_string, z_sign, x_sign) {}

CommuteInfo Reset::get_commute_info() const {
  return {{z_string_, x_string_}, {}};
}

std::string Reset::repr() const {
  return "Reset [" + pauli_pair_str(z_string_, x_string_) + "]";
}

}  // namespace GreedyPauliSimp

}  // namespace Transforms

}  // namespace tket
