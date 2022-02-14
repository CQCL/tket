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

#include "ClassicalOps.hpp"

#include "OpType/OpType.hpp"
#include "Utils/Json.hpp"

namespace tket {

static uint32_t u32_from_boolvec(const std::vector<bool> &x) {
  unsigned n = x.size();
  if (n > 32) {
    throw std::domain_error("Vector of bool exceeds maximum size (32)");
  }
  uint32_t X = 0;
  for (unsigned i = 0; i < n; i++) {
    if (x[i]) X |= (1u << i);
  }
  return X;
}

static nlohmann::json classical_to_json(const Op_ptr &op, const OpType &type) {
  nlohmann::json j_class;
  switch (type) {
    case OpType::MultiBit: {
      const auto &multibit = static_cast<const MultiBitOp &>(*op);
      j_class["op"] = multibit.get_op();
      j_class["n"] = multibit.get_n();
      return j_class;
    }
    case OpType::RangePredicate: {
      const auto &rangop = static_cast<const RangePredicateOp &>(*op);
      j_class["lower"] = rangop.lower();
      j_class["upper"] = rangop.upper();
      j_class["n_i"] = rangop.get_n_i();
      return j_class;
    }
    case OpType::ExplicitModifier: {
      const auto &expmodop = static_cast<const ExplicitModifierOp &>(*op);
      j_class["n_i"] = expmodop.get_n_i();
      j_class["values"] = expmodop.get_values();
      j_class["name"] = expmodop.get_name();
      return j_class;
    }
    case OpType::ExplicitPredicate: {
      const auto &exppredop = static_cast<const ExplicitPredicateOp &>(*op);
      j_class["n_i"] = exppredop.get_n_i();
      j_class["values"] = exppredop.get_values();
      j_class["name"] = exppredop.get_name();
      return j_class;
    }
    case OpType::ClassicalTransform: {
      const auto &classtop = static_cast<const ClassicalTransformOp &>(*op);
      j_class["n_io"] = classtop.get_n_io();
      j_class["values"] = classtop.get_values();
      j_class["name"] = classtop.get_name();
      return j_class;
    }
    case OpType::SetBits: {
      const auto &setop = static_cast<const SetBitsOp &>(*op);
      j_class["values"] = setop.get_values();
      return j_class;
    }
    case OpType::CopyBits: {
      const auto &cop = static_cast<const CopyBitsOp &>(*op);
      j_class["n_i"] = cop.get_n_i();
      return j_class;
    }
    default:
      throw JsonError(
          "Classical op with type " + optypeinfo().at(type).name +
          " cannot be serialized.");
  }
}

static std::shared_ptr<ClassicalOp> classical_from_json(
    const nlohmann::json &j_class, const OpType &type) {
  switch (type) {
    case OpType::MultiBit:
      return std::make_shared<MultiBitOp>(
          classical_from_json(
              j_class.at("op").at("classical"),
              j_class.at("op").at("type").get<OpType>()),
          j_class.at("n").get<unsigned>());
    case OpType::RangePredicate:
      return std::make_shared<RangePredicateOp>(
          j_class.at("n_i").get<unsigned>(),
          j_class.at("lower").get<unsigned>(),
          j_class.at("upper").get<unsigned>());
    case OpType::CopyBits:
      return std::make_shared<CopyBitsOp>(j_class.at("n_i").get<unsigned>());
    case OpType::SetBits:
      return std::make_shared<SetBitsOp>(
          j_class.at("values").get<std::vector<bool>>());
    case OpType::ExplicitModifier:
      return std::make_shared<ExplicitModifierOp>(
          j_class.at("n_i").get<unsigned>(),
          j_class.at("values").get<std::vector<bool>>(),
          j_class.at("name").get<std::string>());
    case OpType::ExplicitPredicate:
      return std::make_shared<ExplicitPredicateOp>(
          j_class.at("n_i").get<unsigned>(),
          j_class.at("values").get<std::vector<bool>>(),
          j_class.at("name").get<std::string>());
    case OpType::ClassicalTransform:
      return std::make_shared<ClassicalTransformOp>(
          j_class.at("n_io").get<unsigned>(),
          j_class.at("values").get<std::vector<uint32_t>>(),
          j_class.at("name").get<std::string>());
    default:
      throw JsonError(
          "Classical op with name " + j_class.at("name").get<std::string>() +
          " cannot be deserialized.");
  }
}

ClassicalOp::ClassicalOp(
    OpType type, unsigned n_i, unsigned n_io, unsigned n_o,
    const std::string &name)
    : Op(type), n_i_(n_i), n_io_(n_io), n_o_(n_o), name_(name), sig_() {
  for (unsigned i = 0; i < n_i; i++) {
    sig_.push_back(EdgeType::Boolean);
  }
  for (unsigned j = 0; j < n_io + n_o; j++) {
    sig_.push_back(EdgeType::Classical);
  }
}

nlohmann::json ClassicalOp::serialize() const {
  nlohmann::json j;
  j["type"] = get_type();
  j["classical"] = classical_to_json(shared_from_this(), get_type());
  return j;
}

Op_ptr ClassicalOp::deserialize(const nlohmann::json &j) {
  return classical_from_json(j.at("classical"), j.at("type").get<OpType>());
}

std::string ClassicalOp::get_name(bool) const { return name_; }

bool ClassicalOp::is_equal(const Op &op_other) const {
  const ClassicalOp &other = dynamic_cast<const ClassicalOp &>(op_other);

  if (n_i_ != other.get_n_i()) return false;
  if (n_io_ != other.get_n_io()) return false;
  if (n_o_ != other.get_n_o()) return false;
  unsigned N = n_i_ + n_io_;
  uint32_t xlim = 1u << N;
  std::vector<bool> v(N);
  for (uint32_t x = 0; x < xlim; x++) {
    for (unsigned i = 0; i < N; i++) {
      v[i] = (x >> i) & 1;
    }
    if (eval(v) != other.eval(v)) return false;
  }
  return true;
}

ClassicalTransformOp::ClassicalTransformOp(
    unsigned n, const std::vector<uint32_t> &values, const std::string &name)
    : ClassicalOp(OpType::ClassicalTransform, 0, n, 0, name), values_(values) {
  if (n > 32) {
    throw std::domain_error("Too many inputs/outputs (maximum is 32)");
  }
}

std::vector<bool> ClassicalTransformOp::eval(const std::vector<bool> &x) const {
  if (x.size() != n_io_) {
    throw std::domain_error("Incorrect input size");
  }
  uint32_t val = values_[u32_from_boolvec(x)];
  std::vector<bool> y(n_io_);
  for (unsigned j = 0; j < n_io_; j++) {
    y[j] = (val >> j) & 1;
  }
  return y;
}

std::string SetBitsOp::get_name(bool) const {
  std::stringstream name;
  name << name_ << "(";
  for (auto v : values_) {
    name << v;  // "0" or "1"
  }
  name << ")";
  return name.str();
}

std::vector<bool> SetBitsOp::eval(const std::vector<bool> &x) const {
  if (!x.empty()) {
    throw std::domain_error("Non-empty input");
  }
  return values_;
}

std::vector<bool> CopyBitsOp::eval(const std::vector<bool> &x) const {
  if (x.size() != n_i_) {
    throw std::domain_error("Incorrect input size");
  }
  return x;
}

std::string RangePredicateOp::get_name(bool) const {
  std::stringstream name;
  name << name_ << "([" << a << "," << b << "])";
  return name.str();
}

std::vector<bool> RangePredicateOp::eval(const std::vector<bool> &x) const {
  if (x.size() != n_i_) {
    throw std::domain_error("Incorrect input size");
  }
  uint32_t X = u32_from_boolvec(x);
  std::vector<bool> y(1);
  y[0] = (X >= a && X <= b);
  return y;
}

bool RangePredicateOp::is_equal(const Op &op_other) const {
  const RangePredicateOp &other =
      dynamic_cast<const RangePredicateOp &>(op_other);

  if (n_i_ != other.n_i_) return false;
  return (a == other.a && b == other.b);
}

ExplicitPredicateOp::ExplicitPredicateOp(
    unsigned n, const std::vector<bool> &values, const std::string &name)
    : PredicateOp(OpType::ExplicitPredicate, n, name), values_(values) {
  if (n > 32) {
    throw std::domain_error("Too many inputs");
  }
}

std::vector<bool> ExplicitPredicateOp::eval(const std::vector<bool> &x) const {
  if (x.size() != n_i_) {
    throw std::domain_error("Incorrect input size");
  }
  std::vector<bool> y(1);
  y[0] = values_[u32_from_boolvec(x)];
  return y;
}

ExplicitModifierOp::ExplicitModifierOp(
    unsigned n, const std::vector<bool> &values, const std::string &name)
    : ModifyingOp(OpType::ExplicitModifier, n, name), values_(values) {
  if (n > 31) {
    throw std::domain_error("Too many inputs");
  }
}

std::vector<bool> ExplicitModifierOp::eval(const std::vector<bool> &x) const {
  if (x.size() != n_i_ + 1) {
    throw std::domain_error("Incorrect input size");
  }
  std::vector<bool> y(1);
  y[0] = values_[u32_from_boolvec(x)];
  return y;
}

MultiBitOp::MultiBitOp(std::shared_ptr<const ClassicalOp> op, unsigned n)
    : ClassicalOp(
          OpType::MultiBit, n * op->get_n_i(), n * op->get_n_io(),
          n * op->get_n_o(), op->get_name()),
      op_(op),
      n_(n) {
  op_signature_t op_sig = op->get_signature();
  const std::size_t op_sig_size = op_sig.size();
  sig_.clear();
  sig_.reserve(n_ * op_sig.size());
  for (unsigned i = 0; i < n_; i++) {
    std::copy_n(op_sig.begin(), op_sig_size, std::back_inserter(sig_));
  }
}

std::string MultiBitOp::get_name(bool) const {
  std::stringstream name;
  name << name_ << " (*" << n_ << ")";
  return name.str();
}

std::vector<bool> MultiBitOp::eval(const std::vector<bool> &x) const {
  if (x.size() != n_i_ + n_io_) {
    throw std::domain_error("Incorrect input size");
  }
  unsigned n_op_inputs = op_->get_n_i() + op_->get_n_io();
  unsigned n_op_outputs = op_->get_n_io() + op_->get_n_o();
  std::vector<bool> y(n_io_ + n_o_);
  for (unsigned i = 0; i < n_; i++) {
    std::vector<bool> x_i(n_op_inputs);
    for (unsigned j = 0; j < n_op_inputs; j++) {
      x_i[j] = x[n_op_inputs * i + j];
    }
    std::vector<bool> y_i = op_->eval(x_i);
    for (unsigned j = 0; j < n_op_outputs; j++) {
      y[n_op_outputs * i + j] = y_i[j];
    }
  }
  return y;
}

bool MultiBitOp::is_equal(const Op &op_other) const {
  const MultiBitOp &other = dynamic_cast<const MultiBitOp &>(op_other);

  if (n_ != other.n_) return false;
  return (*op_ == *(other.op_));
}

std::shared_ptr<ClassicalTransformOp> ClassicalX() {
  static const std::vector<uint32_t> values = {1, 0};
  static const std::shared_ptr<ClassicalTransformOp> op =
      std::make_shared<ClassicalTransformOp>(1, values, "ClassicalX");
  return op;
}

std::shared_ptr<ClassicalTransformOp> ClassicalCX() {
  static const std::vector<uint32_t> values = {0, 3, 2, 1};
  static const std::shared_ptr<ClassicalTransformOp> op =
      std::make_shared<ClassicalTransformOp>(2, values, "ClassicalCX");
  return op;
}

std::shared_ptr<ExplicitPredicateOp> NotOp() {
  static const std::vector<bool> values = {1, 0};
  static const std::shared_ptr<ExplicitPredicateOp> op =
      std::make_shared<ExplicitPredicateOp>(1, values, "NOT");
  return op;
}

std::shared_ptr<ExplicitPredicateOp> AndOp() {
  static const std::vector<bool> values = {0, 0, 0, 1};
  static const std::shared_ptr<ExplicitPredicateOp> op =
      std::make_shared<ExplicitPredicateOp>(2, values, "AND");
  return op;
}

std::shared_ptr<ExplicitPredicateOp> OrOp() {
  static const std::vector<bool> values = {0, 1, 1, 1};
  static const std::shared_ptr<ExplicitPredicateOp> op =
      std::make_shared<ExplicitPredicateOp>(2, values, "OR");
  return op;
}

std::shared_ptr<ExplicitPredicateOp> XorOp() {
  static const std::vector<bool> values = {0, 1, 1, 0};
  static const std::shared_ptr<ExplicitPredicateOp> op =
      std::make_shared<ExplicitPredicateOp>(2, values, "XOR");
  return op;
}

std::shared_ptr<ExplicitModifierOp> AndWithOp() {
  static const std::vector<bool> values = {0, 0, 0, 1};
  static const std::shared_ptr<ExplicitModifierOp> op =
      std::make_shared<ExplicitModifierOp>(1, values, "AND");
  return op;
}

std::shared_ptr<ExplicitModifierOp> OrWithOp() {
  static const std::vector<bool> values = {0, 1, 1, 1};
  static const std::shared_ptr<ExplicitModifierOp> op =
      std::make_shared<ExplicitModifierOp>(1, values, "OR");
  return op;
}

std::shared_ptr<ExplicitModifierOp> XorWithOp() {
  static const std::vector<bool> values = {0, 1, 1, 0};
  static const std::shared_ptr<ExplicitModifierOp> op =
      std::make_shared<ExplicitModifierOp>(1, values, "XOR");
  return op;
}

}  // namespace tket
