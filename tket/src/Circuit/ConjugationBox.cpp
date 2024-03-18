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

#include "tket/Circuit/ConjugationBox.hpp"

#include <memory>
#include <optional>

#include "Utils/Expression.hpp"
#include "tket/Circuit/Circuit.hpp"
#include "tket/Ops/OpJsonFactory.hpp"
#include "tket/Utils/HelperFunctions.hpp"
#include "tket/Utils/Json.hpp"

namespace tket {

ConjugationBox::ConjugationBox(
    const Op_ptr &compute, const Op_ptr &action,
    const std::optional<Op_ptr> uncompute)
    : Box(OpType::ConjugationBox),
      compute_(compute),
      action_(action),
      uncompute_(uncompute) {
  op_signature_t compute_sig = compute_->get_signature();
  op_signature_t action_sig = action_->get_signature();
  unsigned compute_size_ = compute_sig.size();
  unsigned action_size_ = action_sig.size();
  unsigned compute_n_qubits_ =
      std::count(compute_sig.begin(), compute_sig.end(), EdgeType::Quantum);
  unsigned action_n_qubits_ =
      std::count(action_sig.begin(), action_sig.end(), EdgeType::Quantum);
  unsigned uncompute_size_ = 0, uncompute_n_qubits_ = 0;
  if (uncompute_ != std::nullopt) {
    op_signature_t uncompute_sig = uncompute_.value()->get_signature();
    uncompute_size_ = uncompute_sig.size();
    uncompute_n_qubits_ = std::count(
        uncompute_sig.begin(), uncompute_sig.end(), EdgeType::Quantum);
  }
  if (compute_size_ != compute_n_qubits_ || action_size_ != action_n_qubits_ ||
      uncompute_size_ != uncompute_n_qubits_) {
    throw std::invalid_argument(
        "ConjugationBox only supports quantum operations");
  }
  if (compute_size_ != action_size_ ||
      (uncompute_ != std::nullopt && uncompute_size_ != compute_size_)) {
    throw std::invalid_argument(
        "Operations provided to ConjugationBox need to have the same number of "
        "qubits");
  }
  signature_ = op_signature_t(compute_size_, EdgeType::Quantum);
}

ConjugationBox::ConjugationBox(const ConjugationBox &other)
    : Box(other),
      compute_(other.compute_),
      action_(other.action_),
      uncompute_(other.uncompute_) {}

Op_ptr ConjugationBox::dagger() const {
  return std::make_shared<ConjugationBox>(
      compute_, action_->dagger(), uncompute_);
}

Op_ptr ConjugationBox::transpose() const {
  return std::make_shared<ConjugationBox>(
      (uncompute_ == std::nullopt) ? compute_->dagger()->transpose()
                                   : uncompute_.value()->transpose(),
      action_->transpose(), compute_->transpose());
}

void ConjugationBox::generate_circuit() const {
  Circuit circ(signature_.size());
  std::vector<unsigned> args(circ.n_qubits());
  std::iota(args.begin(), args.end(), 0);
  circ.add_op<unsigned>(compute_, args);
  circ.add_op<unsigned>(action_, args);
  if (uncompute_ != std::nullopt) {
    circ.add_op<unsigned>(uncompute_.value(), args);
  } else {
    circ.add_op<unsigned>(compute_->dagger(), args);
  }
  circ_ = std::make_shared<Circuit>(circ);
}

Op_ptr ConjugationBox::symbol_substitution(
    const SymEngine::map_basic_basic &sub_map) const {
  if (uncompute_.has_value()) {
    return std::make_shared<ConjugationBox>(
        compute_->symbol_substitution(sub_map),
        action_->symbol_substitution(sub_map),
        uncompute_.value()->symbol_substitution(sub_map));
  } else {
    return std::make_shared<ConjugationBox>(
        compute_->symbol_substitution(sub_map),
        action_->symbol_substitution(sub_map));
  }
}

SymSet ConjugationBox::free_symbols() const {
  SymSet compute_syms = compute_->free_symbols();
  SymSet action_syms = action_->free_symbols();
  SymSet uncompute_syms;
  if (uncompute_.has_value()) {
    SymSet s = uncompute_.value()->free_symbols();
    uncompute_syms.insert(s.begin(), s.end());
  }
  SymSet sym_set;
  sym_set.insert(compute_syms.begin(), compute_syms.end());
  sym_set.insert(action_syms.begin(), action_syms.end());
  sym_set.insert(uncompute_syms.begin(), uncompute_syms.end());
  return sym_set;
}

bool ConjugationBox::is_equal(const Op &op_other) const {
  const ConjugationBox &other = dynamic_cast<const ConjugationBox &>(op_other);
  if (id_ == other.get_id()) return true;
  // if only one of them has uncompute_, compare the uncompute_ with
  // the other's compute_.dagger().
  return *compute_ == *other.compute_ && *action_ == *other.action_ &&
         ((uncompute_ == std::nullopt && other.uncompute_ == std::nullopt) ||
          (uncompute_ != std::nullopt && other.uncompute_ != std::nullopt &&
           *uncompute_.value() == *other.uncompute_.value()) ||
          (uncompute_ == std::nullopt && other.uncompute_ != std::nullopt &&
           *compute_->dagger() == *other.uncompute_.value()) ||
          (uncompute_ != std::nullopt && other.uncompute_ == std::nullopt &&
           *uncompute_.value() == *other.compute_->dagger()));
}

nlohmann::json ConjugationBox::to_json(const Op_ptr &op) {
  const auto &box = static_cast<const ConjugationBox &>(*op);
  nlohmann::json j = core_box_json(box);
  j["compute"] = box.get_compute();
  j["action"] = box.get_action();
  // set j["uncompute"] to null
  j["uncompute"] = nlohmann::json();
  std::optional<Op_ptr> uncompute = box.get_uncompute();
  if (uncompute != std::nullopt) {
    j["uncompute"] = uncompute.value();
  }
  return j;
}

Op_ptr ConjugationBox::from_json(const nlohmann::json &j) {
  std::optional<Op_ptr> uncompute = std::nullopt;
  if (j.contains("uncompute") && !j.at("uncompute").is_null()) {
    uncompute = j.at("uncompute").get<Op_ptr>();
  }
  ConjugationBox box = ConjugationBox(
      j.at("compute").get<Op_ptr>(), j.at("action").get<Op_ptr>(), uncompute);
  return set_box_id(
      box,
      boost::lexical_cast<boost::uuids::uuid>(j.at("id").get<std::string>()));
}

REGISTER_OPFACTORY(ConjugationBox, ConjugationBox)

}  // namespace tket
