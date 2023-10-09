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

#include "tket/Circuit/PauliExpBoxes.hpp"

#include <iostream>

#include "tket/Circuit/CircUtils.hpp"
#include "tket/Converters/PauliGadget.hpp"
#include "tket/Converters/PhasePoly.hpp"
#include "tket/Diagonalisation/Diagonalisation.hpp"
#include "tket/Ops/OpJsonFactory.hpp"

namespace tket {

PauliExpBox::PauliExpBox(
    const SymPauliTensor &paulis, CXConfigType cx_config_type)
    : Box(OpType::PauliExpBox,
          op_signature_t(paulis.size(), EdgeType::Quantum)),
      paulis_(paulis),
      cx_config_(cx_config_type) {}

PauliExpBox::PauliExpBox(const PauliExpBox &other)
    : Box(other), paulis_(other.paulis_), cx_config_(other.cx_config_) {}

PauliExpBox::PauliExpBox() : PauliExpBox({{}, 0}) {}

bool PauliExpBox::is_clifford() const {
  return equiv_0(4 * paulis_.coeff) || paulis_.string.empty();
}

SymSet PauliExpBox::free_symbols() const { return paulis_.free_symbols(); }

Op_ptr PauliExpBox::dagger() const {
  return std::make_shared<PauliExpBox>(
      SymPauliTensor(paulis_.string, -paulis_.coeff), cx_config_);
}

Op_ptr PauliExpBox::transpose() const {
  SymPauliTensor tr = paulis_;
  tr.transpose();
  return std::make_shared<PauliExpBox>(tr, cx_config_);
}

Op_ptr PauliExpBox::symbol_substitution(
    const SymEngine::map_basic_basic &sub_map) const {
  return std::make_shared<PauliExpBox>(
      this->paulis_.symbol_substitution(sub_map), this->cx_config_);
}

void PauliExpBox::generate_circuit() const {
  Circuit circ = pauli_gadget(paulis_.string, paulis_.coeff, cx_config_);
  circ_ = std::make_shared<Circuit>(circ);
}

bool PauliExpBox::is_equal(const Op &op_other) const {
  const PauliExpBox &other = dynamic_cast<const PauliExpBox &>(op_other);
  if (id_ == other.get_id()) return true;
  return cx_config_ == other.cx_config_ && paulis_.equiv_mod(other.paulis_, 4);
}

nlohmann::json PauliExpBox::to_json(const Op_ptr &op) {
  const auto &box = static_cast<const PauliExpBox &>(*op);
  nlohmann::json j = core_box_json(box);
  j["paulis"] = box.get_paulis();
  j["phase"] = box.get_phase();
  j["cx_config"] = box.get_cx_config();
  return j;
}

Op_ptr PauliExpBox::from_json(const nlohmann::json &j) {
  PauliExpBox box = PauliExpBox(
      SymPauliTensor(
          j.at("paulis").get<std::vector<Pauli>>(), j.at("phase").get<Expr>()),
      j.at("cx_config").get<CXConfigType>());
  return set_box_id(
      box,
      boost::lexical_cast<boost::uuids::uuid>(j.at("id").get<std::string>()));
}

REGISTER_OPFACTORY(PauliExpBox, PauliExpBox)

PauliExpPairBox::PauliExpPairBox(
    const SymPauliTensor &paulis0, const SymPauliTensor &paulis1,
    CXConfigType cx_config_type)
    : Box(OpType::PauliExpPairBox,
          op_signature_t(paulis0.size(), EdgeType::Quantum)),
      paulis0_(paulis0),
      paulis1_(paulis1),
      cx_config_(cx_config_type) {
  if (paulis0.size() != paulis1.size()) {
    throw PauliExpBoxInvalidity(
        "Pauli strings within PauliExpPairBox must be of same length (pad with "
        "identities if necessary)");
  }
}

PauliExpPairBox::PauliExpPairBox(const PauliExpPairBox &other)
    : Box(other),
      paulis0_(other.paulis0_),
      paulis1_(other.paulis1_),
      cx_config_(other.cx_config_) {}

PauliExpPairBox::PauliExpPairBox() : PauliExpPairBox({{}, 0}, {{}, 0}) {}

bool PauliExpPairBox::is_clifford() const {
  auto is_clifford0 = equiv_0(4 * paulis0_.coeff) || paulis0_.string.empty();
  auto is_clifford1 = equiv_0(4 * paulis1_.coeff) || paulis1_.string.empty();
  return is_clifford0 && is_clifford1;
}

SymSet PauliExpPairBox::free_symbols() const {
  return expr_free_symbols({paulis0_.coeff, paulis1_.coeff});
}

Op_ptr PauliExpPairBox::dagger() const {
  return std::make_shared<PauliExpPairBox>(
      SymPauliTensor(paulis1_.string, -paulis1_.coeff),
      SymPauliTensor(paulis0_.string, -paulis0_.coeff), cx_config_);
}

Op_ptr PauliExpPairBox::transpose() const {
  SymPauliTensor tr0 = paulis0_;
  tr0.transpose();
  SymPauliTensor tr1 = paulis1_;
  tr1.transpose();
  return std::make_shared<PauliExpPairBox>(tr1, tr0, cx_config_);
}

Op_ptr PauliExpPairBox::symbol_substitution(
    const SymEngine::map_basic_basic &sub_map) const {
  return std::make_shared<PauliExpPairBox>(
      this->paulis0_.symbol_substitution(sub_map),
      this->paulis1_.symbol_substitution(sub_map), this->cx_config_);
}

void PauliExpPairBox::generate_circuit() const {
  Circuit circ = Circuit(paulis0_.size());
  append_pauli_gadget_pair(circ, paulis0_, paulis1_, cx_config_);
  circ_ = std::make_shared<Circuit>(circ);
}

bool PauliExpPairBox::is_equal(const Op &op_other) const {
  const PauliExpPairBox &other =
      dynamic_cast<const PauliExpPairBox &>(op_other);
  if (id_ == other.get_id()) return true;
  return cx_config_ == other.cx_config_ &&
         paulis0_.equiv_mod(other.paulis0_, 4) &&
         paulis1_.equiv_mod(other.paulis1_, 4);
}

nlohmann::json PauliExpPairBox::to_json(const Op_ptr &op) {
  const auto &box = static_cast<const PauliExpPairBox &>(*op);
  nlohmann::json j = core_box_json(box);
  j["paulis_pair"] = box.get_paulis_pair();
  j["phase_pair"] = box.get_phase_pair();
  j["cx_config"] = box.get_cx_config();
  return j;
}

Op_ptr PauliExpPairBox::from_json(const nlohmann::json &j) {
  const auto [paulis0, paulis1] =
      j.at("paulis_pair")
          .get<std::tuple<std::vector<Pauli>, std::vector<Pauli>>>();
  const auto [phase0, phase1] =
      j.at("phase_pair").get<std::tuple<Expr, Expr>>();
  PauliExpPairBox box = PauliExpPairBox(
      SymPauliTensor(paulis0, phase0), SymPauliTensor(paulis1, phase1),
      j.at("cx_config").get<CXConfigType>());
  return set_box_id(
      box,
      boost::lexical_cast<boost::uuids::uuid>(j.at("id").get<std::string>()));
}

REGISTER_OPFACTORY(PauliExpPairBox, PauliExpPairBox)

PauliExpCommutingSetBox::PauliExpCommutingSetBox(
    const std::vector<SymPauliTensor> &pauli_gadgets,
    CXConfigType cx_config_type)
    : Box(OpType::PauliExpCommutingSetBox),
      pauli_gadgets_(pauli_gadgets),
      cx_config_(cx_config_type) {
  // check at least one gadget
  if (pauli_gadgets.empty()) {
    throw PauliExpBoxInvalidity(
        "PauliExpCommutingSetBox requires at least one Pauli string");
  }
  // check all gadgets have same Pauli string length
  auto n_qubits = pauli_gadgets[0].size();
  for (const auto &gadget : pauli_gadgets) {
    if (gadget.size() != n_qubits) {
      throw PauliExpBoxInvalidity(
          "the Pauli strings within PauliExpCommutingSetBox must all be the "
          "same length");
    }
  }
  if (!this->paulis_commute()) {
    throw PauliExpBoxInvalidity(
        "Pauli strings used to define PauliExpCommutingSetBox must all "
        "commute");
  }
  signature_ = op_signature_t(n_qubits, EdgeType::Quantum);
}

PauliExpCommutingSetBox::PauliExpCommutingSetBox(
    const PauliExpCommutingSetBox &other)
    : Box(other),
      pauli_gadgets_(other.pauli_gadgets_),
      cx_config_(other.cx_config_) {}

PauliExpCommutingSetBox::PauliExpCommutingSetBox()
    : PauliExpCommutingSetBox({{{}, 0}}) {}

bool PauliExpCommutingSetBox::is_clifford() const {
  return std::all_of(
      pauli_gadgets_.begin(), pauli_gadgets_.end(),
      [](const SymPauliTensor &pauli_exp) {
        return equiv_0(4 * pauli_exp.coeff) || pauli_exp.string.empty();
      });
}

bool PauliExpCommutingSetBox::paulis_commute() const {
  for (auto string0 = pauli_gadgets_.begin(); string0 != pauli_gadgets_.end();
       ++string0) {
    for (auto string1 = std::next(string0); string1 != pauli_gadgets_.end();
         ++string1) {
      if (!string0->commutes_with(*string1)) {
        return false;
      }
    }
  }
  return true;
}

SymSet PauliExpCommutingSetBox::free_symbols() const {
  std::vector<Expr> angles;
  for (const auto &pauli_exp : pauli_gadgets_) {
    angles.push_back(pauli_exp.coeff);
  }
  return expr_free_symbols(angles);
}

Op_ptr PauliExpCommutingSetBox::dagger() const {
  std::vector<SymPauliTensor> dagger_gadgets;
  for (const auto &pauli_exp : pauli_gadgets_) {
    dagger_gadgets.emplace_back(pauli_exp.string, -pauli_exp.coeff);
  }
  return std::make_shared<PauliExpCommutingSetBox>(dagger_gadgets, cx_config_);
}

Op_ptr PauliExpCommutingSetBox::transpose() const {
  std::vector<SymPauliTensor> transpose_gadgets;
  for (const auto &pauli_exp : pauli_gadgets_) {
    SymPauliTensor tr = pauli_exp;
    tr.transpose();
    transpose_gadgets.push_back(tr);
  }
  return std::make_shared<PauliExpCommutingSetBox>(
      transpose_gadgets, cx_config_);
}

Op_ptr PauliExpCommutingSetBox::symbol_substitution(
    const SymEngine::map_basic_basic &sub_map) const {
  std::vector<SymPauliTensor> symbol_sub_gadgets;
  for (const auto &pauli_exp : pauli_gadgets_) {
    symbol_sub_gadgets.push_back(pauli_exp.symbol_substitution(sub_map));
  }
  return std::make_shared<PauliExpCommutingSetBox>(
      symbol_sub_gadgets, this->cx_config_);
}

void PauliExpCommutingSetBox::generate_circuit() const {
  Circuit circ = Circuit(pauli_gadgets_[0].size());

  std::list<SpSymPauliTensor> gadgets;
  for (const auto &pauli_gadget : pauli_gadgets_) {
    gadgets.push_back((SpSymPauliTensor)pauli_gadget);
  }
  std::set<Qubit> qubits;
  for (unsigned i = 0; i < pauli_gadgets_[0].size(); i++)
    qubits.insert(Qubit(i));

  Circuit cliff_circ = mutual_diagonalise(gadgets, qubits, cx_config_);
  circ.append(cliff_circ);

  Circuit phase_poly_circ = Circuit(pauli_gadgets_[0].size());

  for (const SpSymPauliTensor &pgp : gadgets) {
    append_single_pauli_gadget(phase_poly_circ, pgp);
  }
  PhasePolyBox ppbox(phase_poly_circ);
  Circuit after_synth_circ = *ppbox.to_circuit();

  circ.append(after_synth_circ);
  circ.append(cliff_circ.dagger());

  circ_ = std::make_shared<Circuit>(circ);
}

bool PauliExpCommutingSetBox::is_equal(const Op &op_other) const {
  const PauliExpCommutingSetBox &other =
      dynamic_cast<const PauliExpCommutingSetBox &>(op_other);
  if (id_ == other.get_id()) return true;
  if (cx_config_ != other.cx_config_) return false;
  return std::equal(
      pauli_gadgets_.begin(), pauli_gadgets_.end(),
      other.pauli_gadgets_.begin(), other.pauli_gadgets_.end(),
      [](const SymPauliTensor &a, const SymPauliTensor &b) {
        return a.equiv_mod(b, 4);
      });
}

nlohmann::json PauliExpCommutingSetBox::to_json(const Op_ptr &op) {
  const auto &box = static_cast<const PauliExpCommutingSetBox &>(*op);
  nlohmann::json j = core_box_json(box);
  j["pauli_gadgets"] = box.get_pauli_gadgets();
  j["cx_config"] = box.get_cx_config();
  return j;
}

Op_ptr PauliExpCommutingSetBox::from_json(const nlohmann::json &j) {
  PauliExpCommutingSetBox box = PauliExpCommutingSetBox(
      j.at("pauli_gadgets").get<std::vector<SymPauliTensor>>(),
      j.at("cx_config").get<CXConfigType>());
  return set_box_id(
      box,
      boost::lexical_cast<boost::uuids::uuid>(j.at("id").get<std::string>()));
}

REGISTER_OPFACTORY(PauliExpCommutingSetBox, PauliExpCommutingSetBox)

}  // namespace tket
