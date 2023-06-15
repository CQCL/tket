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
    const std::vector<Pauli> &paulis, const Expr &t,
    CXConfigType cx_config_type)
    : Box(OpType::PauliExpBox,
          op_signature_t(paulis.size(), EdgeType::Quantum)),
      paulis_(paulis),
      t_(t),
      cx_config_(cx_config_type) {}

PauliExpBox::PauliExpBox(const PauliExpBox &other)
    : Box(other),
      paulis_(other.paulis_),
      t_(other.t_),
      cx_config_(other.cx_config_) {}

PauliExpBox::PauliExpBox() : PauliExpBox({}, 0.) {}

bool PauliExpBox::is_clifford() const {
  return equiv_0(4 * t_) || paulis_.empty();
}

SymSet PauliExpBox::free_symbols() const { return expr_free_symbols(t_); }

Op_ptr PauliExpBox::dagger() const {
  return std::make_shared<PauliExpBox>(paulis_, -t_, cx_config_);
}

Op_ptr PauliExpBox::transpose() const {
  std::vector<Pauli> paulis = get_paulis();
  int number_y_pauli_mod2 =
      std::count(paulis.begin(), paulis.end(), Pauli::Y) % 2;
  // Negate the parameter if odd nr of Paulis (mult with 1 if mod2=0, -1 if
  // mod2=1)
  int t_fac = -number_y_pauli_mod2 * 2 + 1;
  return std::make_shared<PauliExpBox>(paulis_, t_fac * t_, cx_config_);
}

Op_ptr PauliExpBox::symbol_substitution(
    const SymEngine::map_basic_basic &sub_map) const {
  return std::make_shared<PauliExpBox>(
      this->paulis_, this->t_.subs(sub_map), this->cx_config_);
}

void PauliExpBox::generate_circuit() const {
  Circuit circ = pauli_gadget(paulis_, t_, cx_config_);
  circ_ = std::make_shared<Circuit>(circ);
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
      j.at("paulis").get<std::vector<Pauli>>(), j.at("phase").get<Expr>(),
      j.at("cx_config").get<CXConfigType>());
  return set_box_id(
      box,
      boost::lexical_cast<boost::uuids::uuid>(j.at("id").get<std::string>()));
}

REGISTER_OPFACTORY(PauliExpBox, PauliExpBox)

PauliExpPairBox::PauliExpPairBox(
    const std::vector<Pauli> &paulis0, const Expr &t0,
    const std::vector<Pauli> &paulis1, const Expr &t1,
    CXConfigType cx_config_type)
    : Box(OpType::PauliExpPairBox,
          op_signature_t(paulis0.size(), EdgeType::Quantum)),
      paulis0_(paulis0),
      t0_(std::move(t0)),
      paulis1_(paulis1),
      t1_(std::move(t1)),
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
      t0_(other.t0_),
      paulis1_(other.paulis1_),
      t1_(other.t1_),
      cx_config_(other.cx_config_) {}

PauliExpPairBox::PauliExpPairBox() : PauliExpPairBox({}, 0., {}, 0.) {}

bool PauliExpPairBox::is_clifford() const {
  auto is_clifford0 = equiv_0(4 * t0_) || paulis0_.empty();
  auto is_clifford1 = equiv_0(4 * t1_) || paulis1_.empty();
  return is_clifford0 && is_clifford1;
}

SymSet PauliExpPairBox::free_symbols() const {
  return expr_free_symbols({t0_, t1_});
}

Op_ptr PauliExpPairBox::dagger() const {
  return std::make_shared<PauliExpPairBox>(
      paulis1_, -t1_, paulis0_, -t0_, cx_config_);
}

// Get the multiplicative change factor in Pauli angle during transpose ( -1 for
// odd nr of Y, 1 otherwise)
int transpose_angle_factor(const std::vector<Pauli> &paulis) {
  int pauli_odd_number_y = std::count_if(
                               paulis.begin(), paulis.end(),
                               [](auto pauli) { return pauli == Pauli::Y; }) %
                           2;
  // transform 0 -> 1, 1 -> -1
  return -pauli_odd_number_y * 2 + 1;
}

Op_ptr PauliExpPairBox::transpose() const {
  int pauli0_angle_factor = transpose_angle_factor(paulis0_);
  int pauli1_angle_factor = transpose_angle_factor(paulis1_);
  return std::make_shared<PauliExpPairBox>(
      paulis1_, pauli1_angle_factor * t1_, paulis0_, pauli0_angle_factor * t0_,
      cx_config_);
}

Op_ptr PauliExpPairBox::symbol_substitution(
    const SymEngine::map_basic_basic &sub_map) const {
  return std::make_shared<PauliExpPairBox>(
      this->paulis0_, this->t0_.subs(sub_map), this->paulis1_,
      this->t1_.subs(sub_map), this->cx_config_);
}

void PauliExpPairBox::generate_circuit() const {
  Circuit circ = Circuit(paulis0_.size());
  QubitPauliTensor pauli_tensor0(paulis0_);
  QubitPauliTensor pauli_tensor1(paulis1_);
  append_pauli_gadget_pair(
      circ, pauli_tensor0, t0_, pauli_tensor1, t1_, cx_config_);
  circ_ = std::make_shared<Circuit>(circ);
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
      paulis0, phase0, paulis1, phase1, j.at("cx_config").get<CXConfigType>());
  return set_box_id(
      box,
      boost::lexical_cast<boost::uuids::uuid>(j.at("id").get<std::string>()));
}

REGISTER_OPFACTORY(PauliExpPairBox, PauliExpPairBox)

PauliExpCommutingSetBox::PauliExpCommutingSetBox(
    const std::vector<std::pair<std::vector<Pauli>, Expr>> &pauli_gadgets,
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
  auto n_qubits = pauli_gadgets[0].first.size();
  for (const auto &gadget : pauli_gadgets) {
    if (gadget.first.size() != n_qubits) {
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
      [](const std::pair<std::vector<Pauli>, Expr> &pauli_exp) {
        return equiv_0(4 * pauli_exp.second) || pauli_exp.first.empty();
      });
}

bool PauliExpCommutingSetBox::paulis_commute() const {
  std::vector<QubitPauliString> pauli_strings;
  pauli_strings.reserve(pauli_gadgets_.size());
  for (const auto &pauli_gadget : pauli_gadgets_) {
    pauli_strings.emplace_back(pauli_gadget.first);
  }
  for (auto string0 = pauli_strings.begin(); string0 != pauli_strings.end();
       string0++) {
    for (auto string1 = string0 + 1; string1 != pauli_strings.end();
         string1++) {
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
    angles.push_back(pauli_exp.second);
  }
  return expr_free_symbols(angles);
}

Op_ptr PauliExpCommutingSetBox::dagger() const {
  std::vector<std::pair<std::vector<Pauli>, Expr>> dagger_gadgets;
  for (const auto &pauli_exp : pauli_gadgets_) {
    dagger_gadgets.emplace_back(pauli_exp.first, -pauli_exp.second);
  }
  return std::make_shared<PauliExpCommutingSetBox>(dagger_gadgets, cx_config_);
}

Op_ptr PauliExpCommutingSetBox::transpose() const {
  std::vector<std::pair<std::vector<Pauli>, Expr>> transpose_gadgets;
  for (const auto &pauli_exp : pauli_gadgets_) {
    int pauli_angle_factor = transpose_angle_factor(pauli_exp.first);
    transpose_gadgets.emplace_back(
        pauli_exp.first, pauli_angle_factor * pauli_exp.second);
  }
  return std::make_shared<PauliExpCommutingSetBox>(
      transpose_gadgets, cx_config_);
}

Op_ptr PauliExpCommutingSetBox::symbol_substitution(
    const SymEngine::map_basic_basic &sub_map) const {
  std::vector<std::pair<std::vector<Pauli>, Expr>> symbol_sub_gadgets;
  for (const auto &pauli_exp : pauli_gadgets_) {
    symbol_sub_gadgets.emplace_back(
        pauli_exp.first, pauli_exp.second.subs(sub_map));
  }
  return std::make_shared<PauliExpCommutingSetBox>(
      symbol_sub_gadgets, this->cx_config_);
}

void PauliExpCommutingSetBox::generate_circuit() const {
  Circuit circ = Circuit(pauli_gadgets_[0].first.size());

  std::list<std::pair<QubitPauliTensor, Expr>> gadgets;
  for (const auto &pauli_gadget : pauli_gadgets_) {
    gadgets.emplace_back(
        QubitPauliTensor(pauli_gadget.first), pauli_gadget.second);
  }
  std::set<Qubit> qubits;
  for (unsigned i = 0; i < pauli_gadgets_[0].first.size(); i++)
    qubits.insert(Qubit(i));

  Circuit cliff_circ = mutual_diagonalise(gadgets, qubits, cx_config_);
  circ.append(cliff_circ);

  Circuit phase_poly_circ = Circuit(pauli_gadgets_[0].first.size());

  for (const std::pair<QubitPauliTensor, Expr> &pgp : gadgets) {
    append_single_pauli_gadget(phase_poly_circ, pgp.first, pgp.second);
  }
  PhasePolyBox ppbox(phase_poly_circ);
  Circuit after_synth_circ = *ppbox.to_circuit();

  circ.append(after_synth_circ);
  circ.append(cliff_circ.dagger());

  circ_ = std::make_shared<Circuit>(circ);
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
      j.at("pauli_gadgets")
          .get<std::vector<std::pair<std::vector<Pauli>, Expr>>>(),
      j.at("cx_config").get<CXConfigType>());
  return set_box_id(
      box,
      boost::lexical_cast<boost::uuids::uuid>(j.at("id").get<std::string>()));
}

REGISTER_OPFACTORY(PauliExpCommutingSetBox, PauliExpCommutingSetBox)

}  // namespace tket
