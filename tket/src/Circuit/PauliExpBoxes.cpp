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

#include "Circuit/PauliExpBoxes.hpp"

#include <iostream>

#include "Circuit/CircUtils.hpp"
#include "Ops/OpJsonFactory.hpp"
#include "Converters/PauliGadget.hpp"

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
  // Negate the parameter if odd nr of paulis (mult with 1 if mod2=0, -1 if
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
    const unsigned n_qubits, const QubitPauliTensor pauli0, Expr t0,
    const QubitPauliTensor &pauli1, Expr t1, CXConfigType cx_config_type)
    : Box(OpType::PauliExpBox, op_signature_t(n_qubits, EdgeType::Quantum)),
      n_qubits_(n_qubits),
      pauli0_(pauli0),
      t0_(std::move(t0)),
      pauli1_(pauli1),
      t1_(std::move(t1)),
      cx_config_(cx_config_type) {}

PauliExpPairBox::PauliExpPairBox(const PauliExpPairBox &other)
    : Box(other),
      n_qubits_(other.n_qubits_),
      pauli0_(other.pauli0_),
      t0_(other.t0_),
      pauli1_(other.pauli1_),
      t1_(other.t1_),
      cx_config_(other.cx_config_) {}

PauliExpPairBox::PauliExpPairBox() : PauliExpPairBox(0, QubitPauliTensor(), 0., QubitPauliTensor(), 0.) {}

bool PauliExpPairBox::is_clifford() const {
  Expr angle0 = pauli0_.coeff == -1. ? -t0_ : t0_;
  Expr angle1 = pauli0_.coeff == -1. ? -t1_ : t1_;
  auto is_clifford0 = equiv_0(4 * angle0) || pauli0_.string.map.empty();
  auto is_clifford1 = equiv_0(4 * angle1) || pauli1_.string.map.empty();
  return is_clifford0 && is_clifford1;
}

SymSet PauliExpPairBox::free_symbols() const {
  Expr angle0 = pauli0_.coeff == -1. ? -t0_ : t0_;
  Expr angle1 = pauli0_.coeff == -1. ? -t1_ : t1_;
  return expr_free_symbols({t0_, t1_});
}

Op_ptr PauliExpPairBox::dagger() const {
  return std::make_shared<PauliExpPairBox>(
      n_qubits_, pauli1_, -t1_, pauli0_, -t0_, cx_config_);
}

// Get the multiplicative change factor in pauli angle during transpose ( -1 for
// odd nr of Y, 1 otherwise)
int transpose_angle_factor(
    const QubitPauliTensor &pauli) {
  int pauli_odd_number_y =
      std::count_if(
          pauli.string.map.begin(), pauli.string.map.end(),
          [](auto pauli_pair) { return pauli_pair.second == Pauli::Y; }) %
      2;
  // transform 0 -> 1, 1 -> -1
  return -pauli_odd_number_y * 2 + 1;
}

Op_ptr PauliExpPairBox::transpose() const {
  int pauli0_t_fac = transpose_angle_factor(pauli0_);
  int pauli1_t_fac = transpose_angle_factor(pauli1_);
  return std::make_shared<PauliExpPairBox>(
      n_qubits_, pauli1_, pauli1_t_fac * t1_, pauli0_, pauli0_t_fac * t0_,
      cx_config_);
}

Op_ptr PauliExpPairBox::symbol_substitution(
    const SymEngine::map_basic_basic &sub_map) const {
  return std::make_shared<PauliExpPairBox>(
      this->n_qubits_, this->pauli0_, this->t0_.subs(sub_map), this->pauli1_,
      this->t1_.subs(sub_map), this->cx_config_);
}

void PauliExpPairBox::generate_circuit() const {

  Circuit circ = Circuit();
  for (const auto& qubit_pauli: pauli0_.string.map ){
    circ.add_qubit(qubit_pauli.first);
  }
  for (const auto& qubit_pauli: pauli1_.string.map ){
    if (circ.contains_unit(qubit_pauli.first)) continue;
    circ.add_qubit(qubit_pauli.first);
  }
  append_pauli_gadget_pair(circ, pauli0_, t0_, pauli1_, t1_, cx_config_);
  circ_ = std::make_shared<Circuit>(circ);
}

nlohmann::json PauliExpPairBox::to_json(const Op_ptr &op) {
  const auto &box = static_cast<const PauliExpPairBox &>(*op);
  nlohmann::json j = core_box_json(box);
  j["n_qubits"] = box.get_n_qubits();
  j["paulis0"] = box.get_pauli0();
  j["phase0"] = box.get_phase0();
  j["paulis1"] = box.get_pauli1();
  j["phase1"] = box.get_phase1();
  j["cx_config"] = box.get_cx_config();
  return j;
}

Op_ptr PauliExpPairBox::from_json(const nlohmann::json &j) {
  PauliExpPairBox box = PauliExpPairBox(
      j.at("n_qubits").get<unsigned>(),
      j.at("pauli0").get<QubitPauliTensor>(),
      j.at("phase0").get<Expr>(),
      j.at("pauli1").get<QubitPauliTensor>(),
      j.at("phase1").get<Expr>(), j.at("cx_config").get<CXConfigType>());
  return set_box_id(
      box,
      boost::lexical_cast<boost::uuids::uuid>(j.at("id").get<std::string>()));
}

REGISTER_OPFACTORY(PauliExpPairBox, PauliExpPairBox)

}  // namespace tket
