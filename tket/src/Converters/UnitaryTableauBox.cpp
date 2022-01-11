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

#include "UnitaryTableauBox.hpp"

#include "Ops/OpJsonFactory.hpp"

namespace tket {

UnitaryTableauBox::UnitaryTableauBox(const UnitaryTableau& tab)
    : Box(OpType::UnitaryTableauBox), tab_(tab) {
  std::set<Qubit> tab_qbs = tab_.get_qubits();
  std::set<Qubit> qbs;
  for (unsigned i = 0; i < tab_qbs.size(); ++i) {
    qbs.insert(Qubit(i));
  }
  if (tab_qbs != qbs)
    throw NotValid(
        "UnitaryTableauBox requires tableau qubits to have default, "
        "consecutive indexing");
}

UnitaryTableauBox::UnitaryTableauBox(
    const MatrixXb& xx, const MatrixXb& xz, const VectorXb& xph,
    const MatrixXb& zx, const MatrixXb& zz, const VectorXb& zph)
    : Box(OpType::UnitaryTableauBox), tab_(xx, xz, xph, zx, zz, zph) {}

Op_ptr UnitaryTableauBox::dagger() const {
  return std::make_shared<const UnitaryTableauBox>(tab_.dagger());
}

Op_ptr UnitaryTableauBox::transpose() const {
  return std::make_shared<const UnitaryTableauBox>(tab_.transpose());
}

Op_ptr UnitaryTableauBox::symbol_substitution(
    const SymEngine::map_basic_basic&) const {
  return Op_ptr();
}

SymSet UnitaryTableauBox::free_symbols() const { return SymSet(); }

bool UnitaryTableauBox::is_equal(const Op& op_other) const {
  const UnitaryTableauBox& other =
      dynamic_cast<const UnitaryTableauBox&>(op_other);
  return this->tab_ == other.tab_;
}

const UnitaryTableau& UnitaryTableauBox::get_tableau() const { return tab_; }

nlohmann::json UnitaryTableauBox::to_json(const Op_ptr& op) {
  const auto& box = static_cast<const UnitaryTableauBox&>(*op);
  nlohmann::json j = core_box_json(box);
  j["tab"] = box.get_tableau();
  return j;
}

Op_ptr UnitaryTableauBox::from_json(const nlohmann::json& j) {
  UnitaryTableau tab(0);
  j.at("tab").get_to(tab);
  return std::make_shared<const UnitaryTableauBox>(tab);
}

REGISTER_OPFACTORY(UnitaryTableauBox, UnitaryTableauBox)

op_signature_t UnitaryTableauBox::get_signature() const {
  return op_signature_t(tab_.get_qubits().size(), EdgeType::Quantum);
}

void UnitaryTableauBox::generate_circuit() const {
  circ_ = std::make_shared<Circuit>(unitary_tableau_to_circuit(tab_));
}

}  // namespace tket
