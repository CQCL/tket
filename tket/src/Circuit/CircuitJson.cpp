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

#include "Circuit.hpp"
#include "Utils/Json.hpp"

namespace tket {

void to_json(nlohmann::json& j, const Circuit& circ) {
  const auto name = circ.get_name();
  if (name) {
    j["name"] = name.value();
  }
  j["phase"] = circ.get_phase();
  j["qubits"] = circ.all_qubits();
  j["bits"] = circ.all_bits();
  const auto impl = circ.implicit_qubit_permutation();
  // empty maps are mapped to null instead of empty array
  if (impl.empty()) {
    j["implicit_permutation"] = nlohmann::json::array();
  } else {
    j["implicit_permutation"] = impl;
  }
  j["commands"] = nlohmann::json::array();
  for (const Command& com : circ) {
    j["commands"].push_back(com);
  }
}

void from_json(const nlohmann::json& j, Circuit& circ) {
  circ = Circuit();

  if (j.contains("name")) {
    circ.set_name(j["name"].get<std::string>());
  }
  circ.add_phase(j.at("phase").get<Expr>());
  const auto& qubits = j.at("qubits").get<qubit_vector_t>();
  for (const auto& qb : qubits) {
    circ.add_qubit(qb);
  }
  const auto& bits = j.at("bits").get<bit_vector_t>();
  for (const auto& b : bits) {
    circ.add_bit(b);
  }

  for (const auto& j_com : j.at("commands")) {
    const auto& com = j_com.get<Command>();
    circ.add_op(com.get_op_ptr(), com.get_args(), com.get_opgroup());
  }
  const auto& imp_perm = j.at("implicit_permutation").get<qubit_map_t>();
  circ.permute_boundary_output(imp_perm);
}

}  // namespace tket
