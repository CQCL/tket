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

#include "UnitID.hpp"

#include <sstream>

#include "Json.hpp"

namespace tket {

std::string UnitID::repr() const {
  std::stringstream str;
  str << data_->name_;
  if (!data_->index_.empty()) {
    str << "[" << std::to_string(data_->index_[0]);
    for (unsigned i = 1; i < data_->index_.size(); i++) {
      str << ", " << std::to_string(data_->index_[i]);
    }
    str << "]";
  }
  return str.str();
}

void to_json(nlohmann::json& j, const Qubit& qb) { unitid_to_json(j, qb); }
void from_json(const nlohmann::json& j, Qubit& qb) { json_to_unitid(j, qb); }

void to_json(nlohmann::json& j, const Bit& cb) { unitid_to_json(j, cb); }
void from_json(const nlohmann::json& j, Bit& cb) { json_to_unitid(j, cb); }

void to_json(nlohmann::json& j, const Node& node) { unitid_to_json(j, node); }
void from_json(const nlohmann::json& j, Node& node) { json_to_unitid(j, node); }

void to_json(nlohmann::json& j, const qubit_map_t& qm) {
  for (const auto& pair : qm) {
    nlohmann::json qm_j;
    qm_j.push_back(pair.first);
    qm_j.push_back(pair.second);
    j.push_back(qm_j);
  }
}
void from_json(const nlohmann::json& j, qubit_map_t& qm) {
  for (const nlohmann::json& j_pair : j) {
    const auto& qb_pair = j_pair.get<std::pair<Qubit, Qubit>>();
    qm.insert(qb_pair);
  }
}

/* The functions below use the "construct on first use" idiom to return "global"
 * objects. They use pointers to ensure that the objects returned last
 * throughout static initialization and deinitialization. There is no memory
 * leak since the objects are only constructed once and the memory is reclaimed
 * on program termination.
 */

const std::string& q_default_reg() {
  static std::unique_ptr<const std::string> regname =
      std::make_unique<const std::string>("q");
  return *regname;
}

const std::string& c_default_reg() {
  static std::unique_ptr<const std::string> regname =
      std::make_unique<const std::string>("c");
  return *regname;
}

const std::string& node_default_reg() {
  static std::unique_ptr<const std::string> regname =
      std::make_unique<const std::string>("node");
  return *regname;
}

const std::string& c_debug_zero_prefix() {
  static std::unique_ptr<const std::string> regname =
      std::make_unique<const std::string>("tk_DEBUG_ZERO_REG");
  return *regname;
}

const std::string& c_debug_one_prefix() {
  static std::unique_ptr<const std::string> regname =
      std::make_unique<const std::string>("tk_DEBUG_ONE_REG");
  return *regname;
}

const std::string& c_debug_default_name() {
  static std::unique_ptr<const std::string> regname =
      std::make_unique<const std::string>("tket_assert");
  return *regname;
}

}  // namespace tket
