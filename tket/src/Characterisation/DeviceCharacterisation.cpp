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

#include "Characterisation/DeviceCharacterisation.hpp"

#include <optional>

namespace tket {

// simple get key wrapped in std::optional
template <typename T>
static std::optional<typename T::mapped_type> maybe_get(
    const T& error_map, const typename T::key_type& key) {
  const auto& it = error_map.find(key);
  return (it != error_map.end()) ? std::make_optional(it->second)
                                 : std::nullopt;
}

// single-qubit case
gate_error_t DeviceCharacterisation::get_error(const Node& n) const {
  std::optional<gate_error_t> maybe_err = maybe_get(default_node_errors_, n);
  return maybe_err ? *maybe_err : 0.;
}
gate_error_t DeviceCharacterisation::get_error(
    const Node& n, const OpType& op) const {
  std::optional<op_errors_t> maybe_dict = maybe_get(op_node_errors_, n);
  if (maybe_dict) {
    std::optional<gate_error_t> maybe_err = maybe_get(*maybe_dict, op);
    if (maybe_err) {
      return *maybe_err;
    }
  }
  return get_error(n);
}

// two-qubit case
gate_error_t DeviceCharacterisation::get_error(
    const Architecture::Connection& link) const {
  std::optional<gate_error_t> maybe_err = maybe_get(default_link_errors_, link);
  return maybe_err ? *maybe_err : 0.;
}
gate_error_t DeviceCharacterisation::get_error(
    const Architecture::Connection& link, const OpType& op) const {
  std::optional<op_errors_t> maybe_dict = maybe_get(op_link_errors_, link);
  if (maybe_dict) {
    std::optional<gate_error_t> maybe_err = maybe_get(*maybe_dict, op);
    if (maybe_err) {
      return *maybe_err;
    }
  }
  return get_error(link);
}

readout_error_t DeviceCharacterisation::get_readout_error(const Node& n) const {
  std::optional<gate_error_t> maybe_err = maybe_get(default_readout_errors_, n);
  return maybe_err ? *maybe_err : 0.;
}

bool DeviceCharacterisation::operator==(
    const DeviceCharacterisation& other) const {
  return (this->default_node_errors_ == other.default_node_errors_) &&
         (this->default_link_errors_ == other.default_link_errors_) &&
         (this->default_readout_errors_ == other.default_readout_errors_) &&
         (this->op_node_errors_ == other.op_node_errors_) &&
         (this->op_link_errors_ == other.op_link_errors_);
}

void to_json(nlohmann::json& j, const DeviceCharacterisation& dc) {
  j["def_node_errors"] = dc.default_node_errors_;
  j["def_link_errors"] = dc.default_link_errors_;
  j["readouts"] = dc.default_readout_errors_;
  j["op_node_errors"] = dc.op_node_errors_;
  j["op_link_errors"] = dc.op_link_errors_;
}

void from_json(const nlohmann::json& j, DeviceCharacterisation& dc) {
  dc.default_node_errors_ = j.at("def_node_errors").get<avg_node_errors_t>();
  dc.default_link_errors_ = j.at("def_link_errors").get<avg_link_errors_t>();
  dc.default_readout_errors_ = j.at("readouts").get<avg_readout_errors_t>();
  dc.op_node_errors_ = j.at("op_node_errors").get<op_node_errors_t>();
  dc.op_link_errors_ = j.at("op_link_errors").get<op_link_errors_t>();
}

}  // namespace tket
