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

#pragma once

#include <map>
#include <nlohmann/json.hpp>

#include "Architecture/Architecture.hpp"
#include "ErrorTypes.hpp"
#include "OpType/OpType.hpp"
#include "Utils/UnitID.hpp"

/**
 * Defines tket::DeviceCharacterisation, used in NoiseAwarePlacement and in
 * commute_SQ_gates_through_SWAPS as a simple device noise model.
 * This is just a container of errors.
 *
 * This supports single-qubit errors, two-qubit errors and readout errors.
 * Errors can either be OpType-specific, or a default value (average over all
 * possible OpTypes) If an OpType-specific value is provided, this will be used.
 * If not it will fallback to the default value for the given Node or Node pair,
 * which itself falls back to zero error.
 */

namespace tket {

class DeviceCharacterisation {
 public:
  DeviceCharacterisation(
      avg_node_errors_t _node_errors = {}, avg_link_errors_t _link_errors = {},
      avg_readout_errors_t _readout_errors = {})
      : default_node_errors_(_node_errors),
        default_link_errors_(_link_errors),
        default_readout_errors_(_readout_errors),
        op_node_errors_(),
        op_link_errors_() {}
  explicit DeviceCharacterisation(
      op_node_errors_t _node_errors, op_link_errors_t _link_errors = {},
      avg_readout_errors_t _readout_errors = {})
      : default_node_errors_(),
        default_link_errors_(),
        default_readout_errors_(_readout_errors),
        op_node_errors_(_node_errors),
        op_link_errors_(_link_errors) {}

  // get device gate errors, preferring OpType-specific over default values over
  // 0. error
  // single-qubit case
  gate_error_t get_error(const Node& n) const;
  gate_error_t get_error(const Node& n, const OpType& op) const;
  // two-qubit case
  gate_error_t get_error(const Architecture::Connection& link) const;
  gate_error_t get_error(
      const Architecture::Connection& link, const OpType& op) const;
  // readout errors
  readout_error_t get_readout_error(const Node& n) const;

  bool operator==(const DeviceCharacterisation& other) const;

  friend void to_json(nlohmann::json& j, const DeviceCharacterisation& dc);
  friend void from_json(const nlohmann::json& j, DeviceCharacterisation& dc);

 private:
  // default errors per Node
  avg_node_errors_t default_node_errors_;
  avg_link_errors_t default_link_errors_;
  avg_readout_errors_t default_readout_errors_;

  // OpType-specific errors per Node
  op_node_errors_t op_node_errors_;
  op_link_errors_t op_link_errors_;
};

JSON_DECL(DeviceCharacterisation)

}  // namespace tket
