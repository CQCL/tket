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

#include "OpType/OpType.hpp"
#include "Utils/UnitID.hpp"

namespace tket {

// errors are double precision values given by 1. - fidelity (0. error is
// perfect gate)
typedef double gate_error_t;
typedef double readout_error_t;

// default errors per Node (average error)
typedef std::map<Node, gate_error_t> avg_node_errors_t;
typedef std::map<Node, readout_error_t> avg_readout_errors_t;
typedef std::map<std::pair<Node, Node>, gate_error_t> avg_link_errors_t;

// OpType-specific errors per Node
typedef std::map<OpType, gate_error_t> op_errors_t;
typedef std::map<Node, op_errors_t> op_node_errors_t;
typedef std::map<std::pair<Node, Node>, op_errors_t> op_link_errors_t;

}  // namespace tket
