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

#include "OpDesc.hpp"

#include <algorithm>

#include "OpTypeFunctions.hpp"

namespace tket {

OpDesc::OpDesc(OpType type)
    : type_(type),
      info_(optypeinfo().at(type)),
      is_meta_(is_metaop_type(type)),
      is_box_(is_box_type(type)),
      is_gate_(is_gate_type(type)),
      is_flowop_(is_flowop_type(type)),
      is_classical_(is_classical_type(type)),
      is_rotation_(is_rotation_type(type)),
      is_oneway_(is_oneway_type(type)),
      is_clifford_(is_clifford_type(type)),
      is_parameterised_pauli_rotation_(
          is_parameterised_pauli_rotation_type(type)) {}

OpType OpDesc::type() const { return type_; }

std::string OpDesc::name() const { return info_.name; }

std::string OpDesc::latex() const { return info_.latex_name; }

unsigned OpDesc::n_params() const { return info_.n_params(); }

std::optional<op_signature_t> OpDesc::signature() const {
  return info_.signature;
}

OptUInt OpDesc::n_qubits() const {
  if (info_.signature) {
    unsigned n = std::count(
        info_.signature->begin(), info_.signature->end(), EdgeType::Quantum);
    return n;
  }
  return any;
}

OptUInt OpDesc::n_boolean() const {
  if (info_.signature) {
    unsigned n = std::count(
        info_.signature->begin(), info_.signature->end(), EdgeType::Boolean);
    return n;
  }
  return any;
}

OptUInt OpDesc::n_classical() const {
  if (info_.signature) {
    unsigned n = std::count(
        info_.signature->begin(), info_.signature->end(), EdgeType::Classical);
    return n;
  }
  return any;
}

bool OpDesc::is_meta() const { return is_meta_; }

bool OpDesc::is_box() const { return is_box_; }

bool OpDesc::is_gate() const { return is_gate_; }

bool OpDesc::is_flowop() const { return is_flowop_; }

bool OpDesc::is_classical() const { return is_classical_; }

bool OpDesc::is_rotation() const { return is_rotation_; }

unsigned OpDesc::param_mod(unsigned i) const { return info_.param_mod[i]; }

bool OpDesc::is_oneway() const { return is_oneway_; }

bool OpDesc::is_singleq_unitary() const {
  return n_qubits() && n_qubits().value() == 1 && !is_oneway();
}

bool OpDesc::is_clifford_gate() const { return is_clifford_; }

bool OpDesc::is_parameterised_pauli_rotation() const {
  return is_parameterised_pauli_rotation_;
}

}  // namespace tket
