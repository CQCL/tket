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
#include <vector>

#include "OpType/OpType.hpp"

namespace tket {
namespace internal {

// This contains data about number of parameters and qubits for each op type
// representing a valid gate, so that we can just go through all of them
// automatically.
struct GatesData {
  // KEY: number of qubits
  // VALUE: another map: (number of parameters -> list of gates)
  typedef std::map<unsigned, std::map<unsigned, std::vector<OpType>>> InputData;

  InputData input_data;

  // KEY: an op which can take a variable number of qubits
  // VALUE: the minimum number of qubits which must be supplied
  std::map<OpType, unsigned> min_number_of_qubits_for_variable_qubit_type;

  static const GatesData& get();
};

}  // namespace internal
}  // namespace tket
