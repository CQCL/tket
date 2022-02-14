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

#include "Utils/MatrixAnalysis.hpp"

namespace tket {
class Circuit;
namespace tket_sim {
namespace internal {
class GateNodesBuffer;

/** Break up the circuit into individual gates and boxes,
 *  and pass the data for each component one-by-one into the buffer object.
 *  The buffer object is responsible for processing the data
 *  to obtain full (2^n)*(2^n) unitaries and multiplying them as appropriate.
 */
void decompose_circuit(
    const Circuit& circ, GateNodesBuffer& buffer, double abs_epsilon);

}  // namespace internal
}  // namespace tket_sim
}  // namespace tket
