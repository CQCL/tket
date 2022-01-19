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

#include "Circuit/Circuit.hpp"
#include "Utils/UnitID.hpp"

namespace tket {

/**
 * Truncate the circuit so that the last operation on every (non-empty)
 * classical wire is a Measure, and return a classical postprocessing circuit as
 * well as the truncated circuit.
 *
 * The sequential composition of the two circuits returned is equivalent to the
 * original circuit.
 *
 * @param circ circuit to separate
 *
 * @return a pair of circuits (C0, C1) where C0 has no operations following the
 * final Measure on every classical wire; C1 consists of purely classical
 * operations on the same classical units; and C0 followed by C1 is equivalent
 * to the original circuit.
 */
std::pair<Circuit, Circuit> separate_classical(const Circuit &circ);

}  // namespace tket
