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
namespace tket_sim {

// NOTE: this is an identical copy of a file in tket-tests, in the Simulation
// folder. This is deliberate! Please keep both in sync!

enum class MatrixEquivalence { EQUAL, EQUAL_UP_TO_GLOBAL_PHASE };

/** Compare EITHER two state vectors, OR two unitary matrices,
 *  calculated from two circuits. Automatically detects which.
 *  Returns TRUE if the circuits appear to be equivalent
 *  (EITHER with equal unitaries, OR only up to global phase).
 *  Also checks that statevectors have norm 1,
 *  and unitaries really are almost unitary. Throws if not.
 *
 * @param m1 First matrix obtained from circuit (state vector or unitary).
 * @param m2 Second matrix obtained from circuit (state vector or unitary).
 * @param equivalence Whether we demand exact equality, or only compare
 *                              up to global phase.
 * @param tolerance The numerical tolerance for approximate equivalence.
 * @return Whether the two circuits which gave these two matrices
 *                      appear to be equivalent.
 */
bool compare_statevectors_or_unitaries(
    const Eigen::MatrixXcd& m1, const Eigen::MatrixXcd& m2,
    MatrixEquivalence equivalence = MatrixEquivalence::EQUAL,
    double tolerance = 1e-10);

}  // namespace tket_sim
}  // namespace tket
