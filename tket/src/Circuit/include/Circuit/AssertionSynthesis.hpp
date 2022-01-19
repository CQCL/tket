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

#include "Circuit.hpp"
#include "Utils/EigenConfig.hpp"
#include "Utils/MatrixAnalysis.hpp"

namespace tket {

/**
 * Synthesise an assertion circuit from an arbitrary 2x2, 4x4 or 8x8 projector
 * matrix.
 *
 * @param P projector matrix in \ref BasisOrder::ilo
 *
 * @return circuit implementing the projector based assertion and the expected
 * readouts
 */
std::tuple<Circuit, std::vector<bool>> projector_assertion_synthesis(
    const Eigen::MatrixXcd &P);

/**
 * Synthesise an assertion circuit from a list of Paulis strings with +/-1
 * coefficients.
 *
 * @param paulis list of Paulis strings
 *
 * @return circuit implementing the stabiliser based assertion and the expected
 * readouts
 */
std::tuple<Circuit, std::vector<bool>> stabiliser_assertion_synthesis(
    const PauliStabiliserList &paulis);
}  // namespace tket
