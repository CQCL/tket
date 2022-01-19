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
#include "Utils/PauliStrings.hpp"

namespace tket {

/**
 * Captures rules for conjugating a pauli-gadget with single-qubit Clifford
 * gates Maps gate and pauli to the new pauli after the conjugation and whether
 * or not a phase-flip is induced
 */
std::pair<Pauli, bool> conjugate_Pauli(
    OpType op, Pauli p, bool reverse = false);

/**
 * Methods to conjugate a QubitPauliTensor with Clifford gates to
 * change basis
 * Transforms P to P' such that
 * reverse = false : --P'-- = --op--P--opdg--
 * reverse = true  : --P'-- = --opdg--P--op--
 */
void conjugate_PauliTensor(
    QubitPauliTensor &qpt, OpType op, const Qubit &q, bool reverse = false);
void conjugate_PauliTensor(
    QubitPauliTensor &qpt, OpType op, const Qubit &q0, const Qubit &q1);
void conjugate_PauliTensor(
    QubitPauliTensor &qpt, OpType op, const Qubit &q0, const Qubit &q1,
    const Qubit &qb2);

}  // namespace tket
