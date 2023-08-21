// Copyright 2019-2023 Cambridge Quantum Computing
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

#include "DiagUtils.hpp"
#include "tket/Circuit/Boxes.hpp"
#include "tket/Circuit/Circuit.hpp"

namespace tket {

/**
 * Check whether there are any qubits which only requires a single
 * qubit Clifford to make all Paulis I or Z
 */
void check_easy_diagonalise(
    std::list<std::pair<QubitPauliTensor, Expr>> &gadgets,
    std::set<Qubit> &qubits, Circuit &circ);

/**
 * Given two qubits, attempt to find a basis in which a single CX will
 * make the Paulis on one of qubits fully diagonal
 */
std::optional<std::pair<Pauli, Pauli>> check_pair_compatibility(
    const Qubit &qb1, const Qubit &qb2,
    const std::list<std::pair<QubitPauliTensor, Expr>> &gadgets);

/**
 * Diagonalise a qubit greedily by finding the Pauli Gadget with
 * the lowest residual support over the non-diagonal qubits, and apply
 * single qubit Cliffords and CXs to make it a `ZIII...I` string
 */
void greedy_diagonalise(
    const std::list<std::pair<QubitPauliTensor, Expr>> &gadgets,
    std::set<Qubit> &qubits, Conjugations &conjugations, Circuit &circ,
    CXConfigType cx_config);

/**
 * Diagonalise a mutually commuting set of Pauli strings. Modifies the
 * list of Pauli strings in place, and returns the Clifford circuit
 * required to generate the initial set.
 */
Circuit mutual_diagonalise(
    std::list<std::pair<QubitPauliTensor, Expr>> &gadgets,
    std::set<Qubit> qubits, CXConfigType cx_config);
/**
 * Applies Clifford conjugations to a QubitPauliTensor
 */
void apply_conjugations(
    QubitPauliTensor &qps, const Conjugations &conjugations);

/**
 * Given a Pauli tensor P, produces a short Clifford circuit C which maps P to Z
 * on a single qubit, i.e. Z_i C P = C. This can be viewed as the components
 * required to synthesise a single Pauli gadget C^dag RZ(a)_i C = exp(-i pi a
 * P/2) (up to global phase), or as a diagonalisation of a single Pauli string
 * along with CXs to reduce it to a single qubit. Returns the circuit C and the
 * qubit i where the Z ends up.
 */
std::pair<Circuit, Qubit> reduce_pauli_to_z(
    const QubitPauliTensor &pauli, CXConfigType cx_config);

/**
 * Given a pair of anticommuting Pauli tensors P0, P1, produces a short Clifford
 * circuit C which maps P0 to Z and P1 to X on the same qubit, i.e. Z_i C P0 = C
 * = X_i C P1. This can be viewed as the components required to synthesise a
 * pair of noncommuting Pauli gadgets C^dag RX(b)_i RZ(a)_i C = exp(-i pi b
 * P1/2) exp(-i pi a P0/2) (up to global phase). This is not strictly a
 * diagonalisation because anticommuting strings cannot be simultaneously
 * diagonalised. Returns the circuit C and the qubit i where the Z and X end up.
 */
std::pair<Circuit, Qubit> reduce_anticommuting_paulis_to_z_x(
    QubitPauliTensor pauli0, QubitPauliTensor pauli1, CXConfigType cx_config);

/**
 * Given a pair of commuting Pauli tensors P0, P1, produces a short Clifford
 * circuit C which maps P0 and P1 to Z on different qubits, i.e. Z_i C P0 = C =
 * Z_j C P1. This can be viewed as the components required to synthesise a pair
 * of commuting Pauli gadgets C^dag RZ(b)_j RZ(a)_i C = exp(-i pi b P1/2) exp(-i
 * pi a P0/2) (up to global phase), or as a mutual diagonalisation of two Pauli
 * strings along with CXs to reduce them to independent, individual qubits.
 * Returns the circuit C and the qubits i and j where the Zs end up.
 */
std::tuple<Circuit, Qubit, Qubit> reduce_commuting_paulis_to_zi_iz(
    QubitPauliTensor pauli0, QubitPauliTensor pauli1, CXConfigType cx_config);

}  // namespace tket
