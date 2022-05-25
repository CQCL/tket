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
#include "DAGDefs.hpp"
#include "Gate/GatePtr.hpp"
#include "Utils/EigenConfig.hpp"
#include "Utils/PauliStrings.hpp"

namespace tket {

/**
 * Convert a vertex holding a TK1 operation to its corresponding matrix
 *
 * @param circ circuit
 * @param vert vertex in circuit
 *
 * @pre \p vert is an operation of \ref OpType::TK1
 * @pre operation has no symbolic parameters
 *
 * @return corresponding unitary matrix
 */
Eigen::Matrix2cd get_matrix(const Circuit& circ, const Vertex& vert);

/**
 * Convert a one-qubit circuit of TK1 operations to its corresponding matrix
 *
 * @param circ circuit
 *
 * @pre all vertices are operations of \ref OpType::TK1
 * @pre circuit has no symbolic parameters
 *
 * @return corresponding unitary matrix
 */
Eigen::Matrix2cd get_matrix_from_circ(const Circuit& circ);

/**
 * @brief Convert a two-qubit circuit to its corresponding matrix
 *
 * @param circ circuit
 *
 * @pre \p circ is composed of CX, TK2 and single-qubit gates only
 * @pre circuit has no symbolic parameters
 * @post matrix is in \ref BasisOrder::ilo
 *
 * @return corresponding unitary matrix
 */
Eigen::Matrix4cd get_matrix_from_2qb_circ(const Circuit& circ);

/**
 * Convert a 4x4 unitary matrix optimally to a corresponding circuit
 *
 * This will express `U` using a single TK2 gate. See `decompose_TK2` to
 * decompose this gate further into other primitives.
 *
 * @param U unitary matrix
 *
 * @pre \p U is in \ref BasisOrder::ilo
 * @post circuit consists of one TK2 and single-qubit gates
 *
 * @return circuit implementing U with minimal error
 */
Circuit two_qubit_canonical(const Eigen::Matrix4cd& U);

/**
 * Decompose a unitary matrix into a 2-CX circuit following a diagonal operator.
 *
 * Given an arbitrary unitary 4x4 matrix, this method returns a circuit C and a
 * complex number z such that |z|=1 and U = VD where V is the unitary
 * implemented by the circuit and D = diag(z, z*, z*, z). The circuit C consists
 * of CX and TK1 gates only and has at most 2 CX gates.
 *
 * @param U unitary matrix
 *
 * @return circuit and parameter defining diagonal operator
 */
std::pair<Circuit, Complex> decompose_2cx_VD(const Eigen::Matrix4cd& U);

/**
 * Decompose a unitary matrix into a 2-CXs preceding a diagonal operator.
 *
 * Given an arbitrary unitary 4x4 matrix, this method returns a circuit C and a
 * complex number z such that |z|=1 and U = DV where V is the unitary
 * implemented by the circuit and D = diag(z, z*, z*, z). The circuit C consists
 * of CX and TK1 gates only and has at most 2 CX gates.
 *
 * @param U unitary matrix
 *
 * @return circuit and parameter defining diagonal operator
 */
std::pair<Circuit, Complex> decompose_2cx_DV(const Eigen::Matrix4cd& U);

/**
 * Construct a phase gadget
 *
 * @param n_qubits number of qubits
 * @param t phase parameter
 * @param cx_config CX configuration
 *
 * @return phase gadget
 */
Circuit phase_gadget(
    unsigned n_qubits, const Expr& t,
    CXConfigType cx_config = CXConfigType::Snake);

/**
 * Construct a 'Pauli gadget' corresponding to a tensor of Pauli operators.
 *
 * The returned circuit implements the unitary operator
 * \f$ e^{-\frac12 i \pi t \sigma_0 \otimes \sigma_1 \otimes \cdots} \f$
 * where \f$ \sigma_i \in \{I,X,Y,Z\} \f$ are the Pauli operators.
 *
 * @param paulis Pauli operators
 * @param t angle in half-turns
 * @param cx_config CX configuration
 */
Circuit pauli_gadget(
    const std::vector<Pauli>& paulis, const Expr& t,
    CXConfigType cx_config = CXConfigType::Snake);

/**
 * Utility function to replace all CX gates with TK2 and single-qubit gates.
 *
 * @param c circuit to modify
 */
void replace_CX_with_TK2(Circuit& c);

/**
 * Express a gate as a circuit using TK2 as the only multi-qubit gate.
 *
 * @param op operation
 *
 * @return circuit representing the operation
 */
Circuit with_TK2(Gate_ptr op);

/**
 * Express a gate as a circuit using CX as the only multi-qubit gate.
 *
 * @param op operation
 *
 * @return circuit representing the operation
 */
Circuit with_CX(Gate_ptr op);

/**
 * Construct a controlled version of a given circuit.
 *
 * @param c circuit
 * @param n_controls number of controls
 *
 * @return controlled circuit
 *
 * @pre \p c consists only of unitary \ref Gate operations
 * @pre \p c has a single default qubit register
 * @pre \p c has no implicit wireswaps
 * @post the returned circuit has a single default qubit register, whose initial
 *  qubits act as the controls over the other qubits, which correspond to those
 *  of \p c in the same order
 */
Circuit with_controls(const Circuit& c, unsigned n_controls = 1);

}  // namespace tket
