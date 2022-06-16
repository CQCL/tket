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

#include "Characterisation/ErrorTypes.hpp"
#include "Transform.hpp"

namespace tket {

namespace Transforms {

// removes gate-inverse pairs, merges rotations, removes identity rotations,
// and removes redundant gates before measure Expects: Any gates Produces: The
// same gate set
Transform remove_redundancies();

/**
 * Squash all single-qubit gates to TK1.
 */
Transform squash_1qb_to_tk1();

// moves single qubit operations past multiqubit operations they commute with,
// towards front of circuit (hardcoded)
// Expects: Any gates
// Produces: Any gates
Transform commute_through_multis();

// commutes Rz gates through ZZMax, and combines adjacent ZZMax gates
// Expects: ZZMax, Rz, Rx
// Produces: ZZMax, Rz, Rx
Transform commute_and_combine_HQS2();

// squashes each sequence of two-qubit instructions into their canonical 3-CX
// form (Cartan/KAK decomposition) The optional CX gate fidelity cx_fidelity
// is used to produce circuit decompositions that maximise expected circuit
// fidelity, assuming perfect local operations (i.e. with fidelity 1.).
// Expects: CX and any single qubit gates
// Produces: CX and single TK1 qubit gates
Transform two_qubit_squash(double cx_fidelity = 1.);

// 1qb squashing into -Rz-Rx-Rz- or -Rx-Rz-Rx- form
// Expects: Rx, Rz, and any multi-qubit gates
// Produces: Rx, Rz, and any multi-qubit gates
Transform reduce_XZ_chains();

/**
 * general u_squash by converting any chains of p, q gates (p, q in
 * {Rx,Ry,Rz}) to triples -p-q-p- or -q-p-q-
 *
 * if strict = false (default), then the chains will be squashed to either
 * -p-q-p- or -q-p-q-, and the third rotation will be commuted through the
 * next gate whenever possible if strict = true, then all chains will be
 * squashed to -p-q-p-
 *
 * Expects: p, q, and any multi-qubit gates
 * Produces: p, q, and any multi-qubit gates
 */
Transform squash_1qb_to_pqp(
    const OpType& q, const OpType& p, bool strict = false);

// identifies single-qubit chains and squashes them in the target gate set
// Expects: any gates
// Produces: singleqs and any multi-qubit gates
Transform squash_factory(
    const OpTypeSet& singleqs,
    const std::function<Circuit(const Expr&, const Expr&, const Expr&)>&
        tk1_replacement);

// commutes single qubit gates through SWAP gates, leaving them on the
// PhysicalQubit with best fidelity for the given OP Expects: any single qubit
// gates, SWAP gates Produces: any single qubit gates, SWAP gates
Transform commute_SQ_gates_through_SWAPS(const avg_node_errors_t& node_errors);
Transform commute_SQ_gates_through_SWAPS(const op_node_errors_t& node_errors);

/**
 * @brief Absorb Rz gates into NPhasedX where possible
 *
 * For any gate NPhasedX(𝛼, 𝛽) in the circuit, try to change 𝛽 to reduce the
 * number of surrounding Rz gates.
 */
Transform absorb_Rz_NPhasedX();

}  // namespace Transforms

}  // namespace tket
