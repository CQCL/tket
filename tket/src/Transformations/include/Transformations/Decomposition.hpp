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

#include "Architecture/Architecture.hpp"
#include "Transform.hpp"

namespace tket {

namespace Transforms {

/**
 * Decomposes all multi-qubit unitary gates into TK2 and single-qubit gates.
 *
 * Ignores boxes.
 *
 * Expects: any gates
 * Produces: TK2 and any single-qubit gates
 */
Transform decompose_multi_qubits_TK2();

/**
 * Decomposes all multi-qubit unitary gates into CX and single-qubit gates.
 *
 * Ignores boxes.
 *
 * Expects: any gates
 * Produces: CX and any single-qubit gates
 */
Transform decompose_multi_qubits_CX();

/**
 * Decomposes all single-qubit unitary gates into TK1 gates. Ignores boxes.
 */
Transform decompose_single_qubits_TK1();

/**
 * Starting with Rz, Ry and multi-qubit gates, replace all singles with TK1.
 */
Transform decompose_ZYZ_to_TK1();

// converts all single-qubit gates into Rz and Rx gates
// Expects: any gates
// Produces: Rz, Rx and any multi-qubit gates
Transform decompose_ZX();

// converts all single-qubit gates into Rz and Ry gates
// Expects: any gates
// Produces: Rz, Ry and any multi-qubit gates
Transform decompose_ZY();

// converts all single-qubit gates into Rx and Ry gates
// Expects: any gates
// Produces: Rx, Ry and any multi-qubit gates
Transform decompose_XY();

// converts all TK1-gates into rzrx gates
// Expects: TK1-gates and any multi-qubit gates
// Produces: Rz, Rx and any multi-qubit gates
Transform decompose_tk1_to_rzrx();

// replaces CXs with ECR, SX, Rz
// Expects: CX and any single-qubit gates
// Produces: ECR and any single-qubit gates
Transform decompose_CX_to_ECR();

// replaces CXs with ZZMax
// Expects: CX and any single-qubit gates
// Produces: ZZMax and any single-qubit gates
Transform decompose_CX_to_HQS2();

// replaces Rz-Rx-Rz triples with Rz-PhasedX pairs
// Expects: Rz, Rx, and any multi-qubit gates
// Produces: Rz, PhasedX, and any multi-qubit gates
Transform decompose_ZX_to_HQS1();

// converts all CX gates into Molmer-Sorensen gates by recognising exp(-i XX *
// angle * pi/2) and converting rest to exp(-i XX * pi/4) Expects: CX, Rx, and
// other single-qubit gates Produces: Molmer-Sorensen, U2, and other
// single-qubit gates
Transform decompose_MolmerSorensen();

/**
 * @brief A simple struct to store some two-qubit gate fidelities.
 *
 * We currently support CX, ZZMax and ZZPhase.
 *
 * This is meant to be easily extensible when further gate types will be
 * supported.
 *
 */
struct TwoQbFidelities {
  std::optional<double> CX_fidelity;
  std::optional<double> ZZMax_fidelity;
  std::optional<std::function<double(double)>> ZZPhase_fidelity;
};

/**
 * @brief Decomposes each TK2 gate into two-qubit gates.
 *
 * We currently support CX, ZZMax and ZZPhase.
 *
 * Decompose each TK2 gate into two-qubit gates in a noise-aware way. Supported
 * two-qubit gate fidelities will be used to return the optimal decomposition of
 * each TK2 gate, taking noise into consideration.
 *
 * If no fidelities are provided, the decomposition will be exact, using CX
 * gates.
 *
 * If the TK2 angles are symbolic values, the decomposition will be exact
 * (i.e. not noise-aware). It is not possible in general to obtain optimal
 * decompositions for arbitrary symbolic parameters, so consider substituting
 * for concrete values if possible.
 *
 * @param fid The two-qubit gate fidelities (optional).
 * @return Transform
 */
Transform decompose_TK2(const TwoQbFidelities& fid);
Transform decompose_TK2();

/**
 * @brief Synthesise ZZPhase gates from CX and Rz, as well as XX/YYPhase.
 *
 * Expects: CX, Rz, XXPhase, YYPhase Produces: ZZPhase, CX, Rz.
 */
Transform decompose_ZZPhase();

/**
 * Decompose single-qubit Clifford gates to a standard Cliffford gate set.
 *
 * Replaces all single-qubit unitary gates outside the set {Z, X, S, V} that are
 * recognized as Clifford operations with an equivalent sequence of gates from
 * that set.
 */
Transform decompose_cliffords_std();

Transform decompose_ZX_to_cliffords();

// converts all valid CX-Rz-CX strings and CX-Rx-CX strings to 2qb
// PhaseGadgets Expects: CX, Rz, Rx and any other single-qubit gates Produces:
// PhaseGadgets, and any other gates
Transform decompose_PhaseGadgets();

/**
 * Replaces all boxes by their decomposition using Box::to_circuit
 * Expects: any gateset
 * returns potentially all gates
 */
Transform decomp_boxes();

/**
 * Replaces all CX+Rz sub circuits by PhasePolyBox
 * Expects: only CX + Rz + H (and measure + reset + collapse + barrier)
 * @param min_size value for the minimal number of CX in each box, groups with
 * less than min_size CX gates are not converted to a PhasePolyBox, dafault
 * value is 0
 * @return Transformation to perform the conversion
 */
Transform compose_phase_poly_boxes(unsigned min_size = 0);

// converts all SWAP gates to given replacement circuit (not checked to
// preserve unitary) Expects: SWAP gates, replacement circuit Produces:
// Instances of the replacement circuit,
Transform decompose_SWAP(const Circuit& replacement_circuit);

// converts all SWAP gates to 3 CX gates
// providing an Architecture will prefer an orientation that reduces H
// redirection cost -X-    -C-X-C-   -X-C-X-
//  |   =  | | |  =  | | |
// -X-    -X-C-X-   -C-X-C-
// Expects: SWAP gates
// Produces: CX gates
Transform decompose_SWAP_to_CX(const Architecture& arc = Architecture());

// converts all BRIDGE (distance 2 CX) gates to 4 CX gates
// -B-   -C-   -C---C---   ---C---C-
//  R     |     |   |         |   |
// -I- = --- = -X-C-X-C- = -C-X-C-X-
//  D     |       |   |     |   |
// -G-   -X-   ---X---X-   -X---X---
//  E
// Expects: BRIDGE gates
// Produces: CX gates
Transform decompose_BRIDGE_to_CX();

// converts -C- gates to -H-X-H- depending on direction of Architecture edges
//           |              |
//          -X-          -H-C-H-
// Expects: CX gates
// Produces CX and H gates
Transform decompose_CX_directed(const Architecture& arc);

/**
 * @brief Decompose NPhasedX gates into single-qubit PhasedX gates.
 */
Transform decompose_NPhasedX();

/**
 * @brief Turns all PhasedX and NPhasedX gates into global gates
 *
 * Replaces any PhasedX gates with global NPhasedX gates. By default, this
 * transform will squash all single-qubit gates to PhasedX and Rz gates before
 * proceeding further. Existing non-global NPhasedX will not be preserved. This
 * is the recommended setting for best performance. If squashing is disabled,
 * each non-global PhasedX gate will be replaced with two global NPhasedX, but
 * any other gates will be left untouched.
 *
 * @param squash Whether to squash the circuit before globalisation.
 *
 * If squash=true (default), this transform always returns true. For
 * squash=false, it will return true if the circuit was transformed and
 * false otherwise.
 *
 * It is not recommended to use this pass with symbolic expressions, as in
 * certain cases a blow-up in symbolic expression sizes may occur.
 */
Transform globalise_PhasedX(bool squash = true);

}  // namespace Transforms

}  // namespace tket
