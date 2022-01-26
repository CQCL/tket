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

#include "Transform.hpp"

namespace tket {

namespace Transforms {

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

// finds all ZZ gates which are in the form of CX and Rz gates and converts
// into ZZ gates Expects: CX, Rz Produces: ZZ, CX, Rz
Transform decompose_ZZPhase();

// identifies any TK1, Rz,Ry,Rx,u3,u2,u1 gate sequences that can be naively
// decomposed into S,Z,X,V returns clifford sequences in a standard form
// (Z)(X)(S)(V)(S)
// Expects: TK1, U3, U2, U1, Rx, Ry, Rz, and any multi-qubit gates
// Produces: Z, X, S, V, and any remaining non-clifford single-qubit gates,
// and any multi-qubit gates
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
 * returns: H + PhasePolyBox
 */
Transform compose_phase_poly_boxes();

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

}  // namespace Transforms

}  // namespace tket
