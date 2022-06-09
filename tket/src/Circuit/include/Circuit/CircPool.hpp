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
#include "Utils/Expression.hpp"

namespace tket {

namespace CircPool {

/** Equivalent to BRIDGE, using four CX, first CX has control on qubit 0 */
const Circuit &BRIDGE_using_CX_0();

/** Equivalent to BRIDGE, using four CX, first CX has control on qubit 1 */
const Circuit &BRIDGE_using_CX_1();

/** Equivalent to CX, using a TK2 and single-qubit gates */
const Circuit &CX_using_TK2();

/** Equivalent to CX[0,1], using a CX[1,0] and four H gates */
const Circuit &CX_using_flipped_CX();

/** Equivalent to CX, using only ECR, Rx and U3 gates */
const Circuit &CX_using_ECR();

/** Equivalent to CX, using only ZZMax, Rx and Rz gates */
const Circuit &CX_using_ZZMax();

/** Equivalent to CX, using only XXPhase, Rx, Ry and Rz gates */
const Circuit &CX_using_XXPhase_0();

/** Equivalent to CX, using only XXPhase, Rx and Rz gates */
const Circuit &CX_using_XXPhase_1();

/**
 * CX-reduced form of CX/V,S/CX
 *
 *  --C--V--C--      --Z--S--V--X--S--\ /--
 *    |     |    =>             |      \
 *  --X--S--X--      --X--V--S--C--V--/ \--
 */
const Circuit &CX_VS_CX_reduced();

/**
 * CX-reduced form of CX/V,-/CX
 *
 *  --C--V--C--      --X--V--S--V--C--V--S--
 *    |     |    =>                |
 *  --X-----X--      --V-----------X--------
 */
const Circuit &CX_V_CX_reduced();

/**
 * CX-reduced form of CX/-,S/CX (= ZZMax)
 *
 * --C-----C--      --S-----------C--------
 *   |     |    =>                |
 * --X--S--X--      --Z--S--V--S--X--S--V--
 */
const Circuit &CX_S_CX_reduced();

/**
 * CX-reduced form of CX/V,-/S,-/XC
 *
 * --C--V--S--X--        --S-----C--V--S--
 *   |        |     =>           |
 * --X--------C--        --Z--S--X--S-----
 */
const Circuit &CX_V_S_XC_reduced();

/**
 * CX-reduced form of CX/-,S/-,V/XC
 *
 * --C--------X--      --X--V--C--V-----
 *   |        |    =>          |
 * --X--S--V--C--      --V-----X--S--V--
 */
const Circuit &CX_S_V_XC_reduced();

/**
 * CX-reduced form of CX/XC
 *
 * --C--X--      --X--\ /--
 *   |  |    =>    |   \
 * --X--C--      --C--/ \--
 */
const Circuit &CX_XC_reduced();

/** Equivalent to SWAP, using three CX, outer CX have control on qubit 0 */
const Circuit &SWAP_using_CX_0();

/** Equivalent to SWAP, using three CX, outer CX have control on qubit 1 */
const Circuit &SWAP_using_CX_1();

/** A two-qubit circuit with an Rz(1) on each qubit */
const Circuit &two_Rz1();

/** X[1]; CX[0,1] */
const Circuit &X1_CX();

/** Z[0]; CX[0,1] */
const Circuit &Z0_CX();

/**
 * Equivalent to CCX up to phase shift, using three CX
 *
 * Warning: this is not equivalent to CCX up to global phase so cannot be used
 * as a direct substitution except when the phase reversal can be cancelled. Its
 * unitary representation is like CCX but with a -1 at the (5,5) position.
 */
const Circuit &CCX_modulo_phase_shift();

/** Equivalent to CCX, using five CX */
const Circuit &CCX_normal_decomp();

/** Equivalent to CCCX, using 14 CX */
const Circuit &C3X_normal_decomp();

/** Equivalent to CCCCX, using 36 CX */
const Circuit &C4X_normal_decomp();

/** CX[0,1]; CX[2,0]; CCX[0,1,2] */
const Circuit &ladder_down();

/** CX[0,1]; X[0]; X[2]; CCX[0,1,2] */
const Circuit &ladder_down_2();

/** CCX[0,1,2]; CX[2,0]; CX[2,1] */
const Circuit &ladder_up();

/** Just an X gate */
const Circuit &X();

/** Just a CX[0,1] gate */
const Circuit &CX();

/** Just a CCX[0,1,2] gate */
const Circuit &CCX();

/** Just a BRIDGE[0,1,2] gate */
const Circuit &BRIDGE();

/** H[1]; CZ[0,1]; H[1] */
const Circuit &H_CZ_H();

/** Equivalent to CZ, using CX and single-qubit gates */
const Circuit &CZ_using_CX();

/** Equivalent to CY, using CX and single-qubit gates */
const Circuit &CY_using_CX();

/** Equivalent to CH, using CX and single-qubit gates */
const Circuit &CH_using_CX();

/** Equivalent to CV, using CX and single-qubit gates */
const Circuit &CV_using_CX();

/** Equivalent to CVdg, using CX and single-qubit gates */
const Circuit &CVdg_using_CX();

/** Equivalent to CSX, using CX and single-qubit gates */
const Circuit &CSX_using_CX();

/** Equivalent to CSXdg, using CX and single-qubit gates */
const Circuit &CSXdg_using_CX();

/** Equivalent to CSWAP, using CX and single-qubit gates */
const Circuit &CSWAP_using_CX();

/** Equivalent to ECR, using CX, Rx and U3 gates */
const Circuit &ECR_using_CX();

/** Equivalent to ZZMax, using CX, Rz and U3 gates */
const Circuit &ZZMax_using_CX();

/** Equivalent to CRz, using a TK2 and TK1 gates */
Circuit CRz_using_TK2(const Expr &alpha);

/** Equivalent to CRz, using CX and Rz gates */
Circuit CRz_using_CX(const Expr &alpha);

/** Equivalent to CRx, using a TK2 and TK1 gates */
Circuit CRx_using_TK2(const Expr &alpha);

/** Equivalent to CRx, using CX, H and Rx gates */
Circuit CRx_using_CX(const Expr &alpha);

/** Equivalent to CRy, using a TK2 and TK1 gates */
Circuit CRy_using_TK2(const Expr &alpha);

/** Equivalent to CRy, using CX and Ry gates */
Circuit CRy_using_CX(const Expr &alpha);

/** Equivalent to CU1, using a TK2 and TK1 gates */
Circuit CU1_using_TK2(const Expr &lambda);

/** Equivalent to CU1, using CX and U1 gates */
Circuit CU1_using_CX(const Expr &lambda);

/** Equivalent to CU1, using CX, U1 and U3 gates */
Circuit CU3_using_CX(const Expr &theta, const Expr &phi, const Expr &lambda);

/** Equivalent to ISWAP, using a TK2 gate */
Circuit ISWAP_using_TK2(const Expr &alpha);

/** Equivalent to ISWAP, using CX, U3 and Rz gates */
Circuit ISWAP_using_CX(const Expr &alpha);

/** Equivalent to XXPhase, using a TK2 gate */
Circuit XXPhase_using_TK2(const Expr &alpha);

/** Equivalent to XXPhase, using CX and U3 gates */
Circuit XXPhase_using_CX(const Expr &alpha);

/** Equivalent to YYPhase, using a TK2 gate */
Circuit YYPhase_using_TK2(const Expr &alpha);

/** Equivalent to YYPhase, using CX, Rz and U3 gates */
Circuit YYPhase_using_CX(const Expr &alpha);

/** Equivalent to ZZPhase, using a TK2 gate */
Circuit ZZPhase_using_TK2(const Expr &alpha);

/** Equivalent to ZZPhase, using CX and Rz gates */
Circuit ZZPhase_using_CX(const Expr &alpha);

/**
 * @brief Equivalent to XXPhase, using ZZPhase and H gates.
 *
 * @param alpha The gate parameter to the XXPhase gate.
 * @return Circuit Equivalent circuit using ZZPhase.
 */
Circuit XXPhase_using_ZZPhase(const Expr &alpha);

/**
 * @brief Equivalent to YYPhase, using ZZPhase and V/Vdg gates.
 *
 * @param alpha The gate parameter to the YYPhase gate.
 * @return Circuit Equivalent circuit using ZZPhase.
 */
Circuit YYPhase_using_ZZPhase(const Expr &alpha);

/**
 * @brief Equivalent to TK2(0.5, 0, 0), using a single CX gate.

 * Using 1 CX yields an approximate decomposition of the TK2 gate which is
 * equivalent to a TK2(0.5, 0, 0) gate. This is always the optimal
 * 1-CX approximation of any TK2 gate, with respect to the squared trace
 * fidelity metric.
 *
 * @return Circuit Equivalent circuit to TK2(0.5, 0, 0).
 */
Circuit approx_TK2_using_1xCX();

/**
 * @brief Equivalent to TK2(α, β, 0), using 2 CX gates.
 *
 * Using 2 CX gates we can implement any gate of the form TK2(α, β, 0). This is
 * the optimal 2-CX approximation for any TK2(α, β, γ), with respect to the
 * squared trace fidelity metric.
 *
 * Requires 0.5 ≥ α ≥ β ≥ 0.
 *
 * @return Circuit Equivalent circuit to TK2(α, β, 0).
 */
Circuit approx_TK2_using_2xCX(const Expr &alpha, const Expr &beta);

/**
 * @brief Equivalent to TK2(α, β, γ), using 3 CX gates.
 *
 * This is an exact 3 CX decomposition of the TK2(α, β, γ) gate.
 * Prefer using `TK2_using_CX` unless you wish to explicitly use 3 CX or if
 * α, β and γ are not normalised to the Weyl chamber.
 *
 * @return Circuit Equivalent circuit to TK2(α, β, γ).
 */
Circuit TK2_using_3xCX(const Expr &alpha, const Expr &beta, const Expr &gamma);

/**
 * @brief Equivalent to TK2(α, β, γ) with minimal number of CX gates.
 *
 * A TK2-equivalent circuit with as few CX gates as possible (0, 1, 2 or 3 CX).

 * Decomposition is exact. The parameters must be normalised to the Weyl
 * chamber, i.e. it must hold 0.5 ≥ 𝛼 ≥ 𝛽 ≥ |𝛾|.
 *
 * In cases where hardware gate fidelities are known, it might be sensible to
 * use TK2 decompositions that are inexact but less noisy. See DecomposeTK2
 * pass and transform.
 *
 * @return Circuit Equivalent circuit to TK2(α, β, γ).
 */
Circuit TK2_using_CX(const Expr &alpha, const Expr &beta, const Expr &gamma);

/**
 * @brief Equivalent to TK2(α, 0, 0), using 1 ZZPhase gate.
 *
 * Using 1 ZZPhase gate we can implement any gate of the form TK2(α, 0, 0).
 * This is the optimal 1-ZZPhase approximation for any TK2(α, β, γ), with
 * respect to the squared trace fidelity metric.
 *
 * Requires 0.5 ≥ α ≥ 0.
 *
 * @return Circuit Equivalent circuit to TK2(α, β, 0).
 */
Circuit approx_TK2_using_1xZZPhase(const Expr &alpha);

/**
 * @brief Equivalent to TK2(α, β, 0), using 2 ZZPhase gates.
 *
 * Using 2 ZZPhase gates we can implement any gate of the form TK2(α, β, 0).
 * This is the optimal 2-ZZPhase approximation for any TK2(α, β, γ), with
 * respect to the squared trace fidelity metric.
 *
 * Warning: in practice, we would not expect this decomposition to be attractive
 * on real hardware, as the same approximation fidelity can be obtained using
 * 2 ZZMax gates, which would typically have (same or) higher fidelity than
 * variable angle ZZPhase gates.
 *
 * @return Circuit Equivalent circuit to TK2(α, β, 0).
 */
Circuit approx_TK2_using_2xZZPhase(const Expr &alpha, const Expr &beta);

/**
 * @brief Equivalent to TK2(α, β, γ), using 3 ZZPhase gates.
 *
 * This is an exact 3 ZZPhase decomposition of the TK2(α, β, γ) gate.
 *
 * Warning: in practice, we would not expect this decomposition to be attractive
 * on real hardware, as the same approximation fidelity can be obtained using
 * 3 ZZMax gates, which would typically have (same or) higher fidelity than
 * variable angle ZZPhase gates.
 *
 * @return Circuit Equivalent circuit to TK2(α, β, γ).
 */
Circuit TK2_using_ZZPhase(
    const Expr &alpha, const Expr &beta, const Expr &gamma);

/** Equivalent to XXPhase3, using three TK2 gates */
Circuit XXPhase3_using_TK2(const Expr &alpha);

/** Equivalent to 3-qubit MS interaction, using CX and U3 gates */
Circuit XXPhase3_using_CX(const Expr &alpha);

/** Equivalent to ESWAP, using a TK2 and (Clifford) TK1 gates */
Circuit ESWAP_using_TK2(const Expr &alpha);

/** Equivalent to ESWAP, using CX, X, S, Ry and U1 gates */
Circuit ESWAP_using_CX(const Expr &alpha);

/** Equivalent to FSim, using a TK2 and TK1 gates */
Circuit FSim_using_TK2(const Expr &alpha, const Expr &beta);

/** Equivalent to FSim, using CX, X, S, U1 and U3 gates */
Circuit FSim_using_CX(const Expr &alpha, const Expr &beta);

/** Equivalent to PhasedISWAP, using a TK2 and Rz gates */
Circuit PhasedISWAP_using_TK2(const Expr &p, const Expr &t);

/** Equivalent to PhasedISWAP, using CX, U3 and Rz gates */
Circuit PhasedISWAP_using_CX(const Expr &p, const Expr &t);

/** Unwrap NPhasedX, into number_of_qubits PhasedX gates */
Circuit NPhasedX_using_PhasedX(
    unsigned int number_of_qubits, const Expr &alpha, const Expr &beta);

// converts a TK1 gate to a PhasedXRz gate
Circuit tk1_to_PhasedXRz(
    const Expr &alpha, const Expr &beta, const Expr &gamma);

Circuit tk1_to_rzrx(const Expr &alpha, const Expr &beta, const Expr &gamma);

Circuit tk1_to_rzh(const Expr &alpha, const Expr &beta, const Expr &gamma);

Circuit tk1_to_rzsx(const Expr &alpha, const Expr &beta, const Expr &gamma);

Circuit tk1_to_tk1(const Expr &alpha, const Expr &beta, const Expr &gamma);

}  // namespace CircPool

}  // namespace tket
