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

#include <pybind11/pybind11.h>

#include "Circuit/CircPool.hpp"
#include "binder_utils.hpp"
#include "typecast.hpp"

namespace py = pybind11;

namespace tket {

void init_library(py::module &m) {
  /* Circuit library */
  py::module_ library_m = m.def_submodule(
      "_library",
      "Library of reusable circuits and circuit generator functions.");
  library_m.def(
      "_BRIDGE_using_CX_0", &CircPool::BRIDGE_using_CX_0,
      "Equivalent to BRIDGE, using four CX, first CX has control on qubit 0");
  library_m.def(
      "_BRIDGE_using_CX_1", &CircPool::BRIDGE_using_CX_1,
      "Equivalent to BRIDGE, using four CX, first CX has control on qubit 1");
  library_m.def(
      "_CX_using_TK2", &CircPool::CX_using_TK2,
      "Equivalent to CX, using a TK2 and single-qubit gates");
  library_m.def(
      "_TK2_using_CX", &CircPool::TK2_using_CX,
      "Given expressions α, β and γ, return circuit equivalent to "
      "TK2(α, β, γ) using up to 3 CX and single-qubit gates.\n\n"
      "The parameters must be normalised to the Weyl chamber, i.e. it must "
      "hold 0.5 ≥ 𝛼 ≥ 𝛽 ≥ |𝛾|.");
  library_m.def(
      "_approx_TK2_using_1xCX", &CircPool::approx_TK2_using_1xCX,
      "Best approximation of TK2 using 1 CX gate and single-qubit gates, using "
      "squared trace fidelity metric. "
      "No parameter is required for this approximation. The returned circuit "
      "will be equivalent to TK2(0.5, 0, 0).");
  library_m.def(
      "_approx_TK2_using_2xCX", &CircPool::approx_TK2_using_2xCX,
      "Best approximation of TK2 using 2 CX gates and single-qubit gates, "
      "using squared trace fidelity metric. "
      "Given expressions α and β, with 0.5 ≥ α ≥ β ≥ 0, return a circuit "
      "equivalent to TK2(α, β, 0).");
  library_m.def(
      "_TK2_using_3xCX", &CircPool::TK2_using_3xCX,
      "Given expressions α, β and γ, return circuit equivalent to "
      "TK2(α, β, γ) using 3 CX and single-qubit gates.\n\n"
      "Prefer using `_TK2_using_CX` unless you wish to explicitly use 3 CX or "
      "if α, β and γ are not normalised to the Weyl chamber.");
  library_m.def(
      "_CX_using_flipped_CX", &CircPool::CX_using_flipped_CX,
      "Equivalent to CX[0,1], using a CX[1,0] and four H gates");
  library_m.def(
      "_CX_using_ECR", &CircPool::CX_using_ECR,
      "Equivalent to CX, using only ECR, Rx and U3 gates");
  library_m.def(
      "_CX_using_ZZMax", &CircPool::CX_using_ZZMax,
      "Equivalent to CX, using only ZZMax, Rx and Rz gates");
  library_m.def(
      "_CX_using_XXPhase_0", &CircPool::CX_using_XXPhase_0,
      "Equivalent to CX, using only XXPhase, Rx, Ry and Rz gates");

  library_m.def(
      "_CX_using_XXPhase_1", &CircPool::CX_using_XXPhase_1,
      "Equivalent to CX, using only XXPhase, Rx, Ry and Rz gates");
  library_m.def(
      "_CX_VS_CX_reduced", &CircPool::CX_VS_CX_reduced,
      "CX-reduced form of CX/V,S/CX");
  library_m.def(
      "_CX_V_CX_reduced", &CircPool::CX_V_CX_reduced,
      "CX-reduced form of CX/V,-/CX");
  library_m.def(
      "_CX_S_CX_reduced", &CircPool::CX_S_CX_reduced,
      "CX-reduced form of CX/-,S/CX (= ZZMax)");
  library_m.def(
      "_CX_V_S_XC_reduced", &CircPool::CX_V_S_XC_reduced,
      "CX-reduced form of CX/V,-/S,-/XC");
  library_m.def(
      "_CX_S_V_XC_reduced", &CircPool::CX_S_V_XC_reduced,
      "CX-reduced form of CX/-,S/-,V/XC");
  library_m.def(
      "_CX_XC_reduced", &CircPool::CX_XC_reduced, "CX-reduced form of CX/XC");
  library_m.def(
      "_SWAP_using_CX_0", &CircPool::SWAP_using_CX_0,
      "Equivalent to SWAP, using three CX, outer CX have control on qubit 0");
  library_m.def(
      "_SWAP_using_CX_1", &CircPool::SWAP_using_CX_1,
      "Equivalent to SWAP, using three CX, outer CX have control on qubit 1");
  library_m.def(
      "_two_Rz1", &CircPool::two_Rz1,
      "A two-qubit circuit with an Rz(1) on each qubit");
  library_m.def("_X1_CX", &CircPool::X1_CX, "X[1]; CX[0,1]");
  library_m.def("_Z0_CX", &CircPool::Z0_CX, "Z[0]; CX[0,1] ");

  library_m.def(
      "_CCX_modulo_phase_shift", &CircPool::CCX_modulo_phase_shift,
      "Equivalent to CCX up to phase shift, using three CX. Warning: this is "
      "not equivalent to CCX up to global phase so cannot be used as a direct "
      "substitution except when the phase reversal can be cancelled. Its "
      "unitary representation is like CCX but with a -1 at the (5,5) "
      "position.");
  library_m.def(
      "_CCX_normal_decomp", &CircPool::CCX_normal_decomp,
      "Equivalent to CCX, using five CX");
  library_m.def(
      "_C3X_normal_decomp", &CircPool::C3X_normal_decomp,
      "Equivalent to CCCX, using 14 CX");
  library_m.def(
      "_C4X_normal_decomp", &CircPool::C4X_normal_decomp,
      "Equivalent to CCCCX, using 36 CX ");
  library_m.def(
      "_ladder_down", &CircPool::ladder_down, "CX[0,1]; CX[2,0]; CCX[0,1,2]");
  library_m.def(
      "_ladder_down_2", &CircPool::ladder_down_2,
      "CX[0,1]; X[0]; X[2]; CCX[0,1,2]");
  library_m.def(
      "_ladder_up", &CircPool::ladder_up, "CCX[0,1,2]; CX[2,0]; CX[2,1]");
  library_m.def("_X", &CircPool::X, "Just an X gate");
  library_m.def("_CX", &CircPool::CX, "Just a CX[0,1] gate");
  library_m.def("_CCX", &CircPool::CCX, "Just a CCX[0,1,2] gate");
  library_m.def("_BRIDGE", &CircPool::BRIDGE, "Just a BRIDGE[0,1,2] gate");
  library_m.def("_H_CZ_H", &CircPool::H_CZ_H, "H[1]; CZ[0,1]; H[1] ");
  library_m.def(
      "_CZ_using_CX", &CircPool::CZ_using_CX,
      "Equivalent to CZ, using CX and single-qubit gates");
  library_m.def(
      "_CY_using_CX", &CircPool::CY_using_CX,
      "Equivalent to CY, using CX and single-qubit gates");
  library_m.def(
      "_CH_using_CX", &CircPool::CH_using_CX,
      "Equivalent to CH, using CX and single-qubit gates");
  library_m.def(
      "_CV_using_CX", &CircPool::CV_using_CX,
      "Equivalent to CV, using CX and single-qubit gates ");
  library_m.def(
      "_CVdg_using_CX", &CircPool::CVdg_using_CX,
      "Equivalent to CVdg, using CX and single-qubit gates");
  library_m.def(
      "_CSX_using_CX", &CircPool::CSX_using_CX,
      "Equivalent to CSX, using CX and single-qubit gates");
  library_m.def(
      "_CSXdg_using_CX", &CircPool::CSXdg_using_CX,
      "Equivalent to CSXdg, using CX and single-qubit gates");
  library_m.def(
      "_CSWAP_using_CX", &CircPool::CSWAP_using_CX,
      "Equivalent to CSWAP, using CX and single-qubit gates ");
  library_m.def(
      "_ECR_using_CX", &CircPool::ECR_using_CX,
      "Equivalent to ECR, using CX, Rx and U3 gates ");
  library_m.def(
      "_ZZMax_using_CX", &CircPool::ZZMax_using_CX,
      "Equivalent to ZZMax, using CX, Rz and U3 gates ");
  library_m.def(
      "_CRz_using_TK2", &CircPool::CRz_using_TK2,
      "Equivalent to CRz, using a TK2 and TK1 gates");
  library_m.def(
      "_CRz_using_CX", &CircPool::CRz_using_CX,
      "Equivalent to CRz, using CX and Rz gates");
  library_m.def(
      "_CRx_using_TK2", &CircPool::CRx_using_TK2,
      "Equivalent to CRx, using a TK2 and TK1 gates");
  library_m.def(
      "_CRx_using_CX", &CircPool::CRx_using_CX,
      "Equivalent to CRx, using CX, H and Rx gates");
  library_m.def(
      "_CRy_using_TK2", &CircPool::CRy_using_TK2,
      "Equivalent to CRy, using a TK2 and TK1 gates");
  library_m.def(
      "_CRy_using_CX", &CircPool::CRy_using_CX,
      "Equivalent to CRy, using CX and Ry gates");
  library_m.def(
      "_CU1_using_TK2", &CircPool::CU1_using_TK2,
      "Equivalent to CU1, using a TK2 and TK1 gates");
  library_m.def(
      "_CU1_using_CX", &CircPool::CU1_using_CX,
      "Equivalent to CU1, using CX and U1 gates");
  library_m.def(
      "_CU3_using_CX", &CircPool::CU3_using_CX,
      "Equivalent to CU1, using CX, U1 and U3 gates");
  library_m.def(
      "_ISWAP_using_TK2", &CircPool::ISWAP_using_TK2,
      "Equivalent to ISWAP, using a TK2 gate");
  library_m.def(
      "_ISWAP_using_CX", &CircPool::ISWAP_using_CX,
      "Equivalent to ISWAP, using CX, U3 and Rz gates");
  library_m.def(
      "_XXPhase_using_TK2", &CircPool::XXPhase_using_TK2,
      "Equivalent to XXPhase, using a TK2 gate");
  library_m.def(
      "_XXPhase_using_CX", &CircPool::XXPhase_using_CX,
      "Equivalent to XXPhase, using CX and U3 gates ");
  library_m.def(
      "_YYPhase_using_TK2", &CircPool::YYPhase_using_TK2,
      "Equivalent to YYPhase, using a TK2 gate");
  library_m.def(
      "_YYPhase_using_CX", &CircPool::YYPhase_using_CX,
      "Equivalent to YYPhase, using CX, Rz and U3 gates");
  library_m.def(
      "_ZZPhase_using_TK2", &CircPool::ZZPhase_using_TK2,
      "Equivalent to ZZPhase, using a TK2 gate");
  library_m.def(
      "_ZZPhase_using_CX", &CircPool::ZZPhase_using_CX,
      "Equivalent to ZZPhase, using CX and Rz gates");
  library_m.def(
      "_TK2_using_ZZPhase", &CircPool::TK2_using_ZZPhase,
      "Equivalent to TK2, using 3 ZZPhase gates");
  library_m.def(
      "_approx_TK2_using_1xZZPhase", &CircPool::approx_TK2_using_1xZZPhase,
      "Approximate equivalent to TK2, using 1 ZZPhase gate and single-qubit "
      "gates. Only requires the first angle of the TK2 gate.");
  library_m.def(
      "_approx_TK2_using_2xZZPhase", &CircPool::approx_TK2_using_2xZZPhase,
      "Approximate equivalent to TK2, using 2 ZZPhase gates and single-qubit "
      "gates. Only requires the first two angles of the TK2 gate.");
  library_m.def(
      "_XXPhase3_using_TK2", &CircPool::XXPhase3_using_TK2,
      "Equivalent to XXPhase3, using three TK2 gates");
  library_m.def(
      "_XXPhase3_using_CX", &CircPool::XXPhase3_using_CX,
      "Equivalent to 3-qubit MS interaction, using CX and U3 gates");
  library_m.def(
      "_ESWAP_using_TK2", &CircPool::ESWAP_using_TK2,
      "Equivalent to ESWAP, using a TK2 and (Clifford) TK1 gates");
  library_m.def(
      "_ESWAP_using_CX", &CircPool::XXPhase3_using_CX,
      "Equivalent to ESWAP, using CX, X, S, Ry and U1 gates");
  library_m.def(
      "_FSim_using_TK2", &CircPool::FSim_using_TK2,
      "Equivalent to FSim, using a TK2 and TK1 gates");
  library_m.def(
      "_FSim_using_CX", &CircPool::FSim_using_CX,
      "Equivalent to Fsim, using CX, X, S, U1 and U3 gates ");
  library_m.def(
      "_PhasedISWAP_using_TK2", &CircPool::PhasedISWAP_using_TK2,
      "Equivalent to PhasedISWAP, using a TK2 and Rz gates");
  library_m.def(
      "_PhasedISWAP_using_CX", &CircPool::PhasedISWAP_using_CX,
      "Equivalent to PhasedISWAP, using CX, U3 and Rz gates");
  library_m.def(
      "_NPhasedX_using_PhasedX", &CircPool::NPhasedX_using_PhasedX,
      "Unwrap NPhasedX, into number_of_qubits PhasedX gates");

  library_m.def(
      "_TK1_to_PhasedXRz", &CircPool::tk1_to_PhasedXRz,
      "A tk1 equivalent circuit given tk1 parameters in terms of PhasedX, Rz");
  library_m.def(
      "_TK1_to_RzRx", &CircPool::tk1_to_rzrx,
      "A tk1 equivalent circuit given tk1 parameters in terms of Rz, Rx");
  library_m.def(
      "_TK1_to_RzH", &CircPool::tk1_to_rzh,
      "A tk1 equivalent circuit given tk1 parameters in terms of Rz, H");
  library_m.def(
      "_TK1_to_RzSX", &CircPool::tk1_to_rzsx,
      "A tk1 equivalent circuit given tk1 parameters in terms of Rz, Sx");
  library_m.def(
      "_TK1_to_TK1", &CircPool::tk1_to_tk1,
      "A circuit of a single tk1 gate with given parameters");
}
}  // namespace tket
