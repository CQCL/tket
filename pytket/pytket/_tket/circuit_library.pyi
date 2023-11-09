from __future__ import annotations
import pytket._tket.circuit
import sympy
import typing
__all__ = ['BRIDGE', 'BRIDGE_using_CX_0', 'BRIDGE_using_CX_1', 'C3X_normal_decomp', 'C4X_normal_decomp', 'CCX', 'CCX_modulo_phase_shift', 'CCX_normal_decomp', 'CH_using_CX', 'CRx_using_CX', 'CRx_using_TK2', 'CRy_using_CX', 'CRy_using_TK2', 'CRz_using_CX', 'CRz_using_TK2', 'CSWAP_using_CX', 'CSX_using_CX', 'CSXdg_using_CX', 'CS_using_CX', 'CSdg_using_CX', 'CU1_using_CX', 'CU1_using_TK2', 'CU3_using_CX', 'CV_using_CX', 'CVdg_using_CX', 'CX', 'CX_S_CX_reduced', 'CX_S_V_XC_reduced', 'CX_VS_CX_reduced', 'CX_V_CX_reduced', 'CX_V_S_XC_reduced', 'CX_XC_reduced', 'CX_using_ECR', 'CX_using_TK2', 'CX_using_XXPhase_0', 'CX_using_XXPhase_1', 'CX_using_ZZMax', 'CX_using_ZZPhase', 'CX_using_flipped_CX', 'CY_using_CX', 'CZ_using_CX', 'ECR_using_CX', 'ESWAP_using_CX', 'ESWAP_using_TK2', 'FSim_using_CX', 'FSim_using_TK2', 'H_CZ_H', 'ISWAP_using_CX', 'ISWAP_using_TK2', 'NPhasedX_using_PhasedX', 'PhasedISWAP_using_CX', 'PhasedISWAP_using_TK2', 'SWAP_using_CX_0', 'SWAP_using_CX_1', 'TK1_to_PhasedXRz', 'TK1_to_RzH', 'TK1_to_RzRx', 'TK1_to_RzSX', 'TK1_to_TK1', 'TK2_using_3xCX', 'TK2_using_CX', 'TK2_using_CX_and_swap', 'TK2_using_TK2_or_swap', 'TK2_using_ZZMax', 'TK2_using_ZZMax_and_swap', 'TK2_using_ZZPhase', 'TK2_using_ZZPhase_and_swap', 'TK2_using_normalised_TK2', 'X', 'X1_CX', 'XXPhase3_using_CX', 'XXPhase3_using_TK2', 'XXPhase_using_CX', 'XXPhase_using_TK2', 'YYPhase_using_CX', 'YYPhase_using_TK2', 'Z0_CX', 'ZZMax_using_CX', 'ZZPhase_using_CX', 'ZZPhase_using_TK2', 'approx_TK2_using_1xCX', 'approx_TK2_using_1xZZPhase', 'approx_TK2_using_2xCX', 'approx_TK2_using_2xZZPhase', 'ladder_down', 'ladder_down_2', 'ladder_up', 'two_Rz1']
def BRIDGE() -> pytket._tket.circuit.Circuit:
    """
    Just a BRIDGE[0,1,2] gate
    """
def BRIDGE_using_CX_0() -> pytket._tket.circuit.Circuit:
    """
    Equivalent to BRIDGE, using four CX, first CX has control on qubit 0
    """
def BRIDGE_using_CX_1() -> pytket._tket.circuit.Circuit:
    """
    Equivalent to BRIDGE, using four CX, first CX has control on qubit 1
    """
def C3X_normal_decomp() -> pytket._tket.circuit.Circuit:
    """
    Equivalent to CCCX, using 14 CX
    """
def C4X_normal_decomp() -> pytket._tket.circuit.Circuit:
    """
    Equivalent to CCCCX, using 36 CX 
    """
def CCX() -> pytket._tket.circuit.Circuit:
    """
    Just a CCX[0,1,2] gate
    """
def CCX_modulo_phase_shift() -> pytket._tket.circuit.Circuit:
    """
    Equivalent to CCX up to phase shift, using three CX. Warning: this is not equivalent to CCX up to global phase so cannot be used as a direct substitution except when the phase reversal can be cancelled. Its unitary representation is like CCX but with a -1 at the (5,5) position.
    """
def CCX_normal_decomp() -> pytket._tket.circuit.Circuit:
    """
    Equivalent to CCX, using 6 CX
    """
def CH_using_CX() -> pytket._tket.circuit.Circuit:
    """
    Equivalent to CH, using CX and single-qubit gates
    """
def CRx_using_CX(arg0: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Equivalent to CRx, using CX, H and Rx gates
    """
def CRx_using_TK2(arg0: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Equivalent to CRx, using a TK2 and TK1 gates
    """
def CRy_using_CX(arg0: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Equivalent to CRy, using CX and Ry gates
    """
def CRy_using_TK2(arg0: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Equivalent to CRy, using a TK2 and TK1 gates
    """
def CRz_using_CX(arg0: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Equivalent to CRz, using CX and Rz gates
    """
def CRz_using_TK2(arg0: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Equivalent to CRz, using a TK2 and TK1 gates
    """
def CSWAP_using_CX() -> pytket._tket.circuit.Circuit:
    """
    Equivalent to CSWAP, using CX and single-qubit gates 
    """
def CSX_using_CX() -> pytket._tket.circuit.Circuit:
    """
    Equivalent to CSX, using CX and single-qubit gates
    """
def CSXdg_using_CX() -> pytket._tket.circuit.Circuit:
    """
    Equivalent to CSXdg, using CX and single-qubit gates
    """
def CS_using_CX() -> pytket._tket.circuit.Circuit:
    """
    Equivalent to CS, using CX and single-qubit gates
    """
def CSdg_using_CX() -> pytket._tket.circuit.Circuit:
    """
    Equivalent to CSdg, using CX and single-qubit gates
    """
def CU1_using_CX(arg0: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Equivalent to CU1, using CX and U1 gates
    """
def CU1_using_TK2(arg0: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Equivalent to CU1, using a TK2 and TK1 gates
    """
def CU3_using_CX(arg0: sympy.Expr | float, arg1: sympy.Expr | float, arg2: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Equivalent to CU1, using CX, U1 and U3 gates
    """
def CV_using_CX() -> pytket._tket.circuit.Circuit:
    """
    Equivalent to CV, using CX and single-qubit gates 
    """
def CVdg_using_CX() -> pytket._tket.circuit.Circuit:
    """
    Equivalent to CVdg, using CX and single-qubit gates
    """
def CX() -> pytket._tket.circuit.Circuit:
    """
    Just a CX[0,1] gate
    """
def CX_S_CX_reduced() -> pytket._tket.circuit.Circuit:
    """
    CX-reduced form of CX/-,S/CX (= ZZMax)
    """
def CX_S_V_XC_reduced() -> pytket._tket.circuit.Circuit:
    """
    CX-reduced form of CX/-,S/-,V/XC
    """
def CX_VS_CX_reduced() -> pytket._tket.circuit.Circuit:
    """
    CX-reduced form of CX/V,S/CX
    """
def CX_V_CX_reduced() -> pytket._tket.circuit.Circuit:
    """
    CX-reduced form of CX/V,-/CX
    """
def CX_V_S_XC_reduced() -> pytket._tket.circuit.Circuit:
    """
    CX-reduced form of CX/V,-/S,-/XC
    """
def CX_XC_reduced() -> pytket._tket.circuit.Circuit:
    """
    CX-reduced form of CX/XC
    """
def CX_using_ECR() -> pytket._tket.circuit.Circuit:
    """
    Equivalent to CX, using only ECR, Rx and U3 gates
    """
def CX_using_TK2() -> pytket._tket.circuit.Circuit:
    """
    Equivalent to CX, using a TK2 and single-qubit gates
    """
def CX_using_XXPhase_0() -> pytket._tket.circuit.Circuit:
    """
    Equivalent to CX, using only XXPhase, Rx, Ry and Rz gates
    """
def CX_using_XXPhase_1() -> pytket._tket.circuit.Circuit:
    """
    Equivalent to CX, using only XXPhase, Rx, Ry and Rz gates
    """
def CX_using_ZZMax() -> pytket._tket.circuit.Circuit:
    """
    Equivalent to CX, using only ZZMax, Rx and Rz gates
    """
def CX_using_ZZPhase() -> pytket._tket.circuit.Circuit:
    """
    Equivalent to CX, using only ZZPhase, Rx and Rz gates
    """
def CX_using_flipped_CX() -> pytket._tket.circuit.Circuit:
    """
    Equivalent to CX[0,1], using a CX[1,0] and four H gates
    """
def CY_using_CX() -> pytket._tket.circuit.Circuit:
    """
    Equivalent to CY, using CX and single-qubit gates
    """
def CZ_using_CX() -> pytket._tket.circuit.Circuit:
    """
    Equivalent to CZ, using CX and single-qubit gates
    """
def ECR_using_CX() -> pytket._tket.circuit.Circuit:
    """
    Equivalent to ECR, using CX, Rx and U3 gates 
    """
def ESWAP_using_CX(arg0: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Equivalent to ESWAP, using CX, X, S, Ry and U1 gates
    """
def ESWAP_using_TK2(arg0: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Equivalent to ESWAP, using a TK2 and (Clifford) TK1 gates
    """
def FSim_using_CX(arg0: sympy.Expr | float, arg1: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Equivalent to Fsim, using CX, X, S, U1 and U3 gates 
    """
def FSim_using_TK2(arg0: sympy.Expr | float, arg1: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Equivalent to FSim, using a TK2 and TK1 gates
    """
def H_CZ_H() -> pytket._tket.circuit.Circuit:
    """
    H[1]; CZ[0,1]; H[1] 
    """
def ISWAP_using_CX(arg0: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Equivalent to ISWAP, using CX, U3 and Rz gates
    """
def ISWAP_using_TK2(arg0: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Equivalent to ISWAP, using a TK2 gate
    """
def NPhasedX_using_PhasedX(arg0: int, arg1: sympy.Expr | float, arg2: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Unwrap NPhasedX, into number_of_qubits PhasedX gates
    """
def PhasedISWAP_using_CX(arg0: sympy.Expr | float, arg1: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Equivalent to PhasedISWAP, using CX, U3 and Rz gates
    """
def PhasedISWAP_using_TK2(arg0: sympy.Expr | float, arg1: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Equivalent to PhasedISWAP, using a TK2 and Rz gates
    """
def SWAP_using_CX_0() -> pytket._tket.circuit.Circuit:
    """
    Equivalent to SWAP, using three CX, outer CX have control on qubit 0
    """
def SWAP_using_CX_1() -> pytket._tket.circuit.Circuit:
    """
    Equivalent to SWAP, using three CX, outer CX have control on qubit 1
    """
def TK1_to_PhasedXRz(arg0: sympy.Expr | float, arg1: sympy.Expr | float, arg2: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    A tk1 equivalent circuit given tk1 parameters in terms of PhasedX, Rz
    """
def TK1_to_RzH(arg0: sympy.Expr | float, arg1: sympy.Expr | float, arg2: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    A tk1 equivalent circuit given tk1 parameters in terms of Rz, H
    """
def TK1_to_RzRx(arg0: sympy.Expr | float, arg1: sympy.Expr | float, arg2: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    A tk1 equivalent circuit given tk1 parameters in terms of Rz, Rx
    """
def TK1_to_RzSX(arg0: sympy.Expr | float, arg1: sympy.Expr | float, arg2: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    A tk1 equivalent circuit given tk1 parameters in terms of Rz, Sx
    """
def TK1_to_TK1(arg0: sympy.Expr | float, arg1: sympy.Expr | float, arg2: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    A circuit of a single tk1 gate with given parameters
    """
def TK2_using_3xCX(arg0: sympy.Expr | float, arg1: sympy.Expr | float, arg2: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Given expressions α, β and γ, return circuit equivalent to TK2(α, β, γ) using 3 CX and single-qubit gates.
    
    Prefer using `_TK2_using_CX` unless you wish to explicitly use 3 CX or if α, β and γ are not normalised to the Weyl chamber.
    """
def TK2_using_CX(arg0: sympy.Expr | float, arg1: sympy.Expr | float, arg2: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Given expressions α, β and γ, return circuit equivalent to TK2(α, β, γ) using up to 3 CX and single-qubit gates.
    
    The decomposition minimizes the number of CX gates.
    """
def TK2_using_CX_and_swap(arg0: sympy.Expr | float, arg1: sympy.Expr | float, arg2: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Given expressions α, β and γ, return circuit equivalent to TK2(α, β, γ), up to a wire swap that is encoded in the implicit qubit permutation of the Circuit, using up to 3 CX and single-qubit gates.
    
    The decomposition minimizes the number of CX gates.
    """
def TK2_using_TK2_or_swap(arg0: sympy.Expr | float, arg1: sympy.Expr | float, arg2: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Either the exact TK2, or a wire swap encoded in the implicit qubit permutation of the Circuit and single qubit gates.
    """
def TK2_using_ZZMax(arg0: sympy.Expr | float, arg1: sympy.Expr | float, arg2: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Equivalent to TK2, using up to 3 ZZMax gates.
    """
def TK2_using_ZZMax_and_swap(arg0: sympy.Expr | float, arg1: sympy.Expr | float, arg2: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Equivalent to TK2, up to a wire swap that is encoded in the implicit qubit permutation of the Circuit, using up to 3 ZZMax gates.
    """
def TK2_using_ZZPhase(arg0: sympy.Expr | float, arg1: sympy.Expr | float, arg2: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Equivalent to TK2, using 3 ZZPhase gates
    """
def TK2_using_ZZPhase_and_swap(arg0: sympy.Expr | float, arg1: sympy.Expr | float, arg2: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Equivalent to TK2, up to a wire swap that is encoded in the implicit qubit permutation of the Circuit, using up to 3 ZZPhase gates.
    """
def TK2_using_normalised_TK2(arg0: sympy.Expr | float, arg1: sympy.Expr | float, arg2: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    TK2(a, b, c)-equivalent circuit, using a single normalised TK2 and single-qb gates
    """
def X() -> pytket._tket.circuit.Circuit:
    """
    Just an X gate
    """
def X1_CX() -> pytket._tket.circuit.Circuit:
    """
    X[1]; CX[0,1]
    """
def XXPhase3_using_CX(arg0: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Equivalent to 3-qubit MS interaction, using CX and U3 gates
    """
def XXPhase3_using_TK2(arg0: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Equivalent to XXPhase3, using three TK2 gates
    """
def XXPhase_using_CX(arg0: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Equivalent to XXPhase, using CX and U3 gates 
    """
def XXPhase_using_TK2(arg0: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Equivalent to XXPhase, using a TK2 gate
    """
def YYPhase_using_CX(arg0: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Equivalent to YYPhase, using two CX gates and one Ry, one Sdg and one S gate.
    """
def YYPhase_using_TK2(arg0: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Equivalent to YYPhase, using a TK2 gate
    """
def Z0_CX() -> pytket._tket.circuit.Circuit:
    """
    Z[0]; CX[0,1] 
    """
def ZZMax_using_CX() -> pytket._tket.circuit.Circuit:
    """
    Equivalent to ZZMax, using CX, Rz and U3 gates 
    """
def ZZPhase_using_CX(arg0: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Equivalent to ZZPhase, using CX and Rz gates
    """
def ZZPhase_using_TK2(arg0: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Equivalent to ZZPhase, using a TK2 gate
    """
def approx_TK2_using_1xCX() -> pytket._tket.circuit.Circuit:
    """
    Best approximation of TK2 using 1 CX gate and single-qubit gates, using squared trace fidelity metric. No parameter is required for this approximation. The returned circuit will be equivalent to TK2(0.5, 0, 0).
    """
def approx_TK2_using_1xZZPhase(arg0: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Approximate equivalent to TK2, using 1 ZZPhase gate and single-qubit gates. Only requires the first angle of the TK2 gate.
    """
def approx_TK2_using_2xCX(arg0: sympy.Expr | float, arg1: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Best approximation of TK2 using 2 CX gates and single-qubit gates, using squared trace fidelity metric. Given expressions α and β, with 0.5 ≥ α ≥ β ≥ 0, return a circuit equivalent to TK2(α, β, 0).
    """
def approx_TK2_using_2xZZPhase(arg0: sympy.Expr | float, arg1: sympy.Expr | float) -> pytket._tket.circuit.Circuit:
    """
    Approximate equivalent to TK2, using 2 ZZPhase gates and single-qubit gates. Only requires the first two angles of the TK2 gate.
    """
def ladder_down() -> pytket._tket.circuit.Circuit:
    """
    CX[0,1]; CX[2,0]; CCX[0,1,2]
    """
def ladder_down_2() -> pytket._tket.circuit.Circuit:
    """
    CX[0,1]; X[0]; X[2]; CCX[0,1,2]
    """
def ladder_up() -> pytket._tket.circuit.Circuit:
    """
    CCX[0,1,2]; CX[2,0]; CX[2,1]
    """
def two_Rz1() -> pytket._tket.circuit.Circuit:
    """
    A two-qubit circuit with an Rz(1) on each qubit
    """
