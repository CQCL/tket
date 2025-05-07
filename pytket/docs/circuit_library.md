# pytket.circuit_library

```{eval-rst}
.. currentmodule:: pytket._tket.circuit_library
```

The {py:mod}`~pytket.circuit_library` module is a collection of methods for generating small circuits, mostly of the form of representing one gate in terms of others (as used, for example, within the {py:meth}`~.passes.AutoRebase` pass).

They can be split into a few rough categories:

- Decompositions of multi-qubit gates into other multi-qubit gates (permitting any single-qubit gates), some of which may trade approximation error in the decomposition for reduced gate costs;
- Decompositions of single-qubit gates into universal gatesets, intended to be run after decomposing all multi-qubit gates;
- Constant circuits of 1-3 gates.

Parameters of the gate to be decomposed are passed as method arguments.

```{eval-rst}
.. automodule:: pytket.circuit_library
```

```{eval-rst}
.. automodule:: pytket._tket.circuit_library

   Exact multi-qubit gate decompositions
   -------------------------------------

   .. autofunction:: pytket.circuit_library.BRIDGE_using_CX_0
   .. autofunction:: pytket.circuit_library.BRIDGE_using_CX_1
   .. autofunction:: pytket.circuit_library.CX_using_TK2
   .. autofunction:: pytket.circuit_library.TK2_using_CX
   .. autofunction:: pytket.circuit_library.TK2_using_CX_and_swap
   .. autofunction:: pytket.circuit_library.TK2_using_3xCX
   .. autofunction:: pytket.circuit_library.CX_using_flipped_CX
   .. autofunction:: pytket.circuit_library.CX_using_ECR
   .. autofunction:: pytket.circuit_library.CX_using_ZZMax
   .. autofunction:: pytket.circuit_library.CX_using_ISWAPMax
   .. autofunction:: pytket.circuit_library.CX_using_ISWAPMax_and_swap
   .. autofunction:: pytket.circuit_library.CX_using_ZZPhase
   .. autofunction:: pytket.circuit_library.CX_using_XXPhase_0
   .. autofunction:: pytket.circuit_library.CX_using_XXPhase_1
   .. autofunction:: pytket.circuit_library.CX_using_AAMS
   .. autofunction:: pytket.circuit_library.SWAP_using_CX_0
   .. autofunction:: pytket.circuit_library.SWAP_using_CX_1
   .. autofunction:: pytket.circuit_library.CZ_using_CX
   .. autofunction:: pytket.circuit_library.CY_using_CX
   .. autofunction:: pytket.circuit_library.CH_using_CX
   .. autofunction:: pytket.circuit_library.CV_using_CX
   .. autofunction:: pytket.circuit_library.CVdg_using_CX
   .. autofunction:: pytket.circuit_library.CSX_using_CX
   .. autofunction:: pytket.circuit_library.CSXdg_using_CX
   .. autofunction:: pytket.circuit_library.CS_using_CX
   .. autofunction:: pytket.circuit_library.CSdg_using_CX
   .. autofunction:: pytket.circuit_library.CSWAP_using_CX
   .. autofunction:: pytket.circuit_library.ECR_using_CX
   .. autofunction:: pytket.circuit_library.ZZMax_using_CX
   .. autofunction:: pytket.circuit_library.CRz_using_TK2
   .. autofunction:: pytket.circuit_library.CRz_using_CX
   .. autofunction:: pytket.circuit_library.CRx_using_TK2
   .. autofunction:: pytket.circuit_library.CRx_using_CX
   .. autofunction:: pytket.circuit_library.CRy_using_TK2
   .. autofunction:: pytket.circuit_library.CRy_using_CX
   .. autofunction:: pytket.circuit_library.CU1_using_TK2
   .. autofunction:: pytket.circuit_library.CU1_using_CX
   .. autofunction:: pytket.circuit_library.CU3_using_CX
   .. autofunction:: pytket.circuit_library.ISWAP_using_TK2
   .. autofunction:: pytket.circuit_library.ISWAP_using_CX
   .. autofunction:: pytket.circuit_library.ISWAPMax_using_TK2
   .. autofunction:: pytket.circuit_library.ISWAPMax_using_CX
   .. autofunction:: pytket.circuit_library.XXPhase_using_TK2
   .. autofunction:: pytket.circuit_library.XXPhase_using_CX
   .. autofunction:: pytket.circuit_library.XXPhase_using_AAMS
   .. autofunction:: pytket.circuit_library.YYPhase_using_TK2
   .. autofunction:: pytket.circuit_library.YYPhase_using_CX
   .. autofunction:: pytket.circuit_library.YYPhase_using_AAMS
   .. autofunction:: pytket.circuit_library.ZZPhase_using_TK2
   .. autofunction:: pytket.circuit_library.ZZPhase_using_CX
   .. autofunction:: pytket.circuit_library.ZZPhase_using_AAMS
   .. autofunction:: pytket.circuit_library.XXPhase3_using_TK2
   .. autofunction:: pytket.circuit_library.XXPhase3_using_CX
   .. autofunction:: pytket.circuit_library.ESWAP_using_TK2
   .. autofunction:: pytket.circuit_library.ESWAP_using_CX
   .. autofunction:: pytket.circuit_library.FSim_using_TK2
   .. autofunction:: pytket.circuit_library.FSim_using_CX
   .. autofunction:: pytket.circuit_library.PhasedISWAP_using_TK2
   .. autofunction:: pytket.circuit_library.PhasedISWAP_using_CX
   .. autofunction:: pytket.circuit_library.TK2_using_ZZPhase
   .. autofunction:: pytket.circuit_library.TK2_using_ZZPhase_and_swap
   .. autofunction:: pytket.circuit_library.TK2_using_TK2_or_swap
   .. autofunction:: pytket.circuit_library.TK2_using_TK2
   .. autofunction:: pytket.circuit_library.TK2_using_ZZMax
   .. autofunction:: pytket.circuit_library.TK2_using_ZZMax_and_swap
   .. autofunction:: pytket.circuit_library.TK2_using_ISWAPMax
   .. autofunction:: pytket.circuit_library.TK2_using_ISWAPMax_and_swap
   .. autofunction:: pytket.circuit_library.TK2_using_AAMS
   .. autofunction:: pytket.circuit_library.CCX_modulo_phase_shift
   .. autofunction:: pytket.circuit_library.CCX_normal_decomp
   .. autofunction:: pytket.circuit_library.C3X_normal_decomp
   .. autofunction:: pytket.circuit_library.C4X_normal_decomp
   .. autofunction:: pytket.circuit_library.CnX_vchain_decomp
   .. autofunction:: pytket.circuit_library.ladder_down
   .. autofunction:: pytket.circuit_library.ladder_down_2
   .. autofunction:: pytket.circuit_library.ladder_up

   Approximate multi-qubit gate decompositions
   -------------------------------------------

   .. autofunction:: pytket.circuit_library.approx_TK2_using_1xCX
   .. autofunction:: pytket.circuit_library.approx_TK2_using_2xCX
   .. autofunction:: pytket.circuit_library.approx_TK2_using_1xZZPhase
   .. autofunction:: pytket.circuit_library.approx_TK2_using_2xZZPhase

   Single-qubit gate decompositions
   --------------------------------

   .. autofunction:: pytket.circuit_library.NPhasedX_using_PhasedX
   .. autofunction:: pytket.circuit_library.TK2_using_normalised_TK2
   .. autofunction:: pytket.circuit_library.TK1_to_PhasedXRz
   .. autofunction:: pytket.circuit_library.TK1_to_RzRx
   .. autofunction:: pytket.circuit_library.TK1_to_RxRy
   .. autofunction:: pytket.circuit_library.TK1_to_RzH
   .. autofunction:: pytket.circuit_library.TK1_to_RzSX
   .. autofunction:: pytket.circuit_library.TK1_to_RzXSX
   .. autofunction:: pytket.circuit_library.TK1_to_TK1
   .. autofunction:: pytket.circuit_library.TK1_to_U3
   .. autofunction:: pytket.circuit_library.Rx_using_GPI
   .. autofunction:: pytket.circuit_library.Ry_using_GPI
   .. autofunction:: pytket.circuit_library.Rz_using_GPI
   .. autofunction:: pytket.circuit_library.TK1_using_GPI

   Constant circuits
   -----------------

   .. autofunction:: pytket.circuit_library.X
   .. autofunction:: pytket.circuit_library.CX
   .. autofunction:: pytket.circuit_library.CCX
   .. autofunction:: pytket.circuit_library.X1_CX
   .. autofunction:: pytket.circuit_library.Z0_CX
   .. autofunction:: pytket.circuit_library.BRIDGE
   .. autofunction:: pytket.circuit_library.H_CZ_H
   .. autofunction:: pytket.circuit_library.CX_VS_CX_reduced
   .. autofunction:: pytket.circuit_library.CX_V_CX_reduced
   .. autofunction:: pytket.circuit_library.CX_S_CX_reduced
   .. autofunction:: pytket.circuit_library.CX_V_S_XC_reduced
   .. autofunction:: pytket.circuit_library.CX_S_V_XC_reduced
   .. autofunction:: pytket.circuit_library.CX_XC_reduced
```
