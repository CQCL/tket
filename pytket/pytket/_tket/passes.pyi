from typing import Any
from __future__ import annotations
import pytket._tket.architecture
import pytket._tket.circuit
import pytket._tket.mapping
import pytket._tket.placement
import pytket._tket.predicates
import pytket._tket.transform
import pytket._tket.unit_id
import sympy
import typing
__all__ = ['AASRouting', 'Audit', 'BasePass', 'CNotSynthType', 'CXMappingPass', 'CliffordSimp', 'CnXPairwiseDecomposition', 'CommuteThroughMultis', 'ComposePhasePolyBoxes', 'ContextSimp', 'CustomPass', 'CustomRoutingPass', 'DecomposeArbitrarilyControlledGates', 'DecomposeBoxes', 'DecomposeClassicalExp', 'DecomposeMultiQubitsCX', 'DecomposeSingleQubitsTK1', 'DecomposeSwapsToCXs', 'DecomposeSwapsToCircuit', 'DecomposeTK2', 'Default', 'DefaultMappingPass', 'DelayMeasures', 'EulerAngleReduction', 'FlattenRegisters', 'FlattenRelabelRegistersPass', 'FullMappingPass', 'FullPeepholeOptimise', 'GlobalisePhasedX', 'GuidedPauliSimp', 'HamPath', 'KAKDecomposition', 'NaivePlacementPass', 'NormaliseTK2', 'OptimisePhaseGadgets', 'PauliExponentials', 'PauliSimp', 'PauliSquash', 'PeepholeOptimise2Q', 'PlacementPass', 'RebaseCustom', 'RebaseTket', 'Rec', 'RemoveBarriers', 'RemoveDiscarded', 'RemoveImplicitQubitPermutation', 'RemoveRedundancies', 'RenameQubitsPass', 'RepeatPass', 'RepeatUntilSatisfiedPass', 'RepeatWithMetricPass', 'RoundAngles', 'RoutingPass', 'SWAP', 'SafetyMode', 'SequencePass', 'SimplifyInitial', 'SimplifyMeasured', 'SquashCustom', 'SquashRzPhasedX', 'SquashTK1', 'SynthesiseHQS', 'SynthesiseOQC', 'SynthesiseTK', 'SynthesiseTket', 'SynthesiseUMD', 'ThreeQubitSquash', 'ZXGraphlikeOptimisation', 'ZZPhaseToRz']
class BasePass:
    """
    Base class for passes.
    """
    @staticmethod
    def from_dict(arg0: dict) -> BasePass:
        """
        Construct a new Pass instance from a JSON serializable dictionary representation.
        """
    def __getstate__(self) -> tuple:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, arg0: tuple) -> None:
        ...
    def __str__(self) -> str:
        ...
    @typing.overload
    def apply(self, compilation_unit: pytket._tket.predicates.CompilationUnit, safety_mode: SafetyMode = SafetyMode.Default) -> bool:
        """
        Apply to a :py:class:`CompilationUnit`.
        
        :return: True if the pass modified the circuit. Note that in some cases the method may return True even when the circuit is unmodified (but a return value of False definitely implies no modification).
        """
    @typing.overload
    def apply(self, circuit: pytket._tket.circuit.Circuit) -> bool:
        """
        Apply to a :py:class:`Circuit` in-place.
        
        :return: True if pass modified the circuit, else False
        """
    @typing.overload
    def apply(self, circuit: pytket._tket.circuit.Circuit, before_apply: typing.Callable[[pytket._tket.predicates.CompilationUnit, object], None], after_apply: typing.Callable[[pytket._tket.predicates.CompilationUnit, object], None]) -> bool:
        """
        Apply to a :py:class:`Circuit` in-place and invoke callbacks for all nested passes.
        
        
        :param before_apply: Invoked before a pass is applied. The CompilationUnit and a summary of the pass configuration are passed into the callback.
        :param after_apply: Invoked after a pass is applied. The CompilationUnit and a summary of the pass configuration are passed into the callback.
        :return: True if pass modified the circuit, else False
        """
    def to_dict(self) -> dict:
        """
        :return: A JSON serializable dictionary representation of the Pass.
        """
class CNotSynthType:
    """
    Members:
    
      SWAP : swap-based algorithm for CNOT synthesis
    
      HamPath : Hamilton-path-based method for CNOT synthesis; this method will fail if there is no Hamilton path in the given architecture
    
      Rec : recursive Steiner--Gauss method for CNOT synthesis
    """
    HamPath: typing.ClassVar[CNotSynthType]  # value = <CNotSynthType.HamPath: 1>
    Rec: typing.ClassVar[CNotSynthType]  # value = <CNotSynthType.Rec: 2>
    SWAP: typing.ClassVar[CNotSynthType]  # value = <CNotSynthType.SWAP: 0>
    __members__: typing.ClassVar[dict[str, CNotSynthType]]  # value = {'SWAP': <CNotSynthType.SWAP: 0>, 'HamPath': <CNotSynthType.HamPath: 1>, 'Rec': <CNotSynthType.Rec: 2>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class RepeatPass(BasePass):
    """
    Repeat a pass until its `apply()` method returns False, or if `strict_check` is True until it stops modifying the circuit.
    """
    def __init__(self, compilation_pass: BasePass, strict_check: bool = False) -> None:
        """
        Construct from a compilation pass.
        """
    def __str__(self) -> str:
        ...
    def get_pass(self) -> BasePass:
        """
        :return: The underlying compilation pass.
        """
class RepeatUntilSatisfiedPass(BasePass):
    """
    Repeat a compilation pass until a predicate on the circuit is satisfied.
    """
    @typing.overload
    def __init__(self, compilation_pass: BasePass, predicate: pytket._tket.predicates.Predicate) -> None:
        """
        Construct from a compilation pass and a predicate.
        """
    @typing.overload
    def __init__(self, compilation_pass: BasePass, check_function: typing.Callable[[pytket._tket.circuit.Circuit], bool]) -> None:
        """
        Construct from a compilation pass and a user-defined function from :py:class:`Circuit` to `bool`.
        """
    def __str__(self) -> str:
        ...
    def get_pass(self) -> BasePass:
        """
        :return: The underlying compilation pass.
        """
    def get_predicate(self) -> pytket._tket.predicates.Predicate:
        """
        :return: The underlying predicate.
        """
class RepeatWithMetricPass(BasePass):
    """
    Repeat a compilation pass until the given metric stops decreasing.
    """
    def __init__(self, compilation_pass: BasePass, metric: typing.Callable[[pytket._tket.circuit.Circuit], int]) -> None:
        """
        Construct from a compilation pass and a metric function.
        """
    def __str__(self) -> str:
        ...
    def get_metric(self) -> typing.Callable[[pytket._tket.circuit.Circuit], int]:
        """
        :return: The underlying metric.
        """
    def get_pass(self) -> BasePass:
        """
        :return: The underlying compilation pass.
        """
class SafetyMode:
    """
    Members:
    
      Audit : Checks which predicates a circuit satisfies after the application of each base pass
    
      Default : Only check that a circuit satisfies the preconditions of the overall pass at the start and the postconditions at the end
    """
    Audit: typing.ClassVar[SafetyMode]  # value = <SafetyMode.Audit: 0>
    Default: typing.ClassVar[SafetyMode]  # value = <SafetyMode.Default: 1>
    __members__: typing.ClassVar[dict[str, SafetyMode]]  # value = {'Audit': <SafetyMode.Audit: 0>, 'Default': <SafetyMode.Default: 1>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class SequencePass(BasePass):
    """
    A sequence of compilation passes.
    """
    def __init__(self, pass_list: typing.Sequence[BasePass]) -> None:
        """
        Construct from a list of compilation passes arranged in order of application.
        """
    def __str__(self) -> str:
        ...
    def get_sequence(self) -> list[BasePass]:
        """
        :return: The underlying sequence of passes.
        """
def AASRouting(arc: pytket._tket.architecture.Architecture, **kwargs: Any) -> BasePass:
    """
    Construct a pass to relabel :py:class:`Circuit` Qubits to :py:class:`Device` Nodes, and then use architecture-aware synthesis to route the circuit. In the steps of the pass the circuit will be converted to CX, Rz, H gateset. The limited connectivity of the :py:class:`Architecture` is used for the routing. The direction of the edges is ignored. The placement used is GraphPlacement. This pass can take a few parameters for the routing, described below:
    
    - (unsigned) lookahead=1: parameter for the recursive iteration
    - (CNotSynthType) cnotsynthtype=CNotSynthType.Rec: CNOT synthesis type
    
    NB: The circuit needs to have at most as many qubits as the architecture has nodes. The resulting circuit will always have the same number of qubits as the architecture has nodes, even if the input circuit had fewer.
    
    :param arc: target architecture
    :param \**kwargs: parameters for routing (described above)
    :return: a pass to perform the remapping
    """
def CXMappingPass(arc: pytket._tket.architecture.Architecture, placer: pytket._tket.placement.Placement, **kwargs: Any) -> BasePass:
    """
    Construct a pass to convert all gates to CX, relabel :py:class:`Circuit` Qubits to :py:class:`Architecture` Nodes, route to the connectivty graph of a :py:class:`Architecture` and decompose additional routing gates (SWAP and BRIDGE) to CX gates.
    
    :param arc: The Architecture used for connectivity information.
    :param placer: The placement used for relabelling.
    :param \**kwargs: Parameters for routing: (bool)directed_cx=false, (bool)delay_measures=true
    :return: a pass to perform the remapping
    """
def CliffordSimp(allow_swaps: bool = True) -> BasePass:
    """
    An optimisation pass that performs a number of rewrite rules for simplifying Clifford gate sequences, similar to Duncan & Fagan (https://arxiv.org/abs/1901.10114). Given a circuit with CXs and any single-qubit gates, produces a circuit with TK1, CX gates.
    
    :param allow_swaps: dictates whether the rewriting will disregard CX placement or orientation and introduce wire swaps.
    :return: a pass to perform the rewriting
    """
def CnXPairwiseDecomposition() -> BasePass:
    """
    Decompose CnX gates to 2-qubit gates and single qubit gates. For every two CnX gates, reorder their control qubits to improve the chance of gate cancellation
    """
def CommuteThroughMultis() -> BasePass:
    """
    Moves single-qubit operations past multi-qubit operations that they commute with, towards the front of the circuit.
    """
def ComposePhasePolyBoxes(min_size: int = 0) -> BasePass:
    """
    Pass to convert a given :py:class:`Circuit` to the CX, Rz, H gateset and compose phase polynomial boxes from the groups of the CX+Rz gates.
    
    - (unsigned) min_size=0: minimal number of CX gates in each phase polynominal box: groups with a smaller number of CX gates are not affected by this transformation
    
    :param \**kwargs: parameters for composition (described above)
    :return: a pass to perform the composition
    """
def ContextSimp(allow_classical: bool = True, xcirc: pytket._tket.circuit.Circuit | None = None) -> BasePass:
    """
    Applies simplifications enabled by knowledge of qubit state and discarded qubits.
    
    :param allow_classical: allow replacement of measurements on known state with classical set-bit operations
    :param xcirc: 1-qubit circuit implementing an X gate in the transformed circuit (if omitted, an X gate is used)
    :return: a pass to perform the simplification
    """
def CustomPass(transform: typing.Callable[[pytket._tket.circuit.Circuit], pytket._tket.circuit.Circuit], label: str = '') -> BasePass:
    """
    Generate a custom pass from a user-provided circuit transfomation function.
    
    It is the caller's responsibility to provide a valid transform.
    
    :param transform: function taking a :py:class:`Circuit` as an argument and returning a new transformed circuit
    :param label: optional label for the pass
    :return: a pass to perform the transformation
    """
def CustomRoutingPass(arc: pytket._tket.architecture.Architecture, config: typing.Sequence[pytket._tket.mapping.RoutingMethod]) -> BasePass:
    """
    Construct a pass to route to the connectivity graph of an :py:class:`Architecture`. Edge direction is ignored. 
    
    :return: a pass that routes to the given device architecture 
    """
def DecomposeArbitrarilyControlledGates() -> BasePass:
    """
    Decomposes CCX, CnX, CnY, CnZ, and CnRy gates into CX and single-qubit gates.
    """
def DecomposeBoxes(excluded_types: set[pytket._tket.circuit.OpType] = set(), excluded_opgroups: set[str] = set()) -> BasePass:
    """
    Recursively replaces all boxes by their decomposition into circuits.
    
    :param excluded_types: box `OpType`s excluded from decomposition
    :param excluded_opgroups: opgroups excluded from decomposition
    """
def DecomposeClassicalExp() -> BasePass:
    """
    Replaces each :py:class:`ClassicalExpBox` by a sequence of classical gates.
    """
def DecomposeMultiQubitsCX() -> BasePass:
    """
    Converts all multi-qubit gates into CX and single-qubit gates.
    """
def DecomposeSingleQubitsTK1() -> BasePass:
    """
    Converts all single-qubit gates into TK1 gates.
    """
def DecomposeSwapsToCXs(arc: pytket._tket.architecture.Architecture, respect_direction: bool = False) -> BasePass:
    """
    Construct a pass to decompose SWAP and BRIDGE gates to CX gates, constraining connectivity to an :py:class:`Architecture`, optionally taking the directedness of the connectivity graph into account.
    
    :param arc: The architecture to use for connectivity information.
    :param respect_direction: Optionally takes the directedness of the connectivity graph into account.
    :return: a pass to perform the decomposition
    """
def DecomposeSwapsToCircuit(replacement_circuit: pytket._tket.circuit.Circuit) -> BasePass:
    """
    :param replacement_circuit: An equivalent circuit to replace a SWAP gate with in the desired basis.
    :return: a pass to replace all SWAP gates with the given circuit
    """
def DecomposeTK2(allow_swaps: bool = True, **kwargs: Any) -> BasePass:
    """
    Decompose each TK2 gate into two-qubit gates.
    
    Gate fidelities can be passed as keyword arguments to perform noise-aware decompositions. If the fidelities of several gate types are provided, the best will be chosen.
    
    We currently support `CX_fidelity`, `ZZMax_fidelity` and `ZZPhase_fidelity`. If provided, the `CX` and `ZZMax` fidelities must be given by a single floating point fidelity. The `ZZPhase` fidelity is given as a lambda float -> float, mapping a ZZPhase angle parameter to its fidelity. These parameters will be used to return the optimal decomposition of each TK2 gate, taking noise into consideration.
    
    If no fidelities are provided, the TK2 gates will be decomposed exactly using CX gates. For equal fidelities, ZZPhase will be prefered over ZZMax and CX if the decomposition results in fewer two-qubit gates.
    
    All TK2 gate parameters must be normalised, i.e. they must satisfy `NormalisedTK2Predicate`.
    
    Using the `allow_swaps=True` (default) option, qubits will be swapped when convenient to reduce the two-qubit gate count of the decomposed TK2.
    
    If the TK2 angles are symbolic values, the decomposition will be exact (i.e. not noise-aware). It is not possible in general to obtain optimal decompositions for arbitrary symbolic parameters, so consider substituting for concrete values if possible.
    
    :param allow_swaps: Whether to allow implicit wire swaps.
    """
def DefaultMappingPass(arc: pytket._tket.architecture.Architecture, delay_measures: bool = True) -> BasePass:
    """
    Construct a pass to relabel :py:class:`Circuit` Qubits to :py:class:`Architecture` Nodes, and then route to the connectivity graph of the given :py:class:`Architecture`. Edge direction is ignored. Placement used is GraphPlacement.
    
    :param arc: The Architecture used for connectivity information.
    :param delay_measures: Whether to commute measurements to the end of the circuit, defaulting to true.
    :return: a pass to perform the remapping
    """
def DelayMeasures(allow_partial: bool = True) -> BasePass:
    """
    Commutes Measure operations to the end of the circuit. Throws an exception when this is not possible because of gates following the measure which are dependent on either the resulting quantum state or classical values.
    
    :param allow_partial: Whether to allow measurements that cannot be commuted to the end, and delay them as much as possible instead. If false, the pass includes a :py:class:`CommutableMeasuresPredicate` precondition.
    """
def EulerAngleReduction(q: pytket._tket.circuit.OpType, p: pytket._tket.circuit.OpType, strict: bool = False) -> BasePass:
    """
    Uses Euler angle decompositions to squash all chains of P and Q rotations, where P,Q ∈ {Rx,Ry,Rz}. By default (`strict=False`), this pass will try to decompose the chains into pairs of -P-Q- or -Q-P- rotations, commuting any third rotation past multi-qubit gates. If `strict=True`, all chains will be decomposed to P-Q-P triples and no further optimisation is performed.
    
    :param q: The type of the Q rotation (Q ∈ {Rx,Ry,Rz}).
    :param p: The type of the P rotation (P ∈ {Rx,Ry,Rz}, P ≠ Q).
    :param strict: Optionally performs strict P-Q-P Euler decomposition
    :return: a pass that squashes chains of P and Q rotations
    """
def FlattenRegisters() -> BasePass:
    """
    Merges all quantum and classical registers into their respective default registers with contiguous indexing.
    """
def FlattenRelabelRegistersPass(label: str = 'q') -> BasePass:
    """
    Removes empty Quantum wires from the Circuit and relabels all Qubit to a register from passed name. 
    
    :param label: Name to relabel remaining Qubit to, default 'q'.
    :return: A pass that removes empty wires and relabels.
    """
def FullMappingPass(arc: pytket._tket.architecture.Architecture, placer: pytket._tket.placement.Placement, config: typing.Sequence[pytket._tket.mapping.RoutingMethod]) -> BasePass:
    """
    Construct a pass to relabel :py:class:`Circuit` Qubits to :py:class:`Architecture` Nodes, and then route to the connectivity graph of an :py:class:`Architecture`. Edge direction is ignored.
    
    :param arc: The architecture to use for connectivity information. 
    :param placer: The Placement used for relabelling.
    :param config: Parameters for routing, a list of RoutingMethod, each method is checked and run if applicable in turn.
    :return: a pass to perform the remapping
    """
def FullPeepholeOptimise(allow_swaps: bool = True, target_2qb_gate: pytket._tket.circuit.OpType = pytket._tket.circuit.OpType.CX) -> BasePass:
    """
    Performs peephole optimisation including resynthesis of 2- and 3-qubit gate sequences, and converts to a circuit containing only the given 2-qubit gate (which may be CX or TK2) and TK1 gates.
    
    :param allow_swaps: whether to allow implicit wire swaps
    """
def GlobalisePhasedX(squash: bool = True) -> BasePass:
    """
    Turns all PhasedX and NPhasedX gates into global gates
    
    Replaces any PhasedX gates with global NPhasedX gates. By default, this transform will squash all single-qubit gates to PhasedX and Rz gates before proceeding further. Existing non-global NPhasedX will not be preserved. This is the recommended setting for best performance. If squashing is disabled, each non-global PhasedX gate will be replaced with two global NPhasedX, but any other gates will be left untouched.
    
    :param squash: Whether to squash the circuit in pre-processing (default: true).
    
    If squash=true (default), the `GlobalisePhasedX` transform's `apply` method will always return true. For squash=false, `apply()` will return true if the circuit was changed and false otherwise.
    
    It is not recommended to use this pass with symbolic expressions, as in certain cases a blow-up in symbolic expression sizes may occur.
    """
def GuidedPauliSimp(strat: pytket._tket.transform.PauliSynthStrat = pytket._tket.transform.PauliSynthStrat.Sets, cx_config: pytket._tket.circuit.CXConfigType = pytket._tket.circuit.CXConfigType.Snake) -> BasePass:
    """
    Applies the ``PauliSimp`` optimisation pass to any region of the circuit contained within a :py:class:`CircBox`. This can be useful to focus the synthesis to target specific sets of commuting operations, rather than the default greedy approach.
    
    :param strat: A synthesis strategy for the Pauli graph.
    :param cx_config: A configuration of CXs to convert Pauli gadgets into.
    :return: a pass to perform the simplification
    """
@typing.overload
def KAKDecomposition(target_2qb_gate: pytket._tket.circuit.OpType = pytket._tket.circuit.OpType.CX, cx_fidelity: float = 1.0, allow_swaps: bool = True) -> BasePass:
    """
    Squash sequences of two-qubit operations into minimal form.
    
    Pass to squash together sequences of single- and two-qubit gates into minimal form. Can decompose to TK2 or CX gates.
    
    Two-qubit operations can always be expressed in a minimal form of maximum three CXs, or as a single TK2 gate (a result also known as the KAK or Cartan decomposition).
    
    It is in general recommended to squash to TK2 gates, and to then use the `DecomposeTK2` pass for noise-aware decompositions to other gatesets. For backward compatibility, decompositions to CX are also supported. In this case, `cx_fidelity` can be provided to perform approximate decompositions to CX gates.
    
    When decomposing to TK2 gates, any sequence of two or more two-qubit gates on the same set of qubits are replaced by a single TK2 gate. When decomposing to CX, the substitution is only performed if it results in a reduction of the number of CX gates, or if at least one of the two-qubit gates is not a CX.
    
    Using the `allow_swaps=True` (default) option, qubits will be swapped when convenient to further reduce the two-qubit gate count (only applicable when decomposing to CX gates).
    
    :param target_2qb_gate: OpType to decompose to. Either TK2 or CX.
    :param cx_fidelity: Estimated CX gate fidelity, used when target_2qb_gate=CX.
    :param allow_swaps: Whether to allow implicit wire swaps.
    """
@typing.overload
def KAKDecomposition(cx_fidelity: float) -> BasePass:
    ...
def NaivePlacementPass(architecture: pytket._tket.architecture.Architecture) -> BasePass:
    """
    :param architecture: The Architecture used for relabelling.
    :return: a pass to relabel :py:class:`Circuit` Qubits to :py:class:`Architecture` Nodes
    """
def NormaliseTK2() -> BasePass:
    """
    Normalises all TK2 gates.
    
    TK2 gates have three angles in the interval [0, 4], but these can always be normalised to be within the so-called Weyl chamber by adding single-qubit gates.
    
    More precisely, the three angles a, b, c of TK2(a, b, c) are normalised exactly when the two following conditions are met:
     - numerical values must be in the Weyl chamber, ie `1/2 >= a >= b >= |c|`,
     - symbolic values must come before any numerical value in the array.
    
    After this pass, all TK2 angles will be normalised and the circuit will satisfy `NormalisedTK2Predicate`.
    """
def OptimisePhaseGadgets(cx_config: pytket._tket.circuit.CXConfigType = pytket._tket.circuit.CXConfigType.Snake) -> BasePass:
    """
    Construct a pass that synthesises phase gadgets and converts to a circuit containing only CX, TK1 and Phase gates.
    
    :param cx_config: A configuration of CXs to convert phase gadgets into.
    :return: a pass to perform the synthesis
    """
def PauliExponentials(strat: pytket._tket.transform.PauliSynthStrat = pytket._tket.transform.PauliSynthStrat.Sets, cx_config: pytket._tket.circuit.CXConfigType = pytket._tket.circuit.CXConfigType.Snake) -> BasePass:
    """
    Construct a pass that converts a circuit into a graph of Pauli exponential boxes, with information
    
    :param strat: A synthesis strategy for the Pauli graph.
    :param cx_config: A configuration of CXs to convert Pauli gadgets into.
    :return: a pass to perform the simplification
    """
def PauliSimp(strat: pytket._tket.transform.PauliSynthStrat = pytket._tket.transform.PauliSynthStrat.Sets, cx_config: pytket._tket.circuit.CXConfigType = pytket._tket.circuit.CXConfigType.Snake) -> BasePass:
    """
    Construct a pass that converts a circuit into a graph of Pauli gadgets to account for commutation and phase folding, and resynthesises them as either individual gagdets, pairwise constructions, or by diagonalising sets of commuting gadgets.
    
    This pass will not preserve the global phase of the circuit.
    
    :param strat: A synthesis strategy for the Pauli graph.
    :param cx_config: A configuration of CXs to convert Pauli gadgets into.
    :return: a pass to perform the simplification
    """
def PauliSquash(strat: pytket._tket.transform.PauliSynthStrat = pytket._tket.transform.PauliSynthStrat.Sets, cx_config: pytket._tket.circuit.CXConfigType = pytket._tket.circuit.CXConfigType.Snake) -> BasePass:
    """
    Applies :py:meth:`PauliSimp` followed by :py:meth:`FullPeepholeOptimise`.
    
    :param strat: a synthesis strategy for the Pauli graph
    :param cx_config: a configuration of CXs to convert Pauli gadgets into
    :return: a pass to perform the simplification
    """
def PeepholeOptimise2Q(allow_swaps: bool = True) -> BasePass:
    """
    Performs peephole optimisation including resynthesis of 2-qubit gate sequences, and converts to a circuit containing only CX and TK1 gates.
    
    :param allow_swaps: whether to allow implicit wire swaps
    """
def PlacementPass(placer: pytket._tket.placement.Placement) -> BasePass:
    """
    :param placer: The Placement used for relabelling.
    :return: a pass to relabel :py:class:`Circuit` Qubits to :py:class:`Architecture` Nodes
    """
@typing.overload
def RebaseCustom(gateset: set[pytket._tket.circuit.OpType], cx_replacement: pytket._tket.circuit.Circuit, tk1_replacement: typing.Callable[[typing.Union[sympy.Expr, float], typing.Union[sympy.Expr, float], typing.Union[sympy.Expr, float]], pytket._tket.circuit.Circuit]) -> BasePass:
    """
    Construct a custom rebase pass, given user-defined rebases for TK1 and CX. This pass:
    
    1. decomposes multi-qubit gates not in the set of gate types `gateset` to CX gates;
    2. if CX is not in `gateset`, replaces CX gates with `cx_replacement`;
    3. converts any single-qubit gates not in the gate type set to the form :math:`\mathrm{Rz}(a)\mathrm{Rx}(b)\mathrm{Rz}(c)` (in matrix-multiplication order, i.e. reverse order in the circuit);
    4. applies the `tk1_replacement` function to each of these triples :math:`(a,b,c)` to generate replacement circuits.
    
    :param gateset: the allowed operations in the rebased circuit
    :param cx_replacement: the equivalent circuit to replace a CX gate using two qubit gates from the desired basis (can use any single qubit OpTypes)
    :param tk1_replacement: a function which, given the parameters of an Rz(a)Rx(b)Rz(c) triple, returns an equivalent circuit in the desired basis
    :return: a pass that rebases to the given gate set (possibly including conditional and phase operations)
    """
@typing.overload
def RebaseCustom(gateset: set[pytket._tket.circuit.OpType], tk2_replacement: typing.Callable[[typing.Union[sympy.Expr, float], typing.Union[sympy.Expr, float], typing.Union[sympy.Expr, float]], pytket._tket.circuit.Circuit], tk1_replacement: typing.Callable[[typing.Union[sympy.Expr, float], typing.Union[sympy.Expr, float], typing.Union[sympy.Expr, float]], pytket._tket.circuit.Circuit]) -> BasePass:
    """
    Construct a custom rebase pass, given user-defined rebases for TK1 and TK2. This pass:
    
    1. decomposes multi-qubit gates not in the set of gate types `gateset` to TK2 gates;
    2. if TK2 is not in `gateset`, replaces TK2(a,b,c) gates via the `tk2_replacement` function;
    3. converts any single-qubit gates not in the gate type set to TK1;
    4. if TK2 is not in `gateset`. applies the `tk1_replacement` function to each TK1(a,b,c).
    
    :param gateset: the allowed operations in the rebased circuit
    :param tk2_replacement: a function which, given the parameters (a,b,c) of an XXPhase(a)YYPhase(b)ZZPhase(c) triple, returns an equivalent circuit in the desired basis
    :param tk1_replacement: a function which, given the parameters (a,b,c) of an Rz(a)Rx(b)Rz(c) triple, returns an equivalent circuit in the desired basis
    :return: a pass that rebases to the given gate set (possibly including conditional and phase operations)
    """
def RebaseTket() -> BasePass:
    """
    Converts all gates to CX, TK1 and Phase.
    """
def RemoveBarriers() -> BasePass:
    """
    A pass to remove all barrier instructions from the circuit.
    """
def RemoveDiscarded() -> BasePass:
    """
    A pass to remove all operations that have no ``OpType.Output`` or ``OpType.ClOutput`` in their causal future (in other words, all operations whose causal future is discarded).
    """
def RemoveImplicitQubitPermutation() -> BasePass:
    """
    Remove any implicit qubit permutation by appending SWAP gates.
    
    Note that if the circuit contains measurements, they may become mid-circuit measurements in the transformed circuit.
    """
def RemoveRedundancies() -> BasePass:
    """
    Removes gate-inverse pairs, merges rotations, removes identity rotations, and removes redundant gates before measurement. Does not add any new gate types.
    
    When merging rotations with the same op group name, the merged operation keeps the same name.
    """
def RenameQubitsPass(qubit_map: dict[pytket._tket.unit_id.Qubit, pytket._tket.unit_id.Qubit]) -> BasePass:
    """
    Rename some or all qubits. 
    
    :param qubit_map: map from old to new qubit names 
    """
def RoundAngles(n: int, only_zeros: bool = False) -> BasePass:
    """
    Round angles to the nearest :math:`\pi / 2^n`. 
    
    :param n: precision parameter, must be >= 0 and < 32 
    :param only_zeros: if True, only round angles less than :math:`\pi / 2^{n+1}` to zero, leave other angles alone (default False)
    """
def RoutingPass(arc: pytket._tket.architecture.Architecture) -> BasePass:
    """
    Construct a pass to route to the connectivity graph of an :py:class:`Architecture`. Edge direction is ignored. Uses :py:class:`LexiLabellingMethod` and :py:class:`LexiRouteRoutingMethod`.
    
    :return: a pass that routes to the given device architecture
    """
def SimplifyInitial(allow_classical: bool = True, create_all_qubits: bool = False, remove_redundancies: bool = True, xcirc: pytket._tket.circuit.Circuit | None = None) -> BasePass:
    """
    Simplify the circuit using knowledge of qubit state.
    
    :param allow_classical: allow replacement of measurements on known state with classical set-bit operations
    :param create_all_qubits: automatically annotate all qubits as initialized to the zero state
    :param remove_redundancies: apply a :py:meth:`RemoveRedundancies` pass after the initial simplification
    :param xcirc: 1-qubit circuit implementing an X gate in the transformed circuit (if omitted, an X gate is used)
    :return: a pass to perform the simplification
    """
def SimplifyMeasured() -> BasePass:
    """
    A pass to replace all 'classical maps' followed by measure operations whose quantum output is discarded with classical operations following the measure. (A 'classical map' is a quantum operation that acts as a permutation of the computational basis states followed by a diagonal operation.)
    """
def SquashCustom(singleqs: set[pytket._tket.circuit.OpType], tk1_replacement: typing.Callable[[typing.Union[sympy.Expr, float], typing.Union[sympy.Expr, float], typing.Union[sympy.Expr, float]], pytket._tket.circuit.Circuit], always_squash_symbols: bool = False) -> BasePass:
    """
    Squash sequences of single qubit gates from the target gate set into an optimal form given by `tk1_replacement`.
    
    :param singleqs: The types of single qubit gates in the target gate set. This pass will only affect sequences of gates that are already in this set.
    :param tk1_replacement: A function which, given the parameters of an Rz(a)Rx(b)Rz(c) triple, returns an equivalent circuit in the desired basis.
    :param always_squash_symbols: If true, always squash symbolic gates regardless of the blow-up in complexity. Default is false, meaning that symbolic gates are only squashed if doing so reduces the overall symbolic complexity.
    """
def SquashRzPhasedX() -> BasePass:
    """
    Squash single qubit gates into PhasedX and Rz gates. Commute Rz gates to the back if possible.
    """
def SquashTK1() -> BasePass:
    """
    Squash sequences of single-qubit gates to TK1 gates.
    """
def SynthesiseHQS() -> BasePass:
    """
    Optimises and converts a circuit consisting of CX and single-qubit gates into one containing only ZZMax, PhasedX, Rz and Phase.
    """
def SynthesiseOQC() -> BasePass:
    """
    Optimises and converts all gates to ECR, Rz, SX and Phase.
    """
def SynthesiseTK() -> BasePass:
    """
    Optimises and converts all gates to TK2, TK1 and Phase gates.
    """
def SynthesiseTket() -> BasePass:
    """
    Optimises and converts all gates to CX, TK1 and Phase gates.
    """
def SynthesiseUMD() -> BasePass:
    """
    Optimises and converts all gates to XXPhase, PhasedX, Rz and Phase.
    """
def ThreeQubitSquash(allow_swaps: bool = True) -> BasePass:
    """
    Squash three-qubit subcircuits into subcircuits having fewer CX gates, when possible, and apply Clifford simplification.
    
    :param allow_swaps: whether to allow implicit wire swaps
    """
def ZXGraphlikeOptimisation() -> BasePass:
    """
    Attempt to optimise the circuit by simplifying in ZX calculus and extracting a circuit back out. Due to limitations in extraction, may not work if the circuit contains created or discarded qubits. As a resynthesis pass, this will ignore almost all optimisations achieved beforehand and may increase the cost of the circuit.
    """
def ZZPhaseToRz() -> BasePass:
    """
    Converts all ZZPhase gates in a circuit with angle 1 or -1 (half-turns) into two Rz gates each with a parameter value of 1 (half-turns). ZZPhase gates with parameter values other than 1 or -1 (half-turns) are left unchanged.
    
    :return: a pass to convert ZZPhase gates to Rz.
    """
def _sympy_import() -> sympy.Expr | float:
    """
    This function only exists so that sympy gets imported in the resulting .pyi file. It's needed due to a bug in pybind11-stubgen's translation for Callables most likely.
    """
Audit: SafetyMode  # value = <SafetyMode.Audit: 0>
Default: SafetyMode  # value = <SafetyMode.Default: 1>
HamPath: CNotSynthType  # value = <CNotSynthType.HamPath: 1>
Rec: CNotSynthType  # value = <CNotSynthType.Rec: 2>
SWAP: CNotSynthType  # value = <CNotSynthType.SWAP: 0>
