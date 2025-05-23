# Copyright Quantinuum
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from typing import cast

from lark import Lark, Transformer

from pytket.circuit import CXConfigType, OpType
from pytket.passes import (
    BasePass,
    CliffordSimp,
    CommuteThroughMultis,
    ContextSimp,
    DecomposeArbitrarilyControlledGates,
    DecomposeBoxes,
    DecomposeClassicalExp,
    DecomposeMultiQubitsCX,
    DecomposeSingleQubitsTK1,
    DelayMeasures,
    EulerAngleReduction,
    FlattenRegisters,
    FullPeepholeOptimise,
    GuidedPauliSimp,
    KAKDecomposition,
    OptimisePhaseGadgets,
    PauliExponentials,
    PauliSimp,
    PauliSquash,
    PeepholeOptimise2Q,
    RebaseTket,
    RemoveBarriers,
    RemoveDiscarded,
    RemovePhaseOps,
    RemoveRedundancies,
    RepeatPass,
    SequencePass,
    SimplifyInitial,
    SimplifyMeasured,
    SynthesiseTket,
    ThreeQubitSquash,
)
from pytket.transform import PauliSynthStrat

pass_grammar = """
start: comp_pass
comp_pass: ( basic_pass | seq_pass | repeat_pass )
basic_pass:
    | clifford_simp
    | clifford_simp_no_swaps
    | commute_through_multis
    | context_simp
    | context_simp_no_classical
    | decompose_arbitrarily_controlled_gates
    | decompose_boxes
    | decompose_classical_exp
    | decompose_multi_qubits_cx
    | decompose_single_qubits_tk1
    | delay_measures
    | delay_measures_try
    | euler_angle_reduction
    | flatten_registers
    | full_peephole_optimise
    | full_peephole_optimise_no_swaps
    | guided_pauli_simp
    | guided_pauli_simp_default
    | kak_decomposition
    | optimise_phase_gadgets
    | optimise_phase_gadgets_default
    | pauli_simp
    | pauli_simp_default
    | pauli_squash
    | pauli_squash_default
    | peephole_optimise_2q
    | rebase_tket
    | remove_barriers
    | remove_discarded
    | remove_phase_ops
    | remove_redundancies
    | simplify_initial
    | simplify_initial_no_classical
    | simplify_measured
    | synthesise_tket
    | three_qubit_squash
seq_pass: "[" pass_list "]"
pass_list: comp_pass ("," comp_pass)*
repeat_pass: "repeat" "(" comp_pass ")"

clifford_simp: "CliffordSimp"
clifford_simp_no_swaps: "CliffordSimpNoSwaps"
commute_through_multis: "CommuteThroughMultis"
context_simp: "ContextSimp"
context_simp_no_classical: "ContextSimpNoClassical"
decompose_arbitrarily_controlled_gates: "DecomposeArbitrarilyControlledGates"
decompose_boxes: "DecomposeBoxes"
decompose_classical_exp: "DecomposeClassicalExp"
decompose_multi_qubits_cx: "DecomposeMultiQubitsCX"
decompose_single_qubits_tk1: "DecomposeSingleQubitsTK1"
delay_measures: "DelayMeasures"
delay_measures_try: "TryDelayMeasures"
euler_angle_reduction: "EulerAngleReduction" "(" op_type "," op_type ")"
flatten_registers: "FlattenRegisters"
full_peephole_optimise: "FullPeepholeOptimise"
full_peephole_optimise_no_swaps: "FullPeepholeOptimiseNoSwaps"
guided_pauli_simp: "GuidedPauliSimp" "(" pauli_synth_strat "," cx_config_type ")"
guided_pauli_simp_default: "GuidedPauliSimp"
kak_decomposition: "KAKDecomposition"
optimise_phase_gadgets: "OptimisePhaseGadgets" "(" cx_config_type ")"
optimise_phase_gadgets_default: "OptimisePhaseGadgets"
pauli_exponentials: "PauliExponentials" "(" pauli_synth_strat "," cx_config_type ")"
pauli_exponentials_default: "PauliExponentials"
pauli_simp: "PauliSimp" "(" pauli_synth_strat "," cx_config_type ")"
pauli_simp_default: "PauliSimp"
pauli_squash: "PauliSquash" "(" pauli_synth_strat "," cx_config_type ")"
pauli_squash_default: "PauliSquash"
peephole_optimise_2q: "PeepholeOptimise2Q"
rebase_tket: "RebaseTket"
remove_barriers: "RemoveBarriers"
remove_discarded: "RemoveDiscarded"
remove_phase_ops: "RemovePhaseOps"
remove_redundancies: "RemoveRedundancies"
simplify_initial: "SimplifyInitial"
simplify_initial_no_classical: "SimplifyInitialNoClassical"
simplify_measured: "SimplifyMeasured"
synthesise_tket: "SynthesiseTket"
three_qubit_squash: "ThreeQubitSquash"

cx_config_type:
    | cx_config_type_snake
    | cx_config_type_star
    | cx_config_type_tree
    | cx_config_type_multi_q_gate
cx_config_type_snake: "Snake"
cx_config_type_star: "Star"
cx_config_type_tree: "Tree"
cx_config_type_multi_q_gate: "MultiQGate"
op_type: ( op_type_rx | op_type_ry | op_type_rz )
op_type_rx: "Rx"
op_type_ry: "Ry"
op_type_rz: "Rz"
pauli_synth_strat:
    | pauli_synth_strat_individual
    | pauli_synth_strat_pairwise
    | pauli_synth_strat_sets
pauli_synth_strat_individual: "Individual"
pauli_synth_strat_pairwise: "Pairwise"
pauli_synth_strat_sets: "Sets"

%import common.WS_INLINE -> WS
%import common.CR
%import common.LF
_NEWLINE: CR? LF
%ignore WS
%ignore _NEWLINE
"""


class _PassTransformer(Transformer):
    def start(self, t: list[BasePass]) -> BasePass:
        return t[0]

    def comp_pass(self, t: list[BasePass]) -> BasePass:
        return t[0]

    def basic_pass(self, t: list[BasePass]) -> BasePass:
        return t[0]

    def seq_pass(self, t: list[BasePass]) -> BasePass:
        return t[0]

    def pass_list(self, t: list[BasePass]) -> BasePass:
        return SequencePass(t)

    def repeat_pass(self, t: list[BasePass]) -> BasePass:
        return RepeatPass(t[0])

    def clifford_simp(self, t: list) -> BasePass:
        return CliffordSimp()

    def clifford_simp_no_swaps(self, t: list) -> BasePass:
        return CliffordSimp(allow_swaps=False)

    def commute_through_multis(self, t: list) -> BasePass:
        return CommuteThroughMultis()

    def context_simp(self, t: list) -> BasePass:
        return ContextSimp()

    def context_simp_no_classical(self, t: list) -> BasePass:
        return ContextSimp(allow_classical=False)

    def decompose_arbitrarily_controlled_gates(self, t: list) -> BasePass:
        return DecomposeArbitrarilyControlledGates()

    def decompose_boxes(self, t: list) -> BasePass:
        return DecomposeBoxes()

    def decompose_classical_exp(self, t: list) -> BasePass:
        return DecomposeClassicalExp()

    def decompose_multi_qubits_cx(self, t: list) -> BasePass:
        return DecomposeMultiQubitsCX()

    def decompose_single_qubits_tk1(self, t: list) -> BasePass:
        return DecomposeSingleQubitsTK1()

    def delay_measures(self, t: list) -> BasePass:
        return DelayMeasures(False)

    def delay_measures_try(self, t: list) -> BasePass:
        return DelayMeasures(True)

    def euler_angle_reduction(self, t: list[OpType]) -> BasePass:
        return EulerAngleReduction(t[0], t[1])

    def flatten_registers(self, t: list) -> BasePass:
        return FlattenRegisters()

    def full_peephole_optimise(self, t: list) -> BasePass:
        return FullPeepholeOptimise()

    def full_peephole_optimise_no_swaps(self, t: list) -> BasePass:
        return FullPeepholeOptimise(allow_swaps=False)

    def guided_pauli_simp(self, t: list) -> BasePass:
        assert isinstance(t[0], PauliSynthStrat)
        assert isinstance(t[1], CXConfigType)
        return GuidedPauliSimp(strat=t[0], cx_config=t[1])

    def guided_pauli_simp_default(self, t: list) -> BasePass:
        return GuidedPauliSimp()

    def kak_decomposition(self, t: list) -> BasePass:
        return KAKDecomposition()

    def optimise_phase_gadgets(self, t: list) -> BasePass:
        assert isinstance(t[0], CXConfigType)
        return OptimisePhaseGadgets(cx_config=t[0])

    def optimise_phase_gadgets_default(self, t: list) -> BasePass:
        return OptimisePhaseGadgets()

    def pauli_exponentials(self, t: list) -> BasePass:
        assert isinstance(t[0], PauliSynthStrat)
        assert isinstance(t[1], CXConfigType)
        return PauliExponentials(strat=t[0], cx_config=t[1])

    def pauli_exponentials_default(self, t: list) -> BasePass:
        return PauliExponentials()

    def pauli_simp(self, t: list) -> BasePass:
        assert isinstance(t[0], PauliSynthStrat)
        assert isinstance(t[1], CXConfigType)
        return PauliSimp(strat=t[0], cx_config=t[1])

    def pauli_simp_default(self, t: list) -> BasePass:
        return PauliSimp()

    def pauli_squash(self, t: list) -> BasePass:
        assert isinstance(t[0], PauliSynthStrat)
        assert isinstance(t[1], CXConfigType)
        return PauliSquash(strat=t[0], cx_config=t[1])

    def pauli_squash_default(self, t: list) -> BasePass:
        return PauliSquash()

    def peephole_optimise_2q(self, t: list) -> BasePass:
        return PeepholeOptimise2Q()

    def rebase_tket(self, t: list) -> BasePass:
        return RebaseTket()

    def remove_barriers(self, t: list) -> BasePass:
        return RemoveBarriers()

    def remove_discarded(self, t: list) -> BasePass:
        return RemoveDiscarded()

    def remove_phase_ops(self, t: list) -> BasePass:
        return RemovePhaseOps()

    def remove_redundancies(self, t: list) -> BasePass:
        return RemoveRedundancies()

    def simplify_initial(self, t: list) -> BasePass:
        return SimplifyInitial()

    def simplify_initial_no_classical(self, t: list) -> BasePass:
        return SimplifyInitial(allow_classical=False)

    def simplify_measured(self, t: list) -> BasePass:
        return SimplifyMeasured()

    def synthesise_tket(self, t: list) -> BasePass:
        return SynthesiseTket()

    def three_qubit_squash(self, t: list) -> BasePass:
        return ThreeQubitSquash()

    def cx_config_type(self, t: list[CXConfigType]) -> CXConfigType:
        return t[0]

    def cx_config_type_snake(self, t: list) -> CXConfigType:
        return CXConfigType.Snake

    def cx_config_type_star(self, t: list) -> CXConfigType:
        return CXConfigType.Star

    def cx_config_type_tree(self, t: list) -> CXConfigType:
        return CXConfigType.Tree

    def cx_config_type_multi_q_gate(self, t: list) -> CXConfigType:
        return CXConfigType.MultiQGate

    def op_type(self, t: list[OpType]) -> OpType:
        return t[0]

    def op_type_rx(self, t: list) -> OpType:
        return OpType.Rx

    def op_type_ry(self, t: list) -> OpType:
        return OpType.Ry

    def op_type_rz(self, t: list) -> OpType:
        return OpType.Rz

    def pauli_synth_strat(self, t: list[PauliSynthStrat]) -> PauliSynthStrat:
        return t[0]

    def pauli_synth_strat_individual(self, t: list) -> PauliSynthStrat:
        return PauliSynthStrat.Individual

    def pauli_synth_strat_pairwise(self, t: list) -> PauliSynthStrat:
        return PauliSynthStrat.Pairwise

    def pauli_synth_strat_sets(self, t: list) -> PauliSynthStrat:
        return PauliSynthStrat.Sets


parser = Lark(pass_grammar)
transformer = _PassTransformer()


def compilation_pass_from_script(script: str) -> BasePass:
    """Generate a compilation pass from a specification.

    The specification must conform to a simple grammar. For example, the following are
    valid specifications:

    * "RemoveRedundancies"
    * "[RemoveBarriers, RemoveRedundancies]" (a sequence of passes)
    * "repeat(FullPeepholeOptimise)" (repeat a pass until it doesn't change the circuit)

    Sequences and repeats can be nested arbitrarily. Whitespace is ignored.

    Most passes are specified using their Python names. For those that take enums as
    parameters, non-default values can be specified using their Python names:

    * "PauliSimp" (default parameters)
    * "PauliSimp(Pairwise, Tree)"
    * "EulerAngleReduction(Ry, Rz)"

    For some passes with optional boolean parameters the name can be modified as
    follows:

    * "CliffordSimp" (default parameters)
    * "CliffordSimpNoSwaps"
    * "SimplifyInitial" (default parameters)
    * "SimplifyInitialNoClassical"

    There is currently no support for passes requiring more complex parameters such as
    lambdas or circuits.

    The full formal grammar can be inspected using :py:meth:`compilation_pass_grammar`.

    :param script: specification of pass
    """
    tree = parser.parse(script)
    return cast("BasePass", transformer.transform(tree))


def compilation_pass_grammar() -> str:
    """Formal grammar for specifying compilation passes.

    This is the grammar assumed by :py:meth:`compilation_pass_from_script`.

    :return: grammar in extended Backus--Naur form"""
    return pass_grammar
