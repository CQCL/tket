# Copyright 2019-2022 Cambridge Quantum Computing
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

from typing import List
from lark import Lark, Transformer
from pytket.circuit import OpType  # type: ignore
from pytket.passes import BasePass, RepeatPass, SequencePass  # type: ignore
from pytket.passes import (  # type: ignore
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
    PauliSimp,
    PauliSquash,
    PeepholeOptimise2Q,
    RebaseCirq,
    RebaseHQS,
    RebaseProjectQ,
    RebasePyZX,
    RebaseQuil,
    RebaseTket,
    RebaseUMD,
    RebaseUFR,
    RebaseOQC,
    RemoveBarriers,
    RemoveDiscarded,
    RemoveRedundancies,
    SimplifyInitial,
    SimplifyMeasured,
    SynthesiseHQS,
    SynthesiseTket,
    SynthesiseOQC,
    SynthesiseUMD,
    SquashHQS,
    ThreeQubitSquash,
)
from pytket.transform import CXConfigType, PauliSynthStrat  # type: ignore

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
    | rebase_cirq
    | rebase_hqs
    | rebase_oqc
    | rebase_projectq
    | rebase_pyzx
    | rebase_quil
    | rebase_tket
    | rebase_ufr
    | rebase_umd
    | remove_barriers
    | remove_discarded
    | remove_redundancies
    | simplify_initial
    | simplify_initial_no_classical
    | simplify_measured
    | squash_hqs
    | synthesise_hqs
    | synthesise_tket
    | synthesise_oqc
    | synthesise_umd
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
euler_angle_reduction: "EulerAngleReduction" "(" op_type "," op_type ")"
flatten_registers: "FlattenRegisters"
full_peephole_optimise: "FullPeepholeOptimise"
full_peephole_optimise_no_swaps: "FullPeepholeOptimiseNoSwaps"
guided_pauli_simp: "GuidedPauliSimp" "(" pauli_synth_strat "," cx_config_type ")"
guided_pauli_simp_default: "GuidedPauliSimp"
kak_decomposition: "KAKDecomposition"
optimise_phase_gadgets: "OptimisePhaseGadgets" "(" cx_config_type ")"
optimise_phase_gadgets_default: "OptimisePhaseGadgets"
pauli_simp: "PauliSimp" "(" pauli_synth_strat "," cx_config_type ")"
pauli_simp_default: "PauliSimp"
pauli_squash: "PauliSquash" "(" pauli_synth_strat "," cx_config_type ")"
pauli_squash_default: "PauliSquash"
peephole_optimise_2q: "PeepholeOptimise2Q"
rebase_cirq: "RebaseCirq"
rebase_hqs: "RebaseHQS"
rebase_oqc: "RebaseOQC"
rebase_projectq: "RebaseProjectQ"
rebase_pyzx: "RebasePyZX"
rebase_quil: "RebaseQuil"
rebase_tket: "RebaseTket"
rebase_ufr: "RebaseUFR"
rebase_umd: "RebaseUMD"
remove_barriers: "RemoveBarriers"
remove_discarded: "RemoveDiscarded"
remove_redundancies: "RemoveRedundancies"
simplify_initial: "SimplifyInitial"
simplify_initial_no_classical: "SimplifyInitialNoClassical"
simplify_measured: "SimplifyMeasured"
squash_hqs: "SquashHQS"
synthesise_hqs: "SynthesiseHQS"
synthesise_tket: "SynthesiseTket"
synthesise_oqc: "SynthesiseOQC"
synthesise_umd: "SynthesiseUMD"
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


class PassTransformer(Transformer):
    def start(self, t: List) -> BasePass:
        return t[0]

    def comp_pass(self, t: List) -> BasePass:
        return t[0]

    def basic_pass(self, t: list) -> BasePass:
        return t[0]

    def seq_pass(self, t: List) -> BasePass:
        return t[0]

    def pass_list(self, t: List) -> BasePass:
        return SequencePass(t)

    def repeat_pass(self, t: List) -> BasePass:
        return RepeatPass(t[0])

    def clifford_simp(self, t: List) -> BasePass:
        return CliffordSimp()

    def clifford_simp_no_swaps(self, t: List) -> BasePass:
        return CliffordSimp(allow_swaps=False)

    def commute_through_multis(self, t: List) -> BasePass:
        return CommuteThroughMultis()

    def context_simp(self, t: List) -> BasePass:
        return ContextSimp()

    def context_simp_no_classical(self, t: List) -> BasePass:
        return ContextSimp(allow_classical=False)

    def decompose_arbitrarily_controlled_gates(self, t: List) -> BasePass:
        return DecomposeArbitrarilyControlledGates()

    def decompose_boxes(self, t: List) -> BasePass:
        return DecomposeBoxes()

    def decompose_classical_exp(self, t: List) -> BasePass:
        return DecomposeClassicalExp()

    def decompose_multi_qubits_cx(self, t: List) -> BasePass:
        return DecomposeMultiQubitsCX()

    def decompose_single_qubits_tk1(self, t: List) -> BasePass:
        return DecomposeSingleQubitsTK1()

    def delay_measures(self, t: List) -> BasePass:
        return DelayMeasures()

    def euler_angle_reduction(self, t: List) -> BasePass:
        return EulerAngleReduction(t[0], t[1])

    def flatten_registers(self, t: List) -> BasePass:
        return FlattenRegisters()

    def full_peephole_optimise(self, t: List) -> BasePass:
        return FullPeepholeOptimise()

    def full_peephole_optimise_no_swaps(self, t: List) -> BasePass:
        return FullPeepholeOptimise(allow_swaps=False)

    def guided_pauli_simp(self, t: List) -> BasePass:
        return GuidedPauliSimp(strat=t[0], cx_config=t[1])

    def guided_pauli_simp_default(self, t: List) -> BasePass:
        return GuidedPauliSimp()

    def kak_decomposition(self, t: List) -> BasePass:
        return KAKDecomposition()

    def optimise_phase_gadgets(self, t: List) -> BasePass:
        return OptimisePhaseGadgets(cx_config=t[0])

    def optimise_phase_gadgets_default(self, t: List) -> BasePass:
        return OptimisePhaseGadgets()

    def pauli_simp(self, t: List) -> BasePass:
        return PauliSimp(strat=t[0], cx_config=t[1])

    def pauli_simp_default(self, t: List) -> BasePass:
        return PauliSimp()

    def pauli_squash(self, t: List) -> BasePass:
        return PauliSquash(strat=t[0], cx_config=t[1])

    def pauli_squash_default(self, t: List) -> BasePass:
        return PauliSquash()

    def peephole_optimise_2q(self, t: List) -> BasePass:
        return PeepholeOptimise2Q()

    def rebase_cirq(self, t: List) -> BasePass:
        return RebaseCirq()

    def rebase_hqs(self, t: List) -> BasePass:
        return RebaseHQS()

    def rebase_oqc(self, t: List) -> BasePass:
        return RebaseOQC()

    def rebase_projectq(self, t: List) -> BasePass:
        return RebaseProjectQ()

    def rebase_pyzx(self, t: List) -> BasePass:
        return RebasePyZX()

    def rebase_quil(self, t: List) -> BasePass:
        return RebaseQuil()

    def rebase_tket(self, t: List) -> BasePass:
        return RebaseTket()

    def rebase_ufr(self, t: List) -> BasePass:
        return RebaseUFR()

    def rebase_umd(self, t: List) -> BasePass:
        return RebaseUMD()

    def remove_barriers(self, t: List) -> BasePass:
        return RemoveBarriers()

    def remove_discarded(self, t: List) -> BasePass:
        return RemoveDiscarded()

    def remove_redundancies(self, t: List) -> BasePass:
        return RemoveRedundancies()

    def simplify_initial(self, t: List) -> BasePass:
        return SimplifyInitial()

    def simplify_initial_no_classical(self, t: List) -> BasePass:
        return SimplifyInitial(allow_classical=False)

    def simplify_measured(self, t: List) -> BasePass:
        return SimplifyMeasured()

    def squash_hqs(self, t: List) -> BasePass:
        return SquashHQS()

    def synthesise_hqs(self, t: List) -> BasePass:
        return SynthesiseHQS()

    def synthesise_tket(self, t: List) -> BasePass:
        return SynthesiseTket()

    def synthesise_oqc(self, t: List) -> BasePass:
        return SynthesiseOQC()

    def synthesise_umd(self, t: List) -> BasePass:
        return SynthesiseUMD()

    def three_qubit_squash(self, t: List) -> BasePass:
        return ThreeQubitSquash()

    def cx_config_type(self, t: List) -> CXConfigType:
        return t[0]

    def cx_config_type_snake(self, t: List) -> CXConfigType:
        return CXConfigType.Snake

    def cx_config_type_star(self, t: List) -> CXConfigType:
        return CXConfigType.Star

    def cx_config_type_tree(self, t: List) -> CXConfigType:
        return CXConfigType.Tree

    def cx_config_type_multi_q_gate(self, t: List) -> CXConfigType:
        return CXConfigType.MultiQGate

    def op_type(self, t: List) -> OpType:
        return t[0]

    def op_type_rx(self, t: List) -> OpType:
        return OpType.Rx

    def op_type_ry(self, t: List) -> OpType:
        return OpType.Ry

    def op_type_rz(self, t: List) -> OpType:
        return OpType.Rz

    def pauli_synth_strat(self, t: List) -> PauliSynthStrat:
        return t[0]

    def pauli_synth_strat_individual(self, t: List) -> PauliSynthStrat:
        return PauliSynthStrat.Individual

    def pauli_synth_strat_pairwise(self, t: List) -> PauliSynthStrat:
        return PauliSynthStrat.Pairwise

    def pauli_synth_strat_sets(self, t: List) -> PauliSynthStrat:
        return PauliSynthStrat.Sets


parser = Lark(pass_grammar)
transformer = PassTransformer()


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
    return transformer.transform(tree)


def compilation_pass_grammar() -> str:
    """Formal grammar for specifying compilation passes.

    This is the grammar assumed by :py:meth:`complilation_pass_from_script`.

    :return: grammar in extended Backus--Naur form"""
    return pass_grammar
