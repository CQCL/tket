# Copyright 2019-2023 Cambridge Quantum Computing
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

from pytket.circuit import Circuit, OpType, Op, PauliExpBox, Unitary2qBox, Node, Qubit  # type: ignore
from pytket.pauli import Pauli  # type: ignore
from pytket.passes import (  # type: ignore
    SequencePass,
    RemoveRedundancies,
    SynthesiseTket,
    SynthesiseHQS,
    SynthesiseUMD,
    RepeatUntilSatisfiedPass,
    CommuteThroughMultis,
    RepeatPass,
    DecomposeMultiQubitsCX,
    SquashTK1,
    SquashRzPhasedX,
    RepeatWithMetricPass,
    RebaseCustom,
    EulerAngleReduction,
    RoutingPass,
    CXMappingPass,
    PlacementPass,
    NaivePlacementPass,
    RenameQubitsPass,
    FullMappingPass,
    DefaultMappingPass,
    AASRouting,
    DecomposeSwapsToCXs,
    DecomposeSwapsToCircuit,
    PauliSimp,
    ThreeQubitSquash,
    RebaseTket,
    RemoveDiscarded,
    SimplifyMeasured,
    SimplifyInitial,
    RemoveBarriers,
    PauliSquash,
    auto_rebase_pass,
    ZZPhaseToRz,
    CnXPairwiseDecomposition,
)
from pytket.predicates import (  # type: ignore
    GateSetPredicate,
    NoClassicalControlPredicate,
    DirectednessPredicate,
    NoBarriersPredicate,
    CompilationUnit,
)
from pytket.mapping import (  # type: ignore
    LexiLabellingMethod,
    LexiRouteRoutingMethod,
)
from pytket.architecture import Architecture  # type: ignore
from pytket.placement import Placement, GraphPlacement  # type: ignore
from pytket.transform import Transform, PauliSynthStrat, CXConfigType  # type: ignore
from pytket._tket.passes import SynthesiseOQC  # type: ignore
import numpy as np

from pytket.qasm import circuit_to_qasm_str
import pytest  # type: ignore
from sympy import Symbol, Expr  # type: ignore
from typing import Dict, Any, List, Union

Param = Union[float, "Expr"]

circ2 = Circuit(1)
circ2.Rx(0.25, 0)

ots = {OpType.X, OpType.Z}
gsp = GateSetPredicate(ots)
nccp = NoClassicalControlPredicate()


def tk1_to_phasedxrz(a: float, b: float, c: float) -> Circuit:
    circ = Circuit(1)
    circ.Rz(a + c, 0)
    phasedx_op = Op.create(OpType.PhasedX, [b, a])
    circ.add_gate(phasedx_op, [0])
    Transform.RemoveRedundancies().apply(circ)
    return circ


def test_predicate_generation() -> None:
    string_gsp = repr(gsp)
    assert string_gsp.find("X") != -1
    assert string_gsp.find("Z") != -1
    assert repr(nccp) == "NoClassicalControlPredicate"


def test_compilation_unit_generation() -> None:
    pp_list = [gsp, nccp]
    circ = Circuit(2)
    circ.X(0).Z(1)
    cu = CompilationUnit(circ, pp_list)
    assert cu.check_all_predicates()
    cu2 = CompilationUnit(circ2, pp_list)
    assert not cu2.check_all_predicates()


def test_compilerpass_seq() -> None:
    passlist = [SynthesiseTket(), SynthesiseOQC(), SynthesiseUMD(), SynthesiseHQS()]
    seq = SequencePass(passlist)
    circ = Circuit(2)
    circ.X(0).Z(1)
    cu = CompilationUnit(circ)
    cu2 = CompilationUnit(circ2)
    assert seq.apply(cu)
    assert seq.apply(cu2)


def test_rebase_pass_generation() -> None:
    cx = Circuit(2)
    cx.CX(0, 1)
    pz_rebase = RebaseCustom(
        {OpType.CX, OpType.PhasedX, OpType.Rz}, cx, tk1_to_phasedxrz
    )
    circ = Circuit(2)
    circ.X(0).Y(1)
    cu = CompilationUnit(circ)
    assert pz_rebase.apply(cu)
    coms = cu.circuit.get_commands()
    assert str(coms) == "[PhasedX(1, 0) q[0];, PhasedX(1, 0.5) q[1];]"

    passlist = [pz_rebase, SynthesiseTket()]
    seq = SequencePass(passlist)
    assert seq.apply(cu)
    coms = cu.circuit.get_commands()
    assert str(coms) == "[TK1(0, 1, 0) q[0];, TK1(0, 1, 3) q[1];]"


def test_rebase_pass_generation_via_TK2() -> None:
    def tk1(a: Param, b: Param, c: Param) -> "Circuit":
        circ = Circuit(1)
        circ.Rz(c, 0).Rx(b, 0).Rz(a, 0)
        return circ

    def tk2(a: Param, b: Param, c: Param) -> "Circuit":
        circ = Circuit(2)
        circ.add_gate(OpType.ZZPhase, c, [0, 1])
        circ.add_gate(OpType.YYPhase, b, [0, 1])
        circ.add_gate(OpType.XXPhase, a, [0, 1])
        return circ

    rebase_pass = RebaseCustom(
        {
            OpType.Rx,
            OpType.Ry,
            OpType.Rz,
            OpType.XXPhase,
            OpType.YYPhase,
            OpType.ZZPhase,
        },
        tk2,
        tk1,
    )

    circ = Circuit(3).H(0).CX(0, 1).H(1).CX(1, 2)
    rebase_pass.apply(circ)
    assert circ.n_gates_of_type(OpType.XXPhase) == 2
    assert circ.n_gates_of_type(OpType.YYPhase) == 0
    assert circ.n_gates_of_type(OpType.ZZPhase) == 0


def test_custom_combinator_generation() -> None:
    def test_CX_size_threshold(circ: Circuit) -> bool:
        return bool(circ.n_gates_of_type(OpType.CX) == 0)

    seq_pass = SequencePass([RemoveRedundancies(), CommuteThroughMultis()])
    custom_pass = RepeatUntilSatisfiedPass(seq_pass, test_CX_size_threshold)

    circ = Circuit(2)
    circ.CX(0, 1)
    circ.X(1)
    circ.CX(0, 1)
    circ.X(1)
    circ.CX(0, 1)
    circ.X(1)
    circ.CX(0, 1)
    circ.Z(1)
    circ.CX(1, 0)
    circ.Z(1)
    circ.CX(1, 0)

    cu = CompilationUnit(circ)
    assert custom_pass.apply(cu)

    # Test in-place application
    circ1 = cu.circuit
    assert custom_pass.apply(circ)
    assert circ == circ1
    assert not custom_pass.apply(circ)


def test_routing_and_placement_pass() -> None:
    # Qubit interaction graph:
    # 1 -- 0 -- 3
    #  \   |
    #   `  4 -- 2
    circ = Circuit()
    q = circ.add_q_register("q", 5)
    circ.CX(0, 1)
    circ.H(0)
    circ.Z(1)
    circ.CX(0, 3)
    circ.Rx(1.5, 3)
    circ.CX(2, 4)
    circ.X(2)
    circ.CX(1, 4)
    circ.CX(0, 4)
    # Architecture graph:
    # f2
    # |
    # b1 -- b0
    # |
    # b2 -- a0 -- f0
    n0 = Node("b", 0)
    n1 = Node("b", 1)
    n2 = Node("b", 2)
    n3 = Node("a", 0)
    n4 = Node("f", 0)
    n5 = Node("f", 2)
    arc = Architecture([[n0, n1], [n1, n2], [n2, n3], [n3, n4], [n1, n5]])
    pl = Placement(arc)
    routing = RoutingPass(arc)
    placement = PlacementPass(pl)
    nplacement = NaivePlacementPass(arc)
    cu = CompilationUnit(circ.copy())
    assert placement.apply(cu)
    assert routing.apply(cu)
    assert nplacement.apply(cu)
    arcnodes = arc.nodes
    expected_map = {
        q[0]: arcnodes[0],
        q[1]: arcnodes[1],
        q[2]: arcnodes[2],
        q[3]: arcnodes[3],
        q[4]: arcnodes[4],
    }
    assert cu.initial_map == expected_map

    # check composition works ok
    seq_pass = SequencePass([SynthesiseTket(), placement, routing, SynthesiseUMD()])
    cu2 = CompilationUnit(circ.copy())
    assert seq_pass.apply(cu2)
    assert cu2.initial_map == expected_map

    full_pass = FullMappingPass(
        arc, pl, config=[LexiLabellingMethod(), LexiRouteRoutingMethod()]
    )
    cu3 = CompilationUnit(circ.copy())
    assert full_pass.apply(cu3)
    assert cu3.initial_map == expected_map
    assert cu.circuit == cu3.circuit


def test_default_mapping_pass() -> None:
    circ = Circuit()
    q = circ.add_q_register("q", 6)
    circ.CX(0, 1)
    circ.H(0)
    circ.Z(1)
    circ.CX(0, 3)
    circ.Rx(1.5, 3)
    circ.CX(2, 4)
    circ.X(2)
    circ.CX(1, 4)
    circ.CX(0, 4)
    circ.H(5)
    n0 = Node("b", 0)
    n1 = Node("b", 1)
    n2 = Node("b", 2)
    n3 = Node("a", 0)
    n4 = Node("f", 0)
    n5 = Node("g", 7)
    arc = Architecture([[n0, n1], [n1, n2], [n2, n3], [n3, n4], [n4, n5]])
    pl = GraphPlacement(arc)

    nplacement = NaivePlacementPass(arc)
    routing = RoutingPass(arc)
    placement = PlacementPass(pl)
    default = DefaultMappingPass(arc)
    cu_rp = CompilationUnit(circ.copy())
    cu_def = CompilationUnit(circ.copy())

    assert placement.apply(cu_rp)
    assert routing.apply(cu_rp)
    assert nplacement.apply(cu_rp)
    assert default.apply(cu_def)
    assert cu_rp.circuit == cu_def.circuit


def test_default_mapping_pass_phase_poly_aas() -> None:
    circ = Circuit()
    q = circ.add_q_register("q", 5)
    circ.CX(0, 1)
    circ.H(0)
    circ.Z(1)
    circ.CX(0, 3)
    circ.Rx(1.5, 3)
    circ.CX(2, 4)
    circ.X(2)
    circ.CX(1, 4)
    circ.CX(0, 4)
    n0 = Node("a", 0)
    n1 = Node("b", 1)
    n2 = Node("c", 2)
    n3 = Node("d", 3)
    n4 = Node("e", 4)
    arc = Architecture([[n0, n1], [n1, n2], [n2, n3], [n3, n4]])
    default = AASRouting(arc, lookahead=1)
    assert default.apply(circ)


def test_rename_qubits_pass() -> None:
    circ = Circuit()
    qbs = circ.add_q_register("a", 2)
    circ.CX(Qubit("a", 0), Qubit("a", 1))
    qm = {Qubit("a", 0): Qubit("b", 1), Qubit("a", 1): Qubit("b", 0)}
    p = RenameQubitsPass(qm)
    cu = CompilationUnit(circ)
    p.apply(cu)
    newcirc = cu.circuit
    assert set(newcirc.qubits) == set([Qubit("b", i) for i in range(2)])


def gate_count_metric(circ: Circuit) -> int:
    return int(circ.n_gates)


def test_RebaseOQC_and_SynthesiseOQC() -> None:
    oqc_gateset = {OpType.SX, OpType.Rz, OpType.ECR}
    oqc_gateset_pred = GateSetPredicate(oqc_gateset)
    circ = Circuit(3)
    circ.CX(0, 1)
    circ.H(0)
    circ.Z(1)
    circ.CX(0, 2)
    circ.Rx(1.5, 2)
    circ.CX(2, 1)
    circ.X(2)
    circ.CX(1, 0)
    circ.CX(0, 1)
    u = circ.get_unitary()
    # Test SynthesiseOQC
    circ2 = circ.copy()
    SynthesiseOQC().apply(circ2)
    assert oqc_gateset_pred.verify(circ2)

    u_with_oqc = circ2.get_unitary()
    assert np.allclose(u, u_with_oqc)

    RebaseTket().apply(circ2)
    u2 = circ2.get_unitary()
    assert np.allclose(u, u2)

    # Test RebaseOQC
    circ3 = circ.copy()
    u_before_oqc = circ3.get_unitary()
    assert np.allclose(u, u_before_oqc)

    auto_rebase_pass(oqc_gateset).apply(circ3)
    assert oqc_gateset_pred.verify(circ3)
    u_before_rebase_tket = circ3.get_unitary()
    assert np.allclose(u, u_before_rebase_tket)

    RebaseTket().apply(circ3)
    u3 = circ3.get_unitary()
    assert np.allclose(u, u3)


def test_SynthesiseTket_creation() -> None:
    # my_synthesise_tket should act on a CompilationUnit the same as SynthesiseTket
    seq_pass = SequencePass([CommuteThroughMultis(), RemoveRedundancies()])
    repeat_pass = RepeatPass(seq_pass)
    synth_pass = SequencePass(
        [DecomposeMultiQubitsCX(), RemoveRedundancies(), repeat_pass, SquashTK1()]
    )
    small_part = SequencePass([RemoveRedundancies(), repeat_pass, SquashTK1()])
    repeat_synth_pass = RepeatWithMetricPass(small_part, gate_count_metric)
    my_synthesise_tket = SequencePass([synth_pass, repeat_synth_pass])

    circ1 = Circuit(3)
    circ1.X(0).Y(1).CX(0, 1).Z(0).Rx(1.3, 1).CX(0, 1).Rz(0.4, 0).Ry(0.53, 0).H(1).H(
        2
    ).Rx(1.5, 2).Rx(0.5, 2).H(2)
    cu1 = CompilationUnit(circ1)
    my_synthesise_tket.apply(cu1)
    circ2 = cu1.circuit
    assert circ2.n_gates == 2

    cu2 = CompilationUnit(circ1)
    # Blue Peter voice: here's one I made earlier
    SynthesiseTket().apply(cu2)
    circ3 = cu2.circuit
    assert circ3.n_gates == 2
    assert circ2 == circ3

    # now let's run with routing
    arc = Architecture([(0, 2), (1, 2)])
    pl = Placement(arc)
    pl.place(circ1)
    routing_pass = RoutingPass(arc)
    cu3 = CompilationUnit(circ1)
    cu4 = CompilationUnit(circ1)
    cu5 = CompilationUnit(circ1)
    assert routing_pass.apply(cu3)

    full_pass = SequencePass([SynthesiseTket(), routing_pass])
    full_pass2 = SequencePass([my_synthesise_tket, routing_pass])
    assert full_pass.apply(cu4)
    assert full_pass2.apply(cu5)
    assert cu4.circuit == cu5.circuit


def test_directed_cx_pass() -> None:
    circ = Circuit(5)
    circ.CX(0, 1)
    circ.Rx(1.4, 0)
    circ.H(1)
    circ.CX(0, 2)
    circ.CX(0, 3)
    circ.Sdg(4)
    circ.CX(2, 3)
    circ.CX(3, 4)
    circ.CX(4, 3)
    circ.CX(3, 1)
    arc = Architecture([[0, 1], [1, 2], [2, 3], [3, 4]])
    pl = Placement(arc)
    cu1 = CompilationUnit(circ)
    dir_router = CXMappingPass(arc, pl, directed_cx=True)
    assert dir_router.apply(cu1)

    circ2 = cu1.circuit
    dir_pred = DirectednessPredicate(arc)
    assert dir_pred.verify(circ2)


def test_no_barriers_pred() -> None:
    pred = NoBarriersPredicate()
    c = Circuit(1).H(0)
    assert pred.verify(c)
    c.add_barrier([0]).H(0)
    assert not pred.verify(c)


def test_decompose_routing_gates_to_cxs() -> None:
    circ = Circuit(4)
    circ.CX(1, 0)
    circ.SWAP(0, 1)
    circ.SWAP(1, 2)
    circ.CX(2, 3)

    cu = CompilationUnit(circ)

    arc = Architecture([[0, 1], [1, 2], [2, 3]])
    pss = DecomposeSwapsToCXs(arc)

    assert pss.apply(cu)
    circ1 = cu.circuit

    for cmd in circ1.get_commands():
        assert cmd.op.type == OpType.CX


def test_remove_barriers() -> None:
    circ = Circuit(4)
    circ.CX(0, 1)
    circ.CX(1, 2)
    circ.add_barrier([1, 2])
    circ.CX(2, 3)

    cu = CompilationUnit(circ)
    pss = RemoveBarriers()

    assert pss.apply(cu)
    circ1 = cu.circuit

    for cmd in circ1.get_commands():
        assert cmd.op.type == OpType.CX


def test_user_defined_swap_decomp() -> None:
    circ = Circuit(2)
    circ.SWAP(0, 1)

    cu = CompilationUnit(circ)

    repcirc = Circuit(2)
    repcirc.X(0)
    repcirc.CX(0, 1)
    repcirc.CX(1, 0)
    repcirc.CX(0, 1)
    repcirc.X(1)

    pss = DecomposeSwapsToCircuit(repcirc)

    assert pss.apply(cu)
    circ1 = cu.circuit

    assert circ1 == repcirc


def test_pauligraph_synth() -> None:
    circ = Circuit(4, 4)
    pg = PauliExpBox([Pauli.X, Pauli.Z, Pauli.Y, Pauli.I], 0.3)
    circ.add_pauliexpbox(pg, [0, 1, 2, 3])
    circ.measure_all()
    cu = CompilationUnit(circ)

    pss = PauliSimp(PauliSynthStrat.Sets, CXConfigType.Tree)
    assert pss.apply(cu)
    circ1 = cu.circuit
    assert circ1.depth_by_type(OpType.CX) == 4


def test_squash_chains() -> None:
    # XY
    c = Circuit(2)
    c.H(0).H(1)
    c.Rx(0.1, 0).Ry(0.2, 0).Rx(0.3, 0).Ry(0.4, 0).Ry(0.5, 0).Rx(0.6, 0)
    c.CX(0, 1).H(1)
    c.Ry(0.7, 1).Rx(0.8, 1).Rx(0.9, 1).Ry(1.1, 1).Rx(1.2, 1)
    u = c.get_unitary()
    EulerAngleReduction(OpType.Rx, OpType.Ry, strict=True).apply(c)
    u1 = c.get_unitary()
    assert np.allclose(u, u1)
    optypes = [cmd.op.type for cmd in c.get_commands()]
    assert optypes == [
        OpType.H,
        OpType.H,
        OpType.Ry,
        OpType.Rx,
        OpType.Ry,
        OpType.CX,
        OpType.H,
        OpType.Ry,
        OpType.Rx,
        OpType.Ry,
    ]
    # XZ
    c = Circuit(2)
    c.H(0).H(1)
    c.Rx(0.1, 0).Rz(0.2, 0).Rx(0.3, 0).Rz(0.4, 0).Rz(0.5, 0).Rx(0.6, 0)
    c.CX(0, 1).H(1)
    c.Rz(0.7, 1).Rx(0.8, 1).Rx(0.9, 1).Rz(1.1, 1).Rx(1.2, 1)
    u = c.get_unitary()
    EulerAngleReduction(OpType.Rx, OpType.Rz, strict=True).apply(c)
    u1 = c.get_unitary()
    assert np.allclose(u, u1)
    optypes = [cmd.op.type for cmd in c.get_commands()]
    assert optypes == [
        OpType.H,
        OpType.H,
        OpType.Rz,
        OpType.Rx,
        OpType.Rz,
        OpType.CX,
        OpType.H,
        OpType.Rz,
        OpType.Rx,
        OpType.Rz,
    ]
    # YZ
    c = Circuit(2)
    c.H(0).H(1)
    c.Ry(0.1, 0).Rz(0.2, 0).Ry(0.3, 0).Rz(0.4, 0).Rz(0.5, 0).Ry(0.6, 0)
    c.CX(0, 1).H(1)
    c.Rz(0.7, 1).Ry(0.8, 1).Ry(0.9, 1).Rz(1.1, 1).Ry(1.2, 1)
    u = c.get_unitary()
    EulerAngleReduction(OpType.Ry, OpType.Rz, strict=True).apply(c)
    u1 = c.get_unitary()
    assert np.allclose(u, u1)
    optypes = [cmd.op.type for cmd in c.get_commands()]
    assert optypes == [
        OpType.H,
        OpType.H,
        OpType.Rz,
        OpType.Ry,
        OpType.Rz,
        OpType.CX,
        OpType.H,
        OpType.Rz,
        OpType.Ry,
        OpType.Rz,
    ]


def test_apply_pass_with_callbacks() -> None:
    class CallbackHandler:
        def __init__(self) -> None:
            self.pass_names: List[str] = []

        def before_apply(self, cu: CompilationUnit, config: Dict[str, Any]) -> None:
            if "StandardPass" in config:
                self.pass_names.append(config["StandardPass"]["name"])
            else:
                self.pass_names.append(config["pass_class"])

        def after_apply(self, cu: CompilationUnit, config: Dict[str, Any]) -> None:
            return

    def compile(circ: Circuit, handler: CallbackHandler) -> bool:
        p = SequencePass([CommuteThroughMultis(), RemoveRedundancies()])
        return p.apply(circ, handler.before_apply, handler.after_apply)  # type: ignore

    circ = Circuit(5)
    circ.CX(0, 1)
    circ.CX(2, 4)
    circ.CX(0, 1)

    handler = CallbackHandler()
    compile(circ, handler)

    assert circ.n_gates_of_type(OpType.CX) == 1
    assert len(handler.pass_names) == 3
    assert handler.pass_names[0] == "SequencePass"
    assert handler.pass_names[1] == "CommuteThroughMultis"
    assert handler.pass_names[2] == "RemoveRedundancies"


def test_remove_discarded() -> None:
    c = Circuit(3, 2)
    c.H(0).H(1).H(2).CX(0, 1).Measure(0, 0).Measure(1, 1).H(0).H(1)
    c.qubit_discard(Qubit(0))
    c.qubit_discard(Qubit(2))
    assert not c.qubit_is_discarded(Qubit(1))
    assert c.qubit_is_discarded(Qubit(2))
    assert RemoveDiscarded().apply(c)
    assert c.n_gates_of_type(OpType.H) == 3
    assert c.n_gates_of_type(OpType.CX) == 1
    assert c.n_gates_of_type(OpType.Measure) == 2


def test_simplify_measured() -> None:
    c = Circuit(3, 3)
    c.H(0).H(1).Z(2)
    u = np.array(
        [
            [0, 0, 0, np.exp(0.1j)],
            [np.exp(0.2j), 0, 0, 0],
            [0, np.exp(0.3j), 0, 0],
            [0, 0, np.exp(0.4j), 0],
        ],
        dtype=complex,
    )
    ubox = Unitary2qBox(u)
    c.add_unitary2qbox(ubox, 0, 1)
    c.measure_all()
    c.qubit_discard(Qubit(0))
    c.qubit_discard(Qubit(1))
    assert SimplifyMeasured().apply(c)
    assert c.n_gates_of_type(OpType.H) == 2
    assert c.n_gates_of_type(OpType.Z) == 1
    assert c.n_gates_of_type(OpType.Unitary2qBox) == 0
    assert c.n_gates_of_type(OpType.Measure) == 3
    assert c.n_gates_of_type(OpType.ClassicalTransform) == 1


def test_simplify_initial_1() -> None:
    c = Circuit(4)
    c.H(0).X(1).CY(1, 2).CX(0, 1).CX(2, 3).H(1).H(2)
    assert not c.qubit_is_created(Qubit(0))
    c.qubit_create(Qubit(0))
    assert c.qubit_is_created(Qubit(0))
    c.qubit_create(Qubit(1))
    c.qubit_create(Qubit(2))
    assert SimplifyInitial().apply(c)
    assert c.n_gates_of_type(OpType.CY) == 0
    assert c.n_gates_of_type(OpType.CX) == 2


def test_simplify_initial_2() -> None:
    c = Circuit(1, 1).Y(0).measure_all()
    c.qubit_create_all()
    c1 = c.copy()
    assert SimplifyInitial(False).apply(c1)
    assert c1.n_gates_of_type(OpType.Y) == 0
    assert c1.n_gates_of_type(OpType.X) == 1
    assert c1.n_gates_of_type(OpType.Measure) == 1
    assert c1.n_gates_of_type(OpType.SetBits) == 0
    c2 = c.copy()
    xcirc = Circuit(1).Rx(1, 0)
    assert SimplifyInitial(xcirc=xcirc).apply(c2)
    assert c2.n_gates_of_type(OpType.Y) == 0
    assert c2.n_gates_of_type(OpType.Rx) == 1
    assert c2.n_gates_of_type(OpType.Measure) == 0
    assert c2.n_gates_of_type(OpType.SetBits) == 1


def test_simplify_initial_3() -> None:
    c = Circuit(2).X(0).CX(0, 1).CX(1, 0).X(1).CX(1, 0)
    c0 = c.copy()
    assert SimplifyInitial(create_all_qubits=True, remove_redundancies=False).apply(c0)
    c0_cmds = c0.get_commands()
    assert len(c0_cmds) > 0
    assert all(cmd.op.type == OpType.X for cmd in c0_cmds)
    c1 = c.copy()
    assert SimplifyInitial(create_all_qubits=True, remove_redundancies=True).apply(c1)
    c1_cmds = c1.get_commands()
    assert len(c1_cmds) == 0


def test_simplify_initial_symbolic() -> None:
    # Symbolic circuits should be left alone, no error
    c = Circuit(1)
    c.qubit_create(Qubit(0))
    c.Rx(Symbol("a"), 0)  # type: ignore
    c.measure_all()
    c1 = c.copy()
    SimplifyInitial(allow_classical=False).apply(c)
    assert c1 == c


def test_ZZPhaseToRz() -> None:
    c = (
        Circuit(2)
        .ZZPhase(0.6, 0, 1)
        .ZZPhase(1, 0, 1)
        .ZZPhase(-1, 0, 1)
        .ZZPhase(-0.4, 0, 1)
    )
    comp = (
        Circuit(2)
        .ZZPhase(0.6, 0, 1)
        .Rz(1, 0)
        .Rz(1, 1)
        .Rz(1, 0)
        .Rz(1, 1)
        .ZZPhase(-0.4, 0, 1)
    )
    ZZPhaseToRz().apply(c)
    assert comp == c


def test_pauli_squash() -> None:
    c = Circuit(3)
    c.add_pauliexpbox(PauliExpBox([Pauli.Z, Pauli.X, Pauli.Z], 0.8), [0, 1, 2])
    c.add_pauliexpbox(PauliExpBox([Pauli.Y, Pauli.X, Pauli.X], 0.2), [0, 1, 2])
    for strat in [
        PauliSynthStrat.Individual,
        PauliSynthStrat.Pairwise,
        PauliSynthStrat.Sets,
    ]:
        for cx_config in [CXConfigType.Snake, CXConfigType.Star, CXConfigType.Tree]:
            c1 = c.copy()
            assert PauliSquash().apply(c1)
            assert c1.n_gates_of_type(OpType.CX) <= 4


def test_three_qubit_squash() -> None:
    c = Circuit(3)
    for i in range(21):
        c.H(i % 3)
        c.CX(i % 3, (i + 1) % 3)
    c.measure_all()
    assert ThreeQubitSquash().apply(c)
    assert c.n_gates_of_type(OpType.CX) <= 18


def test_cnx_pairwise_decomp() -> None:
    c = Circuit(6).add_gate(OpType.CnX, [], [0, 1, 2, 3, 4, 5])
    c.add_gate(OpType.CnX, [], [1, 2, 3, 4, 5, 0])
    c.add_gate(OpType.CnX, [], [3, 1, 4, 5, 0, 2])
    CnXPairwiseDecomposition().apply(c)
    DecomposeMultiQubitsCX().apply(c)
    assert c.n_gates_of_type(OpType.CX) < 217


def test_rz_phasedX_squash() -> None:
    c = Circuit(2)
    c.Rz(0.3, 0)
    c.Rz(0.7, 1)
    c.ZZMax(0, 1)
    c.ZZMax(1, 0)
    c.add_gate(OpType.PhasedX, [0.2, 1.3], [0])
    c.add_gate(OpType.PhasedX, [0.5, 1.7], [1])
    c.ZZMax(0, 1)
    c.ZZMax(1, 0)

    assert SquashRzPhasedX().apply(c)
    assert c.n_gates_of_type(OpType.Rz) == 2
    cmds = c.get_commands()
    assert cmds[-1].op.type == OpType.Rz
    assert cmds[-2].op.type == OpType.Rz


def test_conditional_phase() -> None:
    # A conditional H cannot be expressed in terms of TK1 because of the global phase.
    c = Circuit(2, 2)
    c.H(0)
    c.Measure(0, 0)
    c.H(1, condition_bits=[0], condition_value=1)
    c.Measure(1, 1)
    target_gateset = {OpType.TK1, OpType.CX}
    rebase = auto_rebase_pass(target_gateset)
    rebase.apply(c)
    cond_cmds = [cmd for cmd in c.get_commands() if cmd.op.type == OpType.Conditional]
    assert len(cond_cmds) > 0
    assert any(cond_cmd.op.op.type not in target_gateset for cond_cmd in cond_cmds)


if __name__ == "__main__":
    test_predicate_generation()
    test_compilation_unit_generation()
    test_compilerpass_seq()
    test_rebase_pass_generation()
    test_routing_and_placement_pass()
    test_default_mapping_pass()
    test_SynthesiseTket_creation()
    test_directed_cx_pass()
    test_decompose_routing_gates_to_cxs()
    test_user_defined_swap_decomp()
    test_squash_chains()
    test_apply_pass_with_callbacks()
    test_remove_barriers()
    test_RebaseOQC_and_SynthesiseOQC()
    test_ZZPhaseToRz()
