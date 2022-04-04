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

from collections import Counter
import numpy as np
from pytket.backends.backend import Backend
from hypothesis import strategies, given, settings, HealthCheck
from hypothesis.strategies._internal import SearchStrategy

from pytket.circuit import Qubit, Circuit, OpType  # type: ignore
from pytket.pauli import Pauli, QubitPauliString  # type: ignore
from pytket.partition import PauliPartitionStrat, GraphColourMethod  # type: ignore
from pytket.transform import Transform  # type: ignore
from pytket.utils.expectations import (
    expectation_from_shots,
    expectation_from_counts,
    get_operator_expectation_value,
)
from pytket.utils.measurements import append_pauli_measurement, _all_pauli_measurements
from pytket.utils.results import (
    counts_from_shot_table,
    get_n_qb_from_statevector,
    probs_from_counts,
    probs_from_state,
    permute_basis_indexing,
    compare_statevectors,
)
from pytket.utils.outcomearray import OutcomeArray
from pytket.utils import QubitPauliOperator, Graph
from pytket.utils.symbolic import (
    circuit_apply_symbolic_statevector,
    circuit_to_symbolic_unitary,
)
import pytest  # type: ignore
import types
from sympy import symbols  # type: ignore
from typing import Any, Callable, Tuple, Dict, List
from simulator import TketSimShotBackend, TketSimBackend  # type: ignore


def test_append_measurements() -> None:
    c = Circuit(4)
    qps = QubitPauliString(
        [Qubit(i) for i in range(4)], [Pauli.Y, Pauli.I, Pauli.X, Pauli.Z]
    )
    append_pauli_measurement(qps, c)
    coms = c.get_commands()
    assert len(coms) == 5
    assert str(coms[0]) == "Measure q[3] --> c[2];"
    assert str(coms[1]) == "Rx(0.5) q[0];"
    assert str(coms[2]) == "H q[2];"
    assert str(coms[3]) == "Measure q[0] --> c[0];"
    assert str(coms[4]) == "Measure q[2] --> c[1];"


def test_append_measurements_err0() -> None:
    c = Circuit(2)
    qps = QubitPauliString(Qubit(2), Pauli.X)
    with pytest.raises(RuntimeError) as ex:
        append_pauli_measurement(qps, c)
    assert "Circuit does not contain unit with id: q[2]" in str(ex.value)


def test_all_paulis() -> None:
    c = Circuit(2)
    qps1 = QubitPauliString()  # This will be ignored, as it imposes no measurements
    qps2 = QubitPauliString(Qubit(0), Pauli.Z)
    qps3 = QubitPauliString(Qubit(1), Pauli.Z)
    qps4 = QubitPauliString(Qubit(0), Pauli.Z)
    qps4[Qubit(1)] = Pauli.Z
    op = QubitPauliOperator({qps1: 1, qps2: 1, qps3: 2, qps4: -1.0j})
    circs = _all_pauli_measurements(op, c)
    assert isinstance(circs, types.GeneratorType)
    assert len(list(circs)) == 3


def test_shots_to_counts() -> None:
    shot_table = np.asarray([[0, 0], [0, 1], [0, 0]])
    counts = counts_from_shot_table(shot_table)
    assert len(counts) == 2
    assert counts[(0, 0)] == 2
    assert counts[(0, 1)] == 1


def test_counts_to_probs() -> None:
    counts: Dict[Tuple[int, ...], int] = {(0, 0): 4, (0, 1): 1, (1, 1): 3}
    probs = probs_from_counts(counts)
    assert len(probs) == 3
    assert probs[(0, 0)] == 0.5
    assert probs[(0, 1)] == 0.125
    assert probs[(1, 1)] == 0.375


def test_state_to_probs() -> None:
    state = np.asarray([0.5 - 0.5j, 0.5 + 0.5j, 1e-5, 0.999e-5])
    probs = probs_from_state(state)
    assert len(probs) == 2
    assert np.isclose(probs[(0, 0)], 0.5)
    assert np.isclose(probs[(0, 1)], 0.5)


def test_n_qb_from_statevector() -> None:
    state = np.asarray([0.5, 0.5, 0.5, 0.5])
    assert get_n_qb_from_statevector(state) == 2
    state = np.asarray([1.0, 0.0])
    assert get_n_qb_from_statevector(state) == 1
    state = np.asarray([1.0])
    assert get_n_qb_from_statevector(state) == 0
    state = np.zeros(128)
    assert get_n_qb_from_statevector(state) == 7


def test_n_qb_from_statevector_err() -> None:
    state = np.asarray([0.5, 0.5, 0.5])
    with pytest.raises(ValueError) as ex:
        get_n_qb_from_statevector(state)
    assert "is not a power of 2" in str(ex.value)


def test_permute_state() -> None:
    state = np.asarray([0, 0.8, 0.6, 0, 0, 0, 0, 0])
    permuted = permute_basis_indexing(state, (1, 2, 0))
    assert (permuted == np.asarray([0, 0, 0.8, 0, 0.6, 0, 0, 0])).all()


def test_permute_state_err1() -> None:
    state = np.asarray([0, 0.8, 0.6, 0, 0, 0, 0, 0])
    with pytest.raises(ValueError) as ex:
        permute_basis_indexing(state, (0, 1))
    assert "Invalid permutation: length does not match number of qubits" in str(
        ex.value
    )


def test_permute_state_err2() -> None:
    state = np.asarray([0, 0.8, 0.6, 0, 0, 0, 0, 0])
    with pytest.raises(ValueError) as ex:
        permute_basis_indexing(state, (0, 1, 3))
    assert "Permutation is not a valid complete permutation." in str(ex.value)


def test_permute_state_err3() -> None:
    state = np.asarray([0, 0.8, 0.6, 0, 0, 0, 0, 0])
    with pytest.raises(ValueError) as ex:
        permute_basis_indexing(state, (0, 1, 0))
    assert "Permutation is not a valid complete permutation." in str(ex.value)


def test_permute_basis_indexing() -> None:
    dimensions = 3
    bases = 1 << dimensions
    matrix = np.arange(bases**2).reshape((bases, bases))
    new_matrix = permute_basis_indexing(matrix, (1, 2, 0))
    assert np.array_equal(new_matrix, matrix[[0, 4, 1, 5, 2, 6, 3, 7], :])


def test_shot_expectation() -> None:
    shot_table = np.asarray([[0, 0, 1], [0, 1, 0], [1, 0, 0], [1, 1, 1]])
    assert expectation_from_shots(shot_table) == -1.0
    shot_table = np.asarray([[0, 0, 0], [0, 1, 1], [1, 0, 1], [1, 1, 0]])
    assert expectation_from_shots(shot_table) == 1.0
    shot_table = np.asarray([[0, 0, 0], [0, 0, 1], [1, 1, 0], [1, 1, 1]])
    assert expectation_from_shots(shot_table) == 0.0


def test_count_expectation() -> None:
    counts: Dict[Tuple[int, ...], int] = {
        (0, 0, 1): 4,
        (0, 1, 0): 7,
        (1, 0, 0): 1,
        (1, 1, 1): 8,
    }
    assert expectation_from_counts(counts) == -1.0
    counts = {(0, 0, 0): 4, (0, 1, 1): 7, (1, 0, 1): 1, (1, 1, 0): 8}
    assert expectation_from_counts(counts) == 1.0
    counts = {(0, 0, 0): 4, (0, 0, 1): 7, (1, 1, 0): 1, (1, 1, 1): 8}
    assert expectation_from_counts(counts) == -0.5


def test_outcomearray() -> None:
    in_listA = [1, 0, 1, 1] * 3
    in_listB = [1, 1, 0, 1] * 3

    list2int = lambda lis, order=True: int(
        "".join(map(str, lis[:: (-1) ** (not order)])), 2
    )
    assert OutcomeArray.from_readouts([in_listA]).to_readouts().tolist() == [in_listA]
    assert all(
        list(tup) == in_listA
        for tup in OutcomeArray.from_readouts([in_listA] * 4).to_readouts()
    )

    outcomeA = OutcomeArray.from_readouts([in_listA])
    assert outcomeA.n_outcomes == 1
    assert outcomeA.to_intlist() == [list2int(in_listA)]  # type: ignore

    outcomeB = OutcomeArray.from_readouts([in_listB])
    outcome2D = OutcomeArray(np.array([outcomeA[0], outcomeB[0]]), outcomeA.width)
    assert outcome2D.n_outcomes == 2
    assert outcome2D.to_intlist() == [list2int(in_listA), list2int(in_listB)]  # type: ignore
    assert outcome2D.to_intlist(False) == [
        list2int(in_listA, False),  # type: ignore
        list2int(in_listB, False),  # type: ignore
    ]

    for big in (True, False):
        intlist = [12, 5]
        readout_array = np.array([[1, 1, 0, 0], [0, 1, 0, 1]])
        if not big:
            readout_array = np.fliplr(readout_array)  # type: ignore

        outcome_from_ints = OutcomeArray.from_ints(intlist, 4, big_endian=big)
        assert np.array_equal(outcome_from_ints.to_readouts(), readout_array)
        assert outcome_from_ints.to_intlist(big) == intlist
        assert outcome_from_ints.to_intlist(not big) != intlist

    outcomeRepeats = OutcomeArray.from_readouts([in_listA, in_listA, in_listB])

    counts = outcomeRepeats.counts()
    assert len(counts) == 2
    assert counts[outcomeA] == 2
    assert counts[outcomeB] == 1

    counts1D = outcomeA.counts()
    assert len(counts1D) == 1
    assert counts1D[outcomeA] == 1

    # 0 width outcomearrays
    readouts: List[List[int]] = [[] for _ in range(10)]
    empty_array = OutcomeArray.from_readouts(readouts)
    assert np.array_equal(empty_array.to_readouts(), readouts)
    assert empty_array.counts() == Counter(
        {OutcomeArray(np.zeros((1, 0), dtype=np.uint8), 0): 10}
    )


def test_small_pauli_partition_expectation() -> None:
    c = Circuit(2)
    c.X(1)
    qps1 = QubitPauliString(Qubit(0), Pauli.Z)
    qps2 = QubitPauliString(Qubit(1), Pauli.Z)
    op = QubitPauliOperator({qps1: 0.5, qps2: 1.0})
    backend = TketSimShotBackend()
    n_shots = 10000
    strats = [
        None,
        PauliPartitionStrat.NonConflictingSets,
        PauliPartitionStrat.CommutingSets,
    ]
    for strat in strats:
        energy = complex(
            get_operator_expectation_value(c, op, backend, n_shots, strat, seed=4)  # type: ignore
        )
        assert np.isclose(energy, -0.5, atol=0.01)


def test_medium_pauli_partition_expectation() -> None:
    c = Circuit(4)
    c.H(1)
    c.H(3)
    c.Z(3)
    qps1 = QubitPauliString({Qubit(0): Pauli.Z, Qubit(1): Pauli.Z, Qubit(2): Pauli.Z})
    qps2 = QubitPauliString(Qubit(0), Pauli.Y)
    qps3 = QubitPauliString({Qubit(1): Pauli.X, Qubit(3): Pauli.X})

    op = QubitPauliOperator({qps1: 0.5, qps2: 0.8, qps3: -10.2})
    backends = [TketSimShotBackend(), TketSimBackend()]
    n_shots_list = [10000, None]
    strats = [
        None,
        PauliPartitionStrat.NonConflictingSets,
        PauliPartitionStrat.CommutingSets,
    ]
    for backend, n_shots in zip(backends, n_shots_list):
        for strat in strats:
            energy = get_operator_expectation_value(
                c, op, backend, n_shots, strat, GraphColourMethod.LargestFirst, seed=456
            )
            assert np.isclose(float(np.real(energy)), 10.2, atol=0.01)  # type: ignore


def test_large_pauli_partition_expectation() -> None:
    c = Circuit(5)
    c.CX(0, 2)
    c.H(4)
    c.V(2)
    qps1 = QubitPauliString({Qubit(0): Pauli.Z, Qubit(1): Pauli.Z})
    qps2 = QubitPauliString({Qubit(0): Pauli.X, Qubit(2): Pauli.X})
    qps3 = QubitPauliString({Qubit(0): Pauli.Y, Qubit(2): Pauli.Y})
    qps4 = QubitPauliString(
        {Qubit(1): Pauli.Z, Qubit(2): Pauli.Z, Qubit(3): Pauli.X, Qubit(4): Pauli.Z}
    )
    qps5 = QubitPauliString({Qubit(3): Pauli.Z, Qubit(4): Pauli.X})
    qps6 = QubitPauliString()
    op = QubitPauliOperator(
        {qps1: 0.3, qps2: -0.7j, qps3: 0.9, qps4: 0.83, qps5: 0.5, qps6: 0.5}
    )
    backends = [TketSimShotBackend(), TketSimBackend()]
    n_shots_list = [10000, None]
    strats = [
        None,
        PauliPartitionStrat.NonConflictingSets,
        PauliPartitionStrat.CommutingSets,
    ]
    for backend, n_shots in zip(backends, n_shots_list):
        energy = [
            get_operator_expectation_value(
                c,
                op,
                backend,
                n_shots,
                strat,
                GraphColourMethod.LargestFirst,
                seed=3,
            )
            for strat in strats
        ]
        assert np.isclose(energy, [1.3, 1.3, 1.3], atol=0.02).all()
        energy2 = [
            get_operator_expectation_value(
                c, op, backend, n_shots, strat, GraphColourMethod.Lazy, seed=3
            )
            for strat in strats
        ]
        assert np.isclose(
            energy2,
            [1.3, 1.3, 1.3],
            atol=0.02,
        ).all()


def test_inversion_pauli_partition_expectation() -> None:
    c = Circuit(4)

    c.H(0)
    c.Vdg(1)
    c.Vdg(2)
    c.Vdg(3)

    qb_list = [Qubit(i) for i in range(4)]
    qps1 = QubitPauliString(qb_list, [Pauli.X, Pauli.Y, Pauli.Y, Pauli.Y])
    qps2 = QubitPauliString(qb_list, [Pauli.Y, Pauli.X, Pauli.Y, Pauli.Y])
    qps3 = QubitPauliString(qb_list, [Pauli.Y, Pauli.Y, Pauli.X, Pauli.Y])
    qps4 = QubitPauliString(qb_list, [Pauli.Y, Pauli.Y, Pauli.Y, Pauli.X])
    qps5 = QubitPauliString(qb_list, [Pauli.X, Pauli.X, Pauli.X, Pauli.Y])
    qps6 = QubitPauliString(qb_list, [Pauli.X, Pauli.X, Pauli.Y, Pauli.X])
    qps7 = QubitPauliString(qb_list, [Pauli.X, Pauli.Y, Pauli.X, Pauli.X])
    qps8 = QubitPauliString(qb_list, [Pauli.Y, Pauli.X, Pauli.X, Pauli.X])
    op = QubitPauliOperator(
        {
            qps1: 0.1,
            qps2: 0.2,
            qps3: 0.3,
            qps4: 0.4,
            qps5: 0.5,
            qps6: 0.6,
            qps7: 0.7,
            qps8: 0.8,
        }
    )
    backend = TketSimShotBackend()
    n_shots = 10000
    strats = [
        None,
        PauliPartitionStrat.NonConflictingSets,
        PauliPartitionStrat.CommutingSets,
    ]
    energy = [
        get_operator_expectation_value(
            c, op, backend, n_shots, strat, GraphColourMethod.Lazy, seed=54
        )
        for strat in strats
    ]
    assert np.isclose(energy, [0.04248, 0.04248, 0.08612], atol=0.01).all()
    energy2 = [
        get_operator_expectation_value(
            c, op, backend, n_shots, strat, GraphColourMethod.LargestFirst, seed=54
        )
        for strat in strats
    ]
    assert np.isclose(energy2, [0.04248, 0.04248, 0.08612], atol=0.01).all()


def test_compare_statevectors() -> None:
    test_vec = np.array([1 + 2 * 1j, 3 + 4 * 1j, 5 + 6 * 1j, 7 + 8 * 1j])
    other_vec = test_vec + (2 - 1.2 * 1j)
    test_vec /= np.sqrt(np.vdot(test_vec, test_vec))  # type: ignore
    other_vec /= np.sqrt(np.vdot(other_vec, other_vec))  # type: ignore

    assert compare_statevectors(test_vec, test_vec)
    assert not compare_statevectors(test_vec, other_vec)
    assert compare_statevectors(other_vec, other_vec)
    phase = np.exp(1j * 0.453)
    assert compare_statevectors(test_vec, phase * test_vec)


def test_dag() -> None:
    c = Circuit(4, 4)
    c.X(0)
    c.H(1)
    c.Rz(0.5, 1)
    c.CX(2, 0)
    c.CRz(0.5, 0, 3)
    c.Measure(3, 3)
    c.Measure(1, 1)
    c.Z(0, condition_bits=[3, 1], condition_value=2)
    c.H(0)
    D = Graph(c)
    Gnx = D.as_nx()
    assert len(Gnx) == 2 * (4 + 4) + c.n_gates
    G = D.get_DAG()
    assert G.directed
    assert G.name == "Circuit"
    Gqc = D.get_qubit_graph()
    assert not Gqc.directed
    assert Gqc.name == "Qubit connectivity"


def test_dag_implicit_perm() -> None:
    # THET-701
    c = Circuit(3).CX(0, 1).CX(1, 0).CX(1, 2).CX(2, 1)
    Transform.OptimiseCliffords().apply(c)
    assert not all(x == y for x, y in c.implicit_qubit_permutation().items())
    G = Graph(c)
    dag = G.get_DAG()
    assert dag.name == "Circuit"


@strategies.composite
def unitary_circuits(draw: Callable[[SearchStrategy[Any]], Any]) -> Circuit:
    # generate example symbolic unitary circuits
    n_qb = draw(strategies.integers(min_value=1, max_value=3))
    # available qubits as integers
    qb_strat = strategies.integers(min_value=0, max_value=n_qb - 1)
    # some symbols to sample from
    syms = symbols("a b c d e")  # type: ignore
    c = Circuit(n_qb)

    optype_dict = {
        typ: (1, 0)
        for typ in (
            OpType.Z,
            OpType.X,
            OpType.Y,
            OpType.S,
            OpType.Sdg,
            OpType.T,
            OpType.Tdg,
            OpType.V,
            OpType.Vdg,
            OpType.SX,
            OpType.SXdg,
            OpType.H,
        )
    }
    optype_dict.update(
        {typ: (1, 1) for typ in (OpType.Rx, OpType.Rz, OpType.Ry, OpType.U1)}
    )
    optype_dict.update({typ: (1, 2) for typ in (OpType.U2, OpType.PhasedX)})
    optype_dict.update({typ: (1, 3) for typ in (OpType.U3,)})

    optype_dict.update(
        {
            typ: (2, 0)
            for typ in (
                OpType.CX,
                OpType.CY,
                OpType.CZ,
                OpType.CH,
                OpType.CV,
                OpType.CVdg,
                OpType.CSX,
                OpType.CSXdg,
                OpType.SWAP,
                OpType.ISWAPMax,
                OpType.Sycamore,
                OpType.ZZMax,
            )
        }
    )
    optype_dict.update(
        {
            typ: (2, 1)
            for typ in (
                OpType.CRz,
                OpType.CRx,
                OpType.CRy,
                OpType.CU1,
                OpType.ISWAP,
                OpType.XXPhase,
                OpType.YYPhase,
                OpType.ZZPhase,
                OpType.ESWAP,
            )
        }
    )
    optype_dict.update({typ: (2, 2) for typ in (OpType.PhasedISWAP, OpType.FSim)})
    optype_dict.update({typ: (2, 3) for typ in (OpType.CU3,)})

    optype_dict.update(
        {typ: (3, 0) for typ in (OpType.CCX, OpType.CSWAP, OpType.BRIDGE)}
    )

    optype_dict.update({OpType.XXPhase3: (3, 1)})

    optype_strat = strategies.sampled_from(list(optype_dict.keys()))
    for _ in range(5):
        typ = draw(optype_strat.filter(lambda x: optype_dict[x][0] <= n_qb))
        params = [
            draw(strategies.sampled_from(syms)) for _ in range(optype_dict[typ][1])
        ]
        qbs = [draw(qb_strat)]
        for _ in range(1, optype_dict[typ][0]):
            qbs.append(draw(qb_strat.filter(lambda x: x not in qbs)))

        c.add_gate(typ, params, qbs)
    return c


def unitary_from_states(circ: Circuit, back: Backend) -> np.ndarray:
    # use statevector simulation to calculate unitary from all input basis states
    nqb = circ.n_qubits
    matdim = 1 << nqb
    outar = np.zeros((matdim, matdim), dtype=np.complex128)

    for i in range(matdim):
        bitstr = f"{i:0{nqb}b}"
        basis_circ = Circuit(nqb)
        for qb, val in enumerate(bitstr):
            if val == "1":
                basis_circ.X(qb)
        basis_circ.append(circ)
        outar[:, i] = back.run_circuit(basis_circ).get_state()

    return outar


# this is a _slow_ test, so examples are kept low
# deadline has to be None because sympy runtime is very unpredictable
@given(circ=unitary_circuits())
@settings(
    deadline=None, max_examples=20, suppress_health_check=[HealthCheck.data_too_large]
)
def test_symbolic_conversion(circ: Circuit) -> None:

    sym_state = circuit_apply_symbolic_statevector(circ)

    sym_unitary = circuit_to_symbolic_unitary(circ)

    free_symbs = circ.free_symbols()
    # bind random values to symbolic variables to test numeric equality
    bind_vals = np.random.rand(len(free_symbs))

    substitutions = [(sym, val) for sym, val in zip(free_symbs, bind_vals)]
    circ.symbol_substitution(dict(substitutions))
    sym_unitary = sym_unitary.subs(substitutions)  # type: ignore
    sym_state = sym_state.subs(substitutions)  # type: ignore

    numeric_unitary = np.array(sym_unitary).astype(np.complex128)
    numeric_state = np.array(sym_state).astype(np.complex128)

    simulated_state_without_optimisation = circ.get_statevector()
    assert np.allclose(
        numeric_state.T, simulated_state_without_optimisation, atol=1e-10
    )

    simulated_unitary_without_optimisation = circ.get_unitary()
    assert np.allclose(
        numeric_unitary, simulated_unitary_without_optimisation, atol=1e-10
    )

    back = TketSimBackend()
    circ = back.get_compiled_circuit(circ, 1)
    result = back.run_circuit(circ)
    simulated_state = result.get_state()
    assert np.allclose(numeric_state.T, simulated_state, atol=1e-10)

    simulated_unitary = unitary_from_states(circ, back)
    assert np.allclose(numeric_unitary, simulated_unitary, atol=1e-10)


if __name__ == "__main__":
    test_append_measurements()
    test_append_measurements_err0()
    test_all_paulis()
    test_shots_to_counts()
    test_counts_to_probs()
    test_state_to_probs()
    test_permute_state()
    test_n_qb_from_statevector_err()
    test_permute_state_err1()
    test_permute_state_err2()
    test_permute_state_err3()
    test_permute_basis_indexing()
    test_shot_expectation()
    test_count_expectation()
    test_small_pauli_partition_expectation()
    test_medium_pauli_partition_expectation()
    test_large_pauli_partition_expectation()
    test_inversion_pauli_partition_expectation()
    test_dag()
    test_dag_implicit_perm()
