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

from typing import Dict, Tuple, List
import numpy as np
import pytest  # type: ignore
from pytket import Circuit
from pytket.pauli import Pauli, QubitPauliString  # type: ignore
from pytket.circuit import fresh_symbol, OpType, Qubit  # type: ignore
from pytket.utils import gen_term_sequence_circuit, QubitPauliOperator
from pytket.transform import Transform  # type: ignore
from pytket.partition import PauliPartitionStrat, GraphColourMethod  # type: ignore


def test_basic_sequence() -> None:
    c = Circuit(2)
    c.X(0).X(1)
    a = fresh_symbol("a")
    q0 = Qubit(0)
    q1 = Qubit(1)
    op = QubitPauliOperator()
    op[QubitPauliString()] = a
    op[QubitPauliString(q0, Pauli.Z)] = 0.1
    op[QubitPauliString(q1, Pauli.Z)] = 1
    qps = QubitPauliString({q0: Pauli.Z, q1: Pauli.Z})
    op[qps] = 0.3
    circ1 = gen_term_sequence_circuit(op, c)
    assert circ1.n_gates_of_type(OpType.CircBox) == 1
    circ2 = circ1.copy()
    Transform.DecomposeBoxes().apply(circ2)
    assert circ2.n_gates_of_type(OpType.X) == 2
    assert circ2.n_gates_of_type(OpType.CX) == 2
    assert circ2.n_gates_of_type(OpType.Rz) == 3
    Transform.UCCSynthesis().apply(circ1)
    assert circ1.n_gates_of_type(OpType.CX) == 2
    circ1.symbol_substitution({a: 3.2})
    circ2.symbol_substitution({a: 3.2})
    assert np.allclose(circ1.get_statevector(), circ2.get_statevector())


def test_nontrivial_sequence() -> None:
    c = Circuit(5)
    c.X(0).X(2).X(4)
    q = [Qubit(i) for i in range(5)]

    qps1 = QubitPauliString([q[0], q[1]], [Pauli.X, Pauli.X])
    qps2 = QubitPauliString([q[0], q[1], q[2]], [Pauli.Z, Pauli.Z, Pauli.Z])
    qps3 = QubitPauliString([q[0], q[1], q[3]], [Pauli.Y, Pauli.Y, Pauli.Y])
    qps4 = QubitPauliString([q[0], q[2], q[3]], [Pauli.X, Pauli.X, Pauli.X])
    qps5 = QubitPauliString([q[0], q[1], q[2]], [Pauli.Y, Pauli.Z, Pauli.X])
    qps6 = QubitPauliString([q[0], q[1], q[4]], [Pauli.Z, Pauli.Y, Pauli.X])
    qps7 = QubitPauliString(
        [q[0], q[1], q[3], q[4]], [Pauli.Z, Pauli.Y, Pauli.Y, Pauli.Y]
    )

    op = QubitPauliOperator(
        {
            qps1: 1.4,
            qps2: 0.3,
            qps3: 0.7,
            qps4: 0.4,
            qps5: 0.99,
            qps6: 0.2,
            qps7: 1.1,
        }
    )

    # Take the same initial circuit and ops, and test by applying
    # all possible Pauli partition strategies and graph colourings.
    #
    # KEY: a PauliPartitionStrat
    # VALUE: another dictionary: GraphColourMethod -> expected results tuple
    expected_results = {
        PauliPartitionStrat.CommutingSets: {
            GraphColourMethod.LargestFirst: (3, 28, 24, 23),
            GraphColourMethod.Lazy: (3, 28, 20, 19),
            GraphColourMethod.Exhaustive: (3, 28, 20, 19),
        },
        PauliPartitionStrat.NonConflictingSets: {
            GraphColourMethod.LargestFirst: (6, 28, 28, 28),
            GraphColourMethod.Lazy: (6, 28, 28, 26),
            GraphColourMethod.Exhaustive: (6, 28, 28, 26),
        },
    }

    # We'll build the results up, and check at the end that they match exactly
    calculated_results: Dict[
        PauliPartitionStrat, Dict[GraphColourMethod, Tuple[int, ...]]
    ] = {}

    for strategy, colour_method_dict in expected_results.items():
        calculated_results[strategy] = {}
        for colour_method in colour_method_dict:
            circ = gen_term_sequence_circuit(op, c, strategy, colour_method)
            counts_list: List[int] = [circ.n_gates_of_type(OpType.CircBox)]

            circ_other = circ.copy()
            Transform.DecomposeBoxes().apply(circ_other)
            counts_list.append(circ_other.n_gates_of_type(OpType.CX))

            Transform.UCCSynthesis().apply(circ)
            counts_list.append(circ.n_gates_of_type(OpType.CX))

            Transform.OptimiseCliffords().apply(circ)
            counts_list.append(circ.n_gates_of_type(OpType.CX))

            calculated_results[strategy][colour_method] = tuple(counts_list)
            assert np.allclose(circ.get_statevector(), circ_other.get_statevector())

    assert calculated_results == expected_results


if __name__ == "__main__":
    test_basic_sequence()
    test_nontrivial_sequence()
