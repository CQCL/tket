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

import numpy as np
import pytest  # type: ignore

from pytket.circuit import (  # type: ignore
    ProjectorAssertionBox,
    StabiliserAssertionBox,
    Circuit,
    Qubit,
)

from pytket.passes import (  # type: ignore
    DecomposeBoxes,
)

from pytket.pauli import PauliStabiliser, Pauli  # type: ignore
from simulator import TketSimShotBackend  # type: ignore


def test_assertion_init() -> None:
    circ = Circuit(5)
    P = np.asarray(
        [
            [0.5, 0, 0, 0.5],
            [0, 0, 0, 0],
            [0, 0, 0, 0],
            [0.5, 0, 0, 0.5],
        ]
    )
    p_box = ProjectorAssertionBox(P)
    circ.add_assertion(p_box, [0, 1], name="|bell>")
    circ.add_assertion(p_box, [Qubit(0), Qubit(1)], name="|bell> 2")

    # Projector that requires a ancilla

    P = np.asarray(
        [
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 0],
        ]
    )
    p_box = ProjectorAssertionBox(P)
    with pytest.raises(RuntimeError) as errorinfo:
        circ.add_assertion(p_box, [0, 1])
    assert "ancilla" in str(errorinfo.value)
    circ.add_assertion(p_box, [0, 1], ancilla=2)

    stabilisers = StabiliserAssertionBox(
        [
            PauliStabiliser([Pauli.X, Pauli.X], 1),
            PauliStabiliser([Pauli.Y, Pauli.Y], -1),
        ]
    )
    stabilisers_2 = StabiliserAssertionBox(["XX", "-YY"])

    assert stabilisers.get_stabilisers() == stabilisers_2.get_stabilisers()
    circ.add_assertion(stabilisers, [0, 1], 2, name="|stab1>")
    circ.add_assertion(stabilisers_2, [0, 1], 2, name="|stab2>")
    circ.add_assertion(stabilisers_2, [Qubit(0), Qubit(1)], Qubit(2), name="|stab3>")

    # Test the circuit can be decomposed

    assert DecomposeBoxes().apply(circ)


def test_assertion() -> None:
    # P =|00><00| tensor I + |111><111|
    P = np.zeros((8, 8))
    P[0, 0] = 1
    P[1, 1] = 1
    P[7, 7] = 1

    p_box = ProjectorAssertionBox(P)

    circ = Circuit(6)
    circ.X(0)
    circ.H(2)
    circ.X(1)
    circ.X(3)
    circ.add_assertion(p_box, [4, 5, 1], name="|001>")
    circ.add_assertion(p_box, [0, 1, 3], name="|111>")
    circ.add_assertion(p_box, [1, 3, 4], name="|110>")

    b = TketSimShotBackend()
    circ_compiled = b.get_compiled_circuit(circ)
    h = b.process_circuit(circ_compiled, n_shots=5)
    r = b.get_result(h)

    # TketSimShotBackend produces incorrect results with OpType::Measure
    # so we can't check if the assertions are successful
    assert "|001>" in r.get_debug_info()
    assert "|111>" in r.get_debug_info()
    assert "|110>" in r.get_debug_info()

    # TketSimShotBackend doesn't support RESET hence we can't
    # process circuits with stabiliser assertions and some
    # of the projector based assertions


if __name__ == "__main__":
    test_assertion_init()
    test_assertion()
