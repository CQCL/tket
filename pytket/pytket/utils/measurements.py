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

from collections.abc import Iterable

from pytket.circuit import Bit, Circuit
from pytket.pauli import Pauli, QubitPauliString

from .operators import QubitPauliOperator


def append_pauli_measurement(pauli_string: QubitPauliString, circ: Circuit) -> None:
    """Appends measurement instructions to a given circuit, measuring each qubit in a
    given basis.

    :param pauli_string: The pauli string to measure
    :param circ: Circuit to add measurement to.
    """
    measured_qbs = []
    for qb, p in pauli_string.map.items():
        if p == Pauli.I:
            continue
        measured_qbs.append(qb)
        if p == Pauli.X:
            circ.H(qb)
        elif p == Pauli.Y:
            circ.Rx(0.5, qb)
    for b_idx, qb in enumerate(measured_qbs):
        unit = Bit(b_idx)
        circ.add_bit(unit, False)
        circ.Measure(qb, unit)


def _all_pauli_measurements(
    operator: QubitPauliOperator, circ: Circuit
) -> Iterable[Circuit]:
    """For each term in the operator, yields a copy of the given circuit with the
    appropriate measurements on each qubit. The trivial term is omitted.

    :param operator: The operator
    :param circ: The circuit generating the desired state
    :return: List of circuits in order of term from the operator
    """
    for pauli_string in operator._dict.keys():  # noqa: SLF001, SIM118
        if not pauli_string.map:
            continue
        copy = circ.copy()
        append_pauli_measurement(pauli_string, copy)
        yield copy
