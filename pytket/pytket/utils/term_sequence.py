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
from typing import TYPE_CHECKING, cast

from pytket import Circuit
from pytket.circuit import CircBox, PauliExpBox
from pytket.partition import (
    GraphColourMethod,
    PauliPartitionStrat,
    term_sequence,
)

if TYPE_CHECKING:
    from .._tket.unit_id import UnitID

from .operators import QubitPauliOperator


def gen_term_sequence_circuit(
    operator: QubitPauliOperator,
    reference_state: Circuit,
    partition_strat: PauliPartitionStrat = PauliPartitionStrat.CommutingSets,
    colour_method: GraphColourMethod = GraphColourMethod.Lazy,
) -> Circuit:
    """
    Sequences the terms of a :py:class:`~.QubitPauliOperator` :math:`P` to generate
    a circuit approximating :math:`e^{-i \\frac{\\pi}{2} P}`. This method
    performs Trotterisation on :math:`P` with a single Trotter step.

    This method uses a given partitioning strategy and a graph colouring
      method for term sequencing.

    The resulting Circuit will contain a sequence of CircBoxes. Each CircBox
    corresponds to a set of Pauli strings. Each exponentiated Pauli string
    in the set is realised as a PauliExpBox.

    The ordering of terms prioritises reducing the two qubit gate count of the
    circuit when the PauliSimp or GuidedPauliSimp passes are applied rather
    than minimising the trotter error.

    :param operator: The operator terms to sequence
    :param reference_state: reference state to add sequenced terms to.
    :param partition_strat: a Partition strategy
    :param colour_method: a graph colouring method
    """
    qps_list = list(operator._dict.keys())  # noqa: SLF001
    qps_list_list = term_sequence(qps_list, partition_strat, colour_method)
    n_qbs = reference_state.n_qubits
    circ = reference_state.copy()
    qbs = circ.qubits
    for out_qps_list in qps_list_list:
        circ_to_box = Circuit(n_qbs)
        for qps in out_qps_list:
            coeff = operator[qps]
            qps_map = qps.map
            if qps_map:
                qubits = []
                paulis = []
                for qb, pauli in qps_map.items():
                    qubits.append(qb)
                    paulis.append(pauli)
                pbox = PauliExpBox(paulis, coeff)
                circ_to_box.add_pauliexpbox(pbox, qubits)
            else:
                circ_to_box.add_phase(-coeff / 2)
        cbox = CircBox(circ_to_box)
        unit_ids = cast("list[UnitID]", qbs)
        circ.add_circbox(cbox, unit_ids)
    return circ
