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

from typing import List, Dict, Any
import numpy as np
from pytket.circuit import Circuit, BasisOrder, OpType, Qubit  # type: ignore
from pytket.utils.results import (
    permute_basis_indexing,
    permute_rows_cols_in_unitary,
    permute_qubits_in_statevector,
)


class TketSimWrapper:
    """tket_sim always uses ILO BasisOrder. This wrapper allows DLO order also."""

    def __init__(self, circuit: Circuit):
        self._circ = circuit
        self._n_qubits = circuit.n_qubits

        # Dictionary mapping qubits to integers 0,1,2,... (in sorted order)
        self._qmap = {qb: index for index, qb in enumerate(sorted(circuit.qubits))}

        # Implicit permutation: self._qmap_perm lists the input qubits indices in the
        # order of their corresponding output qubit indices.
        implicit_perm = circuit.implicit_qubit_permutation()
        qmap_perm = [0] * self._n_qubits
        for qb in circuit.qubits:
            qmap_perm[self._qmap[implicit_perm[qb]]] = self._qmap[qb]
        self._qmap_perm = tuple(qmap_perm)

    def get_state(self, basis: BasisOrder = BasisOrder.ilo) -> np.ndarray:
        """One can directly use circ.get_statevector() if the BasisOrder is ILO."""
        statevector = self._circ.get_statevector()
        if basis == BasisOrder.dlo:
            # tketsim always uses BasisOrder.ilo, so we must convert.
            rev_perm = tuple(range(self._n_qubits - 1, -1, -1))
            statevector = permute_basis_indexing(statevector, rev_perm)
        statevector = permute_qubits_in_statevector(statevector, self._qmap_perm)
        return statevector  # type: ignore
