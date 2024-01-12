# Copyright 2019-2024 Cambridge Quantum Computing
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

from typing import Optional, Tuple
from pytket.circuit import Circuit
from pytket.passes import ContextSimp
from pytket.transform import separate_classical


def prepare_circuit(
    circ: Circuit, allow_classical: bool = True, xcirc: Optional[Circuit] = None
) -> Tuple[Circuit, Circuit]:
    """
    Prepare a circuit for processing by a backend device.

    This method first makes all inputs into Create operations (assuming an initial all-
    zero state) and all outputs into Discard operations (so that the circuit can no
    longer be usefully extended or appended to another circuit). It then attempts to
    apply various simplifications that take advantage of the known initial state and the
    fact that any unmeasured state is discarded. Finally, it separates the circuit into
    two circuits, the first of which is to be run on the backend (after any further
    compilation has been applied), and the second of which is a pure-classical circuit
    (on the same bits) which encodes classical post-processing of the measurement
    results. This post-processing is applied automatically when you pass the classical
    circuit as the `ppcirc` argument to `BackendResult.get_counts()` or
    `BackendResult.get_shots()`.

    The original circuit is not modified by this method.

    :param circ: input circuit
    :param allow_classical: allow insertion of mid-circuit classical operations?
    :param xcirc: 1-qubit circuit implementing an X gate in the transformed circuit (if
        omitted, an X gate is used)
    :return: (c0, ppcirc) where c0 is the simplified circuit and ppcirc should be passed
        to `BackendResult.get_counts()` or `BackendResult.get_shots()` when retrieving
        the final results.
    """
    c = circ.copy()
    c.qubit_create_all()
    c.qubit_discard_all()
    ContextSimp(allow_classical, xcirc).apply(c)
    return separate_classical(c)
