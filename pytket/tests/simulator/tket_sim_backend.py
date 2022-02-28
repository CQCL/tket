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

"""Methods to allow tket circuits to be run with the tket_sim simulator
"""
### WARNING: TketSimBackend accepts Measure gates but does not apply them
### TketSimBackend will not work properly if there are extra bits that are unwritten to.

from .tket_sim_wrapper import TketSimWrapper

from typing import TYPE_CHECKING, List, Optional, Sequence, Union, cast
from uuid import uuid4

import numpy as np
from pytket.circuit import BasisOrder, Circuit, OpType  # type: ignore
from pytket.backends.backend import Backend, KwargTypes
from pytket.backends.backend_exceptions import CircuitNotRunError
from pytket.backends.backendresult import BackendResult
from pytket.backends.resulthandle import ResultHandle, _ResultIdTuple
from pytket.backends.status import CircuitStatus, StatusEnum
from pytket.passes import (  # type: ignore
    BasePass,
    SequencePass,
    SynthesiseTket,
    FullPeepholeOptimise,
    DecomposeBoxes,
    SimplifyInitial,
    auto_rebase_pass,
)
from pytket.predicates import (  # type: ignore
    GateSetPredicate,
    NoClassicalControlPredicate,
    NoFastFeedforwardPredicate,
    NoMidMeasurePredicate,
    Predicate,
)
from pytket.utils import prepare_circuit
from pytket.utils.outcomearray import OutcomeArray
from pytket.utils.results import probs_from_state


_GATE_SET = {
    OpType.SWAP,
    OpType.CX,
    OpType.CZ,
    OpType.Rz,
    OpType.Rx,
    OpType.S,
    OpType.T,
    OpType.S,
    OpType.X,
    OpType.H,
}


class TketSimBackend(Backend):
    """Backend for running simulations with tket_sim."""

    _supports_shots = False
    _supports_counts = False
    _supports_state = True
    _supports_expectation = False
    _persistent_handles = False

    @property
    def _result_id_type(self) -> _ResultIdTuple:
        return (str,)

    @property
    def required_predicates(self) -> List[Predicate]:
        # We don't include a GateSetPredicate; don't list allowed gates.
        # It's only for testing, so let it raise an exception naturally
        # if an unknown gate is passed in.
        return [
            NoClassicalControlPredicate(),
            NoFastFeedforwardPredicate(),
        ]

    def rebase_pass(self) -> BasePass:
        return auto_rebase_pass(_GATE_SET)

    def default_compilation_pass(self, optimisation_level: int = 1) -> BasePass:
        assert optimisation_level in range(3)
        if optimisation_level == 0:
            return SequencePass([DecomposeBoxes(), self.rebase_pass()])
        elif optimisation_level == 1:
            return SequencePass(
                [DecomposeBoxes(), SynthesiseTket(), self.rebase_pass()]
            )
        else:
            return SequencePass(
                [DecomposeBoxes(), FullPeepholeOptimise(), self.rebase_pass()]
            )

    def process_circuits(
        self,
        circuits: Sequence[Circuit],
        n_shots: Optional[Union[int, Sequence[int]]] = None,
        valid_check: bool = True,
        **kwargs: KwargTypes,
    ) -> List[ResultHandle]:
        circuits = list(circuits)
        if valid_check:
            self._check_all_circuits(circuits)

        handle_list = []
        for circuit in circuits:
            state = circuit.get_statevector()
            handle = ResultHandle(str(uuid4()))
            res = BackendResult(q_bits=sorted(circuit.qubits), state=state)
            self._cache[handle] = {"result": res}
            handle_list.append(handle)
        return handle_list

    def circuit_status(self, handle: ResultHandle) -> CircuitStatus:
        if handle in self._cache:
            return CircuitStatus(StatusEnum.COMPLETED)
        raise CircuitNotRunError(handle)


class TketSimShotBackend(TketSimBackend):
    """Backend for running simulations on the QuEST simulator.
    Naively samples from the statevector to produce shots.
    """

    _supports_shots = True
    _supports_counts = True
    _supports_state = False
    _supports_contextual_optimisation = True

    def default_compilation_pass(self, optimisation_level: int = 1) -> BasePass:
        assert optimisation_level in range(3)
        if optimisation_level == 0:
            return SequencePass([DecomposeBoxes(), self.rebase_pass()])
        elif optimisation_level == 1:
            return SequencePass(
                [
                    DecomposeBoxes(),
                    SynthesiseTket(),
                    self.rebase_pass(),
                    SimplifyInitial(allow_classical=False, create_all_qubits=True),
                ]
            )
        else:
            return SequencePass(
                [
                    DecomposeBoxes(),
                    FullPeepholeOptimise(),
                    self.rebase_pass(),
                    SimplifyInitial(allow_classical=False, create_all_qubits=True),
                ]
            )

    def process_circuits(
        self,
        circuits: Sequence[Circuit],
        n_shots: Optional[Union[int, Sequence[int]]] = None,
        valid_check: bool = True,
        **kwargs: KwargTypes,
    ) -> List[ResultHandle]:
        circuits = list(circuits)
        n_shots_list: List[int] = []
        if hasattr(n_shots, "__iter__"):
            if any(n is None or n < 1 for n in cast(Sequence[Optional[int]], n_shots)):
                raise ValueError(
                    "Shots are artificially generated, specify a positive value "
                    "for all list entries of n_shots."
                )
            n_shots_list = cast(List[int], n_shots)
            if len(n_shots_list) != len(circuits):
                raise ValueError("The length of n_shots and circuits must match")
        else:
            if n_shots is None:
                raise ValueError("Parameter n_shots is required for this backend")
            # convert n_shots to a list
            n_shots_list = [cast(int, n_shots)] * len(circuits)

        if valid_check:
            self._check_all_circuits(circuits)

        postprocess = kwargs.get("postprocess", False)

        handle_list = []
        for circuit, n_shots_circ in zip(circuits, n_shots_list):
            qubits = sorted(circuit.qubits, reverse=True)
            bits = circuit.bits
            if postprocess:
                c0, ppcirc = prepare_circuit(circuit, allow_classical=False)
            else:
                c0, ppcirc = circuit, None
            sim = TketSimWrapper(c0)
            state = sim.get_state(basis=BasisOrder.dlo)
            choices, probs = zip(*probs_from_state(state).items())
            np.random.seed(cast(int, kwargs.get("seed")))
            sample_indices = np.random.choice(len(choices), p=probs, size=n_shots_circ)
            q_to_b = circuit.qubit_to_bit_map
            readouts = []
            for i in sample_indices:
                choice = choices[i]
                readout = [0] * len(bits)
                for j, val in enumerate(choice):
                    b = q_to_b.get(qubits[j])
                    if b is not None:
                        readout[bits.index(b)] = val
                readouts.append(readout)
            shots = OutcomeArray.from_readouts(readouts)
            result = BackendResult(
                q_bits=qubits, c_bits=bits, shots=shots, ppcirc=ppcirc
            )
            handle = ResultHandle(str(uuid4()))
            self._cache[handle] = {"result": result}
            handle_list.append(handle)
        return handle_list
