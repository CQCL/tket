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

import os
from string import Template
from typing import Dict, List, NamedTuple, Optional, Tuple

from pyqir.generator import SimpleModule, BasicQisBuilder, types  # type: ignore
from pytket import Circuit, OpType, Bit, Qubit
from pytket.circuit import Op  # type: ignore


PyQIRGate = NamedTuple(
    "PyQIRGate",
    [
        ("function", List[type(types)]), 
    ]
)

GateSet = NamedTuple(
    "GateSet",
    [
        ("template", Template),
        ("gateset", Dict[str, PyQIRGate])
    ]
)

QUANTINUUM_GATES = GateSet(
    template=Template('__quantinuum__qis__${name}__body'),
    gateset={
        "h": PyQIRGate(function=[types.QUBIT]),
        "x": PyQIRGate(function=[types.QUBIT]),
        "y": PyQIRGate(function=[types.QUBIT]),
        "z": PyQIRGate(function=[types.QUBIT]),
        "rx": PyQIRGate(function=[types.DOUBLE, types.QUBIT]),
        "ry": PyQIRGate(function=[types.DOUBLE, types.QUBIT]),
        "rz": PyQIRGate(function=[types.DOUBLE, types.QUBIT]),
        "phx": PyQIRGate(function=[types.DOUBLE, types.DOUBLE, types.QUBIT]),
        "cnot": PyQIRGate(function=[types.QUBIT, types.QUBIT]),
        "zzmax": PyQIRGate(function=[types.QUBIT, types.QUBIT]),
        "zzph": PyQIRGate(function=[types.DOUBLE, types.QUBIT, types.QUBIT]),
        "mz": PyQIRGate(function=[types.QUBIT, types.RESULT]),
    }
)


# Natively in pyqir
QUANTINUUM_NOPARAM_1Q_COMMANDS = {
    "h": OpType.H,
    "m": OpType.Measure,
    "reset": OpType.Reset,
    "s": OpType.S,
    "t": OpType.T,
    "x": OpType.X,
    "y": OpType.Y,
    "z": OpType.Z,
}


NOPARAM_2Q_COMMANDS = {
    "cx": OpType.CX,
    "cz": OpType.CZ,
}


PARAM_1Q_COMMANDS = {
    "rx": OpType.Rx,
    "ry": OpType.Ry,
    "rz": OpType.Rz,  # Also belongs to H machine gate set 
}


# To be included for the H machines gate set
NOPARAM_INCLUDED = {
    "zzmax": OpType.ZZMax
}

PARAM_INCLUDED = {
    "U1q": OpType.PhasedX,
    "rzz": OpType.ZZPhase,
}


_tk_to_qir_noparams_1q = dict(((item[1], item[0]) for item in QUANTINUUM_NOPARAM_1Q_COMMANDS.items()))
_tk_to_qir_noparams_2q = dict(((item[1], item[0]) for item in NOPARAM_2Q_COMMANDS.items()))
_tk_to_qir_params_1q = dict(((item[1], item[0]) for item in PARAM_1Q_COMMANDS.items()))


class QIRUnsupportedError(Exception):
    pass


class ExtendedModule:
    """Module extensions to account for H series gate set."""

    
    def __init__(self, name: str, num_qubits: int, num_results: int, gateset: GateSet) -> None:
        self.module = SimpleModule(name, num_qubits, num_results)
        for k, v in gateset.gateset.items():
            self.__setattr__(
                str(k),
                self.module.add_external_function(
                    gateset.template.substitute(name=k),
                    types.Function(v.function, types.VOID)
                )
            )


def _get_optype_and_params(op: Op) -> Tuple[OpType, Optional[List[float]]]:
    optype = op.type
    params = op.params if optype in _tk_to_qir_params else None
    if optype == OpType.TK1:
        # convert to U3
        optype = OpType.U3
        params = [op.params[1], op.params[0] - 0.5, op.params[2] + 0.5]
    return (optype, params)


def _to_qis_qubit(qubits: List[Qubit], mod: SimpleModule) -> SimpleModule.qubits:
    return mod.qubits[qubits[0].index[0]]


def _to_qis_results(bits: List[Bit], mod: SimpleModule) -> SimpleModule.results:
    return mod.results[bits[0].index[0]]


def circuit_from_qir(input_file) -> None:
    pass


def circuit_to_qir_str(circ: Circuit, root: str) -> str:
    """A method to generate a QIR string from a pytket circuit."""
    if any(
        circ.n_gates_of_type(typ)
        for typ in (
            OpType.RangePredicate,
            OpType.MultiBit,
            OpType.ExplicitPredicate,
            OpType.ExplicitModifier,
            OpType.SetBits,
        )
    ):
        raise QIRUnsupportedError("Complex classical gates not supported.")
    module = SimpleModule(root, num_qubits=circ.n_qubits, num_results=len(circ.bits))
    qis = BasicQisBuilder(module.builder)

    for command in circ:
        op = command.op
        qubits = _to_qis_qubit(command.qubits, module)
        optype, params = _get_optype_and_params(op)
        if optype == OpType.Measure:
            results = _to_qis_results(command.bits, module)
            add_gate = getattr(qis, _tk_to_qir_noparams[optype])
            add_gate(qubits, results)

        elif optype in _tk_to_qir_noparams:
            add_gate = getattr(qis, _tk_to_qir_noparams[optype])
            add_gate(qubits)
        elif optype in _tk_to_qir_params:
            assert params
            add_gate = getattr(qis, _tk_to_qir_params[optype])
            add_gate(params[0], qubits)
        else:
            raise QIRUnsupportedError(
                "Cannot print command of type: {}".format(op.get_name())
            )
    return str(module.ir())


def circuit_to_qir(circ: Circuit, output_file: str) -> None:
    """A method to generate a qir file from a tket circuit."""
    root, ext = os.path.splitext(os.path.basename(output_file))
    if ext != ".ll":
        raise ValueError("The file extension should be '.ll'.")
    circ_qir_str = circuit_to_qir_str(circ, root)
    with open(output_file, "w") as out:
        out.write(circ_qir_str)
