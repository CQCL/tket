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
from typing import Callable, Dict, List, NamedTuple, Optional, Tuple, Union

from pyqir.generator import SimpleModule, BasicQisBuilder, types  # type: ignore
from pytket import Circuit, OpType, Bit, Qubit
from pytket.circuit import Op  # type: ignore


CustomPyQIRGate = NamedTuple(
    "CustomPyQIRGate",
    [
        ("functions", List[Union[types.DOUBLE, types.QUBIT, types.RESULT]]),
    ],
)

GateSet = NamedTuple(
    "GateSet",
    [
        ("name", str),
        ("template", Template),
        ("gateset", Dict[str, CustomPyQIRGate]),
        ("tk_to_gateset", Callable),
    ],
)


_TK_TO_QUANTINUUM = {
    OpType.H: "h",
    OpType.X: "x",
    OpType.Y: "y",
    OpType.Z: "z",
    OpType.CX: "cnot",
    OpType.ZZMax: "zzmax",
    OpType.Measure: "mz",
    OpType.Rx: "rx",
    OpType.Ry: "ry",
    OpType.Rz: "rz",
    OpType.PhasedX: "phx",
    OpType.ZZPhase: "zzph",
}

QUANTINUUM_GATES = GateSet(
    name="Quantinuum",
    template=Template("__quantinuum__qis__${name}__body"),
    gateset={
        "h": CustomPyQIRGate(functions=[types.QUBIT]),
        "x": CustomPyQIRGate(functions=[types.QUBIT]),
        "y": CustomPyQIRGate(functions=[types.QUBIT]),
        "z": CustomPyQIRGate(functions=[types.QUBIT]),
        "rx": CustomPyQIRGate(functions=[types.DOUBLE, types.QUBIT]),
        "ry": CustomPyQIRGate(functions=[types.DOUBLE, types.QUBIT]),
        "rz": CustomPyQIRGate(functions=[types.DOUBLE, types.QUBIT]),
        "phx": CustomPyQIRGate(functions=[types.DOUBLE, types.DOUBLE, types.QUBIT]),
        "cnot": CustomPyQIRGate(functions=[types.QUBIT, types.QUBIT]),
        "zzmax": CustomPyQIRGate(functions=[types.QUBIT, types.QUBIT]),
        "zzph": CustomPyQIRGate(functions=[types.DOUBLE, types.QUBIT, types.QUBIT]),
        "mz": CustomPyQIRGate(functions=[types.QUBIT, types.RESULT]),
    },
    tk_to_gateset=lambda optype: {**_TK_TO_QUANTINUUM}[optype],
)

_TK_TO_PYQIR = {
    OpType.H: "h",
    OpType.X: "x",
    OpType.Y: "y",
    OpType.Z: "z",
    OpType.S: "s",
    OpType.Sdg: "s_adj",
    OpType.T: "t",
    OpType.Tdg: "t_adj",
    OpType.Reset: "reset",
    OpType.CX: "cx",
    OpType.CZ: "cz",
    OpType.Measure: "m",
    OpType.Rx: "rx",
    OpType.Ry: "ry",
    OpType.Rz: "rz",
}


def _tk_to_pyqir(optype: OpType):
    return _TK_TO_PYQIR[optype]


class QIRUnsupportedError(Exception):
    pass


class ExtendedModule:
    """Module extensions to account for any input gate set."""

    def __init__(
        self, name: str, num_qubits: int, num_results: int, gateset: GateSet
    ) -> None:
        self.module = SimpleModule(name, num_qubits, num_results)
        self.gateset = gateset
        for k, v in gateset.gateset.items():
            self.__setattr__(
                str(k),
                self.module.add_external_function(
                    gateset.template.substitute(name=k),
                    types.Function(v.functions, types.VOID),
                ),
            )


def _get_optype_and_params(op: Op) -> Tuple[OpType, Optional[List[float]]]:
    optype = op.type
    params = op.params
    if optype == OpType.TK1:
        # convert to U3
        optype = OpType.U3
        params = [op.params[1], op.params[0] - 0.5, op.params[2] + 0.5]
    return (optype, params)


def _to_qis_qubits(qubits: List[Qubit], mod: SimpleModule) -> List:
    return (mod.qubits[qubit.index[0]] for qubit in qubits)


def _to_qis_results(bits: List[Bit], mod: SimpleModule) -> SimpleModule.results:
    return mod.results[bits[0].index[0]]


def circuit_from_qir(input_file) -> None:
    pass


def circuit_to_qir_str(circ: Circuit, root: str, gateset: GateSet) -> str:
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
    module = ExtendedModule(
        name=root,
        num_qubits=circ.n_qubits,
        num_results=len(circ.bits),
        gateset=gateset
    ).module
    qis = BasicQisBuilder(module.builder)

    for command in circ:
        op = command.op
        qubits = _to_qis_qubits(command.qubits, module)
        print("qubits {}".format(qubits))
        optype, params = _get_optype_and_params(op)
        if optype == OpType.Measure:
            results = _to_qis_results(command.bits, module)
            add_gate = getattr(qis, _tk_to_qir_noparams_1q[optype])
            add_gate(*qubits, results)

        if optype in _tk_to_qir_noparams_1q:
            add_gate = getattr(qis, _tk_to_qir_noparams_1q[optype])
            add_gate(*qubits)
        elif optype in _tk_to_qir_noparams_2q:
            add_gate = getattr(qis, _tk_to_qir_noparams_2q[optype])
            add_gate(*qubits)
        elif optype in _tk_to_qir_params_1q:
            assert params
            add_gate = getattr(qis, _tk_to_qir_params_1q[optype])
            add_gate(params[0], *qubits)
        else:
            raise QIRUnsupportedError(
                "Cannot print command of type: {}".format(op.get_name())
            )
    return str(module.ir())


def circuit_to_qir(circ: Circuit, output_file: str, gateset: GateSet) -> None:
    """A method to generate a qir file from a tket circuit."""
    root, ext = os.path.splitext(os.path.basename(output_file))
    if ext != ".ll":
        raise ValueError("The file extension should be '.ll'.")
    circ_qir_str = circuit_to_qir_str(circ, root, gateset)
    with open(output_file, "w") as out:
        out.write(circ_qir_str)
