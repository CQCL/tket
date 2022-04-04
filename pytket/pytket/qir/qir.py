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
from typing import List, Optional, Tuple

from pyqir.generator import SimpleModule, BasicQisBuilder  # type: ignore
from pytket import Circuit, OpType, Bit, Qubit
from pytket.circuit import Op # type: ignore


NOPARAM_COMMANDS = {
    "cx": OpType.CX,
    "cz": OpType.CZ,
    "h": OpType.H,
    "m": OpType.Measure,
    "reset": OpType.Reset,
    "s": OpType.S,
    "s_adj": OpType.Sdg,
    "t": OpType.T,
    "t_adj": OpType.Tdg,
    "x": OpType.X,
    "y": OpType.Y,
    "z": OpType.Z,
}


PARAM_COMMANDS = {
    "rx": OpType.Rx,
    "ry": OpType.Ry,
    "rz": OpType.Rz,
}


_tk_to_qir_noparams = dict(((item[1], item[0]) for item in NOPARAM_COMMANDS.items()))
_tk_to_qir_params = dict(((item[1], item[0]) for item in PARAM_COMMANDS.items()))


class QIRUnsupportedError(Exception):
    pass


def _get_optype_and_params(op: Op) -> Tuple[OpType, Optional[List[float]]]:
    optype = op.type
    params = op.params if optype in _tk_to_qir_params else None
    if optype == OpType.TK1:
        # convert to U3
        optype = OpType.U3
        params = [op.params[1], op.params[0] - 0.5, op.params[2] + 0.5]
    return (optype, params)


def _to_qis_qubit(qubits: List[Qubit], mod: SimpleModule):
    return mod.qubits[qubits[0].index[0]]


def _to_qis_results(bits: List[Bit], mod: SimpleModule):
    return mod.results[bits[0].index[0]]


def circuit_from_qir(input_file):
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
    return module.ir()


def circuit_to_qir(circ: Circuit, output_file: str) -> None:
    """A method to generate a qir file from a tket circuit."""
    root, ext = os.path.splitext(os.path.basename(output_file))
    if ext != ".ll":
        raise ValueError("The file extension should be '.ll'.")
    circ_qir_str = circuit_to_qir_str(circ, root)
    with open(output_file, "w") as out:
        out.write(circ_qir_str)
