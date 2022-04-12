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
from pytket.circuit.logic_exp import BitWiseOp


CustomPyQIRGate = NamedTuple(
    "CustomPyQIRGate",
    [
        ("optype", str,),
        (
            "functions",
            List[Union[type[types.DOUBLE], type[types.QUBIT], type[types.RESULT]]],
        ),
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
    BitWiseOp.AND: "and",
    BitWiseOp.OR: "or",
    BitWiseOp.XOR: "xor",
}

QUANTINUUM_GATES = GateSet(
    name="Quantinuum",
    template=Template("__quantinuum__${optype}__${opname}__body"),
    gateset={
        "h": CustomPyQIRGate(
            optype="qis",
            functions=[types.QUBIT]
        ),
        "x": CustomPyQIRGate(
            optype="qis",
            functions=[types.QUBIT]
        ),
        "y": CustomPyQIRGate(
            optype="qis",
            functions=[types.QUBIT]
        ),
        "z": CustomPyQIRGate(
            optype="qis",
            functions=[types.QUBIT]
        ),
        "rx": CustomPyQIRGate(
            optype="qis",
            functions=[types.DOUBLE, types.QUBIT]
        ),
        "ry": CustomPyQIRGate(
            optype="qis",
            functions=[types.DOUBLE, types.QUBIT]
        ),
        "rz": CustomPyQIRGate(
            optype="qis",
            functions=[types.DOUBLE, types.QUBIT]
        ),
        "phx": CustomPyQIRGate(
            optype="qis",
            functions=[types.DOUBLE, types.DOUBLE, types.QUBIT]
        ),
        "cnot": CustomPyQIRGate(
            optype="qis",
            functions=[types.QUBIT, types.QUBIT]
        ),
        "zzmax": CustomPyQIRGate(
            optype="qis",
            functions=[types.QUBIT, types.QUBIT]
        ),
        "zzph": CustomPyQIRGate(
            optype="qis",
            functions=[types.DOUBLE, types.QUBIT, types.QUBIT]
        ),
        "mz": CustomPyQIRGate(
            optype="qis",
            functions=[types.QUBIT, types.RESULT]
        ),
        "and": CustomPyQIRGate(
            optype="cis",
            functions=[types.RESULT, types.RESULT, types.RESULT]
        ),
        "or": CustomPyQIRGate(
            optype="cis",
            functions=[types.RESULT, types.RESULT, types.RESULT]
        ),
        "xor": CustomPyQIRGate(
            optype="cis",
            functions=[types.RESULT, types.RESULT, types.RESULT]
        )
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
    if optype == OpType.ExplicitPredicate:
        params = None
        if op.get_name() == "AND":
            optype = BitWiseOp.AND
        elif op.get_name() == "OR":
            optype = BitWiseOp.OR
        elif op.get_name() == "XOR":
            optype = BitWiseOp.XOR
    else:
        params = op.params
        if optype == OpType.TK1:
            # convert to U3
            optype = OpType.U3
            params = [op.params[1], op.params[0] - 0.5, op.params[2] + 0.5]
    return (optype, params)


def _to_qis_qubits(qubits: List[Qubit], mod: SimpleModule) -> List[type[types.QUBIT]]:
    return [mod.qubits[qubit.index[0]] for qubit in qubits]


def _to_qis_results(bits: List[Bit], mod: SimpleModule) -> Optional[type[types.RESULT]]:
    if bits:
        return mod.results[bits[0].index[0]]
    return None


def circuit_from_qir(input_file) -> None:
    pass


def circuit_to_qir_str(
    circ: Circuit, module: Union[ExtendedModule, SimpleModule]
) -> str:
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
    for command in circ:
        op = command.op
        optype, params = _get_optype_and_params(op)
        if isinstance(module, ExtendedModule):
            mod = module.module
            qubits = _to_qis_qubits(command.qubits, mod)
            results = _to_qis_results(command.bits, mod)
            try:
                gate = module.gateset.tk_to_gateset(optype)
            except KeyError:
                raise KeyError(
                    "Gate not defined in {:} gate set.".format(module.gateset.name)
                )
            get_gate = getattr(module, gate)
            if params:
                mod.builder.call(get_gate, [*params, *qubits])
            elif results:
                mod.builder.call(get_gate, [*qubits, results])
            else:
                mod.builder.call(get_gate, qubits)
        elif isinstance(module, SimpleModule):
            mod = module
            qis = BasicQisBuilder(module.builder)
            qubits = _to_qis_qubits(command.qubits, module)
            results = _to_qis_results(command.bits, module)
            try:
                pyqir_gate = _tk_to_pyqir(optype)
            except KeyError:
                raise KeyError("Gate not defined in PyQIR gate set.")
            get_gate = getattr(qis, pyqir_gate)
            if params:
                get_gate(*params, *qubits)
            elif results:
                get_gate(*qubits, results)
            else:
                get_gate(*qubits)
    return str(mod.ir())


def circuit_to_qir(
    circ: Circuit, output_file: str, gateset: Optional[GateSet] = None
) -> None:
    """A method to generate a qir file from a tket circuit."""
    root, ext = os.path.splitext(os.path.basename(output_file))
    if ext != ".ll":
        raise ValueError("The file extension should be '.ll'.")
    if gateset is not None:
        module = ExtendedModule(
            name=root,
            num_qubits=circ.n_qubits,
            num_results=len(circ.bits),
            gateset=gateset,
        )
    else:
        module = SimpleModule(
            name=root, num_qubits=circ.n_qubits, num_results=len(circ.bits)
        )
    circ_qir_str = circuit_to_qir_str(circ, module)
    with open(output_file, "w") as out:
        out.write(circ_qir_str)
