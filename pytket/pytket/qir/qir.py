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
import platform
from string import Template
from typing import Any, Callable, Dict, List, NamedTuple, Optional, Tuple, Union

from pytket import Circuit, OpType, Bit, Qubit
from pytket.circuit import Op  # type: ignore
from pytket.circuit.logic_exp import BitWiseOp
if platform.machine() == "x86_64":
    from pyqir.parser import QirModule, QirCallInstr  # type: ignore
    from pyqir.parser._native import PyQirInstruction  # type: ignore
    from pyqir.generator import SimpleModule, BasicQisBuilder, types  # type: ignore


CustomPyQIRGate = NamedTuple(
    "CustomPyQIRGate",
    [
        (
            "opnat",
            str,
        ),
        (
            "function_signature",
            List[Union[type[types.DOUBLE], type[types.QUBIT], type[types.RESULT]]],
        ),
        (
            "return_type",
            Union[
                type[types.DOUBLE],
                type[types.QUBIT],
                type[types.RESULT],
                type[types.VOID],
            ],
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
    template=Template("__quantinuum__${opnat}__${opname}__body"),
    gateset={
        "h": CustomPyQIRGate(
            opnat="qis", function_signature=[types.QUBIT], return_type=types.VOID
        ),
        "x": CustomPyQIRGate(
            opnat="qis", function_signature=[types.QUBIT], return_type=types.VOID
        ),
        "y": CustomPyQIRGate(
            opnat="qis", function_signature=[types.QUBIT], return_type=types.VOID
        ),
        "z": CustomPyQIRGate(
            opnat="qis", function_signature=[types.QUBIT], return_type=types.VOID
        ),
        "rx": CustomPyQIRGate(
            opnat="qis",
            function_signature=[types.DOUBLE, types.QUBIT],
            return_type=types.VOID,
        ),
        "ry": CustomPyQIRGate(
            opnat="qis",
            function_signature=[types.DOUBLE, types.QUBIT],
            return_type=types.VOID,
        ),
        "rz": CustomPyQIRGate(
            opnat="qis",
            function_signature=[types.DOUBLE, types.QUBIT],
            return_type=types.VOID,
        ),
        "phx": CustomPyQIRGate(
            opnat="qis",
            function_signature=[types.DOUBLE, types.DOUBLE, types.QUBIT],
            return_type=types.VOID,
        ),
        "cnot": CustomPyQIRGate(
            opnat="qis",
            function_signature=[types.QUBIT, types.QUBIT],
            return_type=types.VOID,
        ),
        "zzmax": CustomPyQIRGate(
            opnat="qis",
            function_signature=[types.QUBIT, types.QUBIT],
            return_type=types.VOID,
        ),
        "zzph": CustomPyQIRGate(
            opnat="qis",
            function_signature=[types.DOUBLE, types.QUBIT, types.QUBIT],
            return_type=types.VOID,
        ),
        "mz": CustomPyQIRGate(
            opnat="qis",
            function_signature=[types.QUBIT, types.RESULT],
            return_type=types.VOID,
        ),
        "and": CustomPyQIRGate(
            opnat="cis",
            function_signature=[types.RESULT, types.RESULT],
            return_type=types.RESULT,
        ),
        "or": CustomPyQIRGate(
            opnat="cis",
            function_signature=[types.RESULT, types.RESULT],
            return_type=types.RESULT,
        ),
        "xor": CustomPyQIRGate(
            opnat="cis",
            function_signature=[types.RESULT, types.RESULT],
            return_type=types.RESULT,
        ),
    },
    tk_to_gateset=lambda optype: {**_TK_TO_QUANTINUUM}[optype],
)

_TK_TO_PYQIR = {
    OpType.H: "__h__body",
    OpType.X: "__x__body",
    OpType.Y: "__y__body",
    OpType.Z: "__z__body",
    OpType.S: "__s__body",
    OpType.Sdg: "__s__adj",
    OpType.T: "__t__body",
    OpType.Tdg: "__t__adj",
    OpType.Reset: "__reset__body",
    OpType.CX: "__cnot__body",
    OpType.CZ: "__cz__body",
    OpType.Measure: "__mz__body",
    OpType.Rx: "__rx__body",
    OpType.Ry: "__ry__body",
    OpType.Rz: "__rz__body",
}


_PYQIR_TO_TKET = dict(((item[1], item[0]) for item in _TK_TO_PYQIR.items()))


class QIRParser:
    """A parser class to return a pytket circuit from a QIR file."""

    def __init__(self, file_path: str) -> None:
        self.module = QirModule(file_path)

    def get_required_qubits(self) -> int:
        interop_funcs = self.module.get_funcs_by_attr("InteropFriendly")
        qubits = interop_funcs[0].get_attribute_value("requiredQubits")
        if qubits is not None:
            return int(qubits)
        return 0

    def get_required_results(self) -> int:
        interop_funcs = self.module.get_funcs_by_attr("InteropFriendly")
        results = interop_funcs[0].get_attribute_value("requiredResults")
        if results is not None:
            return int(results)
        return 0

    def get_optype(self, instr: PyQirInstruction) -> OpType:
        call_func_name = instr.call_func_name
        for k, v in _PYQIR_TO_TKET.items():
            if k in call_func_name:
                return v

    def get_params(self, instr: PyQirInstruction) -> Optional[List[float]]:
        params = instr.call_func_params[0].constant.float_double_value
        if params is not None:
            return [params]
        return None

    def get_operation(self, instr: QirCallInstr) -> Op:
        optype = self.get_optype(instr.instr)
        params = self.get_params(instr.instr)
        if params is not None:
            return Op.create(optype, params)
        return Op.create(optype)

    def get_qubit_indices(self, instr: PyQirInstruction) -> List[int]:
        func_params = instr.call_func_params
        params: List = []
        for param in func_params:
            if param.constant.is_qubit:
                params.append(param.constant.qubit_static_id)
            else:
                params.append(param.constant.result_static_id)
        return params

    def to_circuit(self) -> Circuit:
        qubits = self.get_required_qubits()
        bits = self.get_required_results()
        circuit = Circuit(qubits, bits)
        entry_block = self.module.functions[0].get_block_by_name("entry")
        if entry_block is None:
            raise NotImplementedError("The QIR file does not contain an entry block.")
        instrs = entry_block.instructions
        for instr in instrs:
            op = self.get_operation(instr)
            unitids = self.get_qubit_indices(instr.instr)
            circuit.add_gate(op, unitids)
        return circuit


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
                    gateset.template.substitute(opnat=v.opnat, opname=k),
                    types.Function(v.function_signature, v.return_type),
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


def _to_qis_bits(
    args: List[Bit], mod: SimpleModule
) -> Optional[List[type[types.RESULT]]]:
    if args:
        return [mod.results[bit.index[0]] for bit in args[:-1]]
    return None


def circuit_from_qir(input_file: Union[str, "os.PathLike[Any]"]) -> None:
    ext = os.path.splitext(input_file)[-1]
    if ext not in [".ll", ".bc"]:
        raise TypeError("Can only convert .bc or .ll files")
    if isinstance(input_file, os.PathLike):
        input_file = str(input_file)
    qir_parser = QIRParser(input_file)
    return qir_parser.to_circuit()


def circuit_to_qir_str(
    circ: Circuit, module: Union[ExtendedModule, SimpleModule]
) -> str:
    """A method to generate a QIR string from a pytket circuit."""
    if any(
        circ.n_gates_of_type(typ)
        for typ in (
            OpType.RangePredicate,
            OpType.MultiBit,
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
            bits: Optional[List[type[types.RESULT]]] = None
            if type(optype) == BitWiseOp:
                bits = _to_qis_bits(command.args, mod)
            try:
                gate = module.gateset.tk_to_gateset(optype)
            except KeyError:
                raise KeyError(
                    "Gate not defined in {:} gate set.".format(module.gateset.name)
                )
            get_gate = getattr(module, gate)
            if bits is not None:
                mod.builder.call(get_gate, bits)
            elif params:
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
