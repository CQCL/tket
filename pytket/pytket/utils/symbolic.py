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

"""Collection of methods to calculate symbolic statevectors and unitaries,
for symbolic circuits. This uses the sympy.physics.quantum module and produces
sympy objects. The implementations are slow and scale poorly, so this is
only suitable for very small (up to 5 qubit) circuits."""

from collections.abc import Callable
from typing import cast

import numpy as np
import sympy
from sympy import (
    BlockDiagMatrix,
    BlockMatrix,
    Expr,
    I,
    Identity,
    ImmutableMatrix,
    Matrix,
    Mul,
    diag,
    eye,
    zeros,
)
from sympy.physics.quantum import gate as symgate
from sympy.physics.quantum import represent
from sympy.physics.quantum.qapply import qapply
from sympy.physics.quantum.qubit import Qubit, matrix_to_qubit
from sympy.physics.quantum.tensorproduct import matrix_tensor_product

from pytket.circuit import Circuit, Op, OpType

# gates that have an existing definition in sympy
_FIXED_GATE_MAP: dict[OpType, type[symgate.Gate]] = {
    OpType.H: symgate.HadamardGate,
    OpType.S: symgate.PhaseGate,
    OpType.CX: symgate.CNotGate,
    OpType.SWAP: symgate.SwapGate,
    OpType.T: symgate.TGate,
    OpType.X: symgate.XGate,
    OpType.Y: symgate.YGate,
    OpType.Z: symgate.ZGate,
}

ParamsType = list[Expr | float]
# Make sure the return matrix is Immutable https://github.com/sympy/sympy/issues/18733
SymGateFunc = Callable[[ParamsType], ImmutableMatrix]
SymGateMap = dict[OpType, SymGateFunc]

# Begin matrix definitions for symbolic OpTypes
# matches internal TKET definitions
# see OpType documentation


def _symb_controlled(target: SymGateFunc) -> SymGateFunc:
    return lambda x: ImmutableMatrix(BlockDiagMatrix(Identity(2), target(x)))


def _symb_rz(params: ParamsType) -> ImmutableMatrix:
    return ImmutableMatrix(
        [
            [sympy.exp(-I * (sympy.pi / 2) * params[0]), 0],
            [0, sympy.exp(I * (sympy.pi / 2) * params[0])],
        ]
    )


def _symb_rx(params: ParamsType) -> ImmutableMatrix:
    costerm = sympy.cos((sympy.pi / 2) * params[0])
    sinterm = -I * sympy.sin((sympy.pi / 2) * params[0])
    return ImmutableMatrix(
        [
            [costerm, sinterm],
            [sinterm, costerm],
        ]
    )


def _symb_ry(params: ParamsType) -> ImmutableMatrix:
    costerm = sympy.cos((sympy.pi / 2) * params[0])
    sinterm = sympy.sin((sympy.pi / 2) * params[0])
    return ImmutableMatrix(
        [
            [costerm, -sinterm],
            [sinterm, costerm],
        ]
    )


def _symb_u3(params: ParamsType) -> ImmutableMatrix:
    theta, phi, lam = params
    costerm = sympy.cos((sympy.pi / 2) * theta)
    sinterm = sympy.sin((sympy.pi / 2) * theta)
    return ImmutableMatrix(
        [
            [costerm, -sinterm * sympy.exp(I * sympy.pi * lam)],
            [
                sinterm * sympy.exp(I * sympy.pi * phi),
                costerm * sympy.exp(I * sympy.pi * (phi + lam)),
            ],
        ]
    )


def _symb_u2(params: ParamsType) -> ImmutableMatrix:
    return _symb_u3([0.5, *params])


def _symb_u1(params: ParamsType) -> ImmutableMatrix:
    return _symb_u3([0.0, 0.0, *params])


def _symb_tk1(params: ParamsType) -> ImmutableMatrix:
    return _symb_rz([params[0]]) * _symb_rx([params[1]]) * _symb_rz([params[2]])


def _symb_tk2(params: ParamsType) -> ImmutableMatrix:
    return (
        _symb_xxphase([params[0]])
        * _symb_yyphase([params[1]])
        * _symb_zzphase([params[2]])
    )


def _symb_iswap(params: ParamsType) -> ImmutableMatrix:
    alpha = params[0]
    costerm = sympy.cos((sympy.pi / 2) * alpha)
    sinterm = sympy.sin((sympy.pi / 2) * alpha)
    return ImmutableMatrix(
        [
            [1, 0, 0, 0],
            [0, costerm, I * sinterm, 0],
            [0, I * sinterm, costerm, 0],
            [0, 0, 0, 1],
        ]
    )


def _symb_phasediswap(params: ParamsType) -> ImmutableMatrix:
    p, alpha = params
    costerm = sympy.cos((sympy.pi / 2) * alpha)
    sinterm = I * sympy.sin((sympy.pi / 2) * alpha)
    phase = sympy.exp(2 * I * sympy.pi * p)
    return ImmutableMatrix(
        [
            [1, 0, 0, 0],
            [0, costerm, sinterm * phase, 0],
            [0, sinterm / phase, costerm, 0],
            [0, 0, 0, 1],
        ]
    )


def _symb_xxphase(params: ParamsType) -> ImmutableMatrix:
    alpha = params[0]
    c = sympy.cos((sympy.pi / 2) * alpha)
    s = -I * sympy.sin((sympy.pi / 2) * alpha)
    return ImmutableMatrix(
        [
            [c, 0, 0, s],
            [0, c, s, 0],
            [0, s, c, 0],
            [s, 0, 0, c],
        ]
    )


def _symb_yyphase(params: ParamsType) -> ImmutableMatrix:
    alpha = params[0]
    c = sympy.cos((sympy.pi / 2) * alpha)
    s = I * sympy.sin((sympy.pi / 2) * alpha)
    return ImmutableMatrix(
        [
            [c, 0, 0, s],
            [0, c, -s, 0],
            [0, -s, c, 0],
            [s, 0, 0, c],
        ]
    )


def _symb_zzphase(params: ParamsType) -> ImmutableMatrix:
    alpha = params[0]
    t = sympy.exp(I * (sympy.pi / 2) * alpha)
    return ImmutableMatrix(diag(1 / t, t, t, 1 / t))


def _symb_xxphase3(params: ParamsType) -> ImmutableMatrix:
    xxphase2 = _symb_xxphase(params)
    res1 = matrix_tensor_product(xxphase2, eye(2))
    res2 = Matrix(
        BlockMatrix(
            [
                [xxphase2[:2, :2], zeros(2), xxphase2[:2, 2:], zeros(2)],
                [zeros(2), xxphase2[:2, :2], zeros(2), xxphase2[:2, 2:]],
                [xxphase2[2:, :2], zeros(2), xxphase2[2:, 2:], zeros(2)],
                [zeros(2), xxphase2[2:, :2], zeros(2), xxphase2[2:, 2:]],
            ]
        )
    )
    res3 = matrix_tensor_product(eye(2), xxphase2)
    return ImmutableMatrix(res1 * res2 * res3)


def _symb_phasedx(params: ParamsType) -> ImmutableMatrix:
    alpha, beta = params

    return _symb_rz([beta]) * _symb_rx([alpha]) * _symb_rz([-beta])


def _symb_eswap(params: ParamsType) -> ImmutableMatrix:
    alpha = params[0]
    c = sympy.cos((sympy.pi / 2) * alpha)
    s = -I * sympy.sin((sympy.pi / 2) * alpha)
    t = sympy.exp(-I * (sympy.pi / 2) * alpha)

    return ImmutableMatrix(
        [
            [t, 0, 0, 0],
            [0, c, s, 0],
            [0, s, c, 0],
            [0, 0, 0, t],
        ]
    )


def _symb_fsim(params: ParamsType) -> ImmutableMatrix:
    alpha, beta = params
    c = sympy.cos(sympy.pi * alpha)
    s = -I * sympy.sin(sympy.pi * alpha)
    t = sympy.exp(-I * sympy.pi * beta)

    return ImmutableMatrix(
        [
            [1, 0, 0, 0],
            [0, c, s, 0],
            [0, s, c, 0],
            [0, 0, 0, t],
        ]
    )


def _symb_gpi(params: ParamsType) -> ImmutableMatrix:
    t = sympy.exp(I * sympy.pi * params[0])

    return ImmutableMatrix(
        [
            [0, 1 / t],
            [t, 0],
        ]
    )


def _symb_gpi2(params: ParamsType) -> ImmutableMatrix:
    t = sympy.exp(I * sympy.pi * params[0])
    c = 1 / sympy.sqrt(2)

    return c * ImmutableMatrix(
        [
            [1, -I / t],
            [-I * t, 1],
        ]
    )


def _symb_aams(params: ParamsType) -> ImmutableMatrix:
    alpha, beta, gamma = params
    c = sympy.cos(sympy.pi / 2 * alpha)
    s = sympy.sin(sympy.pi / 2 * alpha)
    s1 = -I * sympy.exp(I * sympy.pi * (-beta - gamma)) * s
    s2 = -I * sympy.exp(I * sympy.pi * (-beta + gamma)) * s
    s3 = -I * sympy.exp(I * sympy.pi * (beta - gamma)) * s
    s4 = -I * sympy.exp(I * sympy.pi * (beta + gamma)) * s

    return ImmutableMatrix(
        [
            [c, 0, 0, s1],
            [0, c, s2, 0],
            [0, s3, c, 0],
            [s4, 0, 0, c],
        ]
    )


# end symbolic matrix definitions


class SymGateRegister:
    """Static class holding mapping from OpType to callable generating symbolic matrix.
    Allows users to add their own definitions, or override existing definitions."""

    _g_map: SymGateMap = {  # noqa: RUF012
        OpType.Rx: _symb_rx,
        OpType.Ry: _symb_ry,
        OpType.Rz: _symb_rz,
        OpType.TK1: _symb_tk1,
        OpType.TK2: _symb_tk2,
        OpType.U1: _symb_u1,
        OpType.U2: _symb_u2,
        OpType.U3: _symb_u3,
        OpType.CRx: _symb_controlled(_symb_rx),
        OpType.CRy: _symb_controlled(_symb_ry),
        OpType.CRz: _symb_controlled(_symb_rz),
        OpType.CU1: _symb_controlled(_symb_u1),
        OpType.CU3: _symb_controlled(_symb_u3),
        OpType.ISWAP: _symb_iswap,
        OpType.PhasedISWAP: _symb_phasediswap,
        OpType.XXPhase: _symb_xxphase,
        OpType.YYPhase: _symb_yyphase,
        OpType.ZZPhase: _symb_zzphase,
        OpType.XXPhase3: _symb_xxphase3,
        OpType.PhasedX: _symb_phasedx,
        OpType.ESWAP: _symb_eswap,
        OpType.FSim: _symb_fsim,
        OpType.GPI: _symb_gpi,
        OpType.GPI2: _symb_gpi2,
        OpType.AAMS: _symb_aams,
    }

    @classmethod
    def register_func(cls, typ: OpType, f: SymGateFunc, replace: bool = False) -> None:
        """Register a callable for an optype.

        :param typ: OpType to register
        :param f: Callable for generating symbolic matrix.
        :param replace: Whether to replace existing entry, defaults to False
        """
        if typ not in cls._g_map or replace:
            cls._g_map[typ] = f

    @classmethod
    def get_func(cls, typ: OpType) -> SymGateFunc:
        """Get registered callable."""
        return cls._g_map[typ]

    @classmethod
    def is_registered(cls, typ: OpType) -> bool:
        """Check if type has a callable registered."""
        return typ in cls._g_map


def _op_to_sympy_gate(op: Op, targets: list[int]) -> symgate.Gate:
    # convert Op to sympy gate
    if op.type in _FIXED_GATE_MAP:
        return _FIXED_GATE_MAP[op.type](*targets)
    if op.is_gate():
        # check if symbolic definition is needed
        float_params = all(isinstance(p, float) for p in op.params)
    else:
        raise ValueError(
            f"Circuit can only contain unitary gates, operation {op} not valid."
        )

    # pytket matrix basis indexing is in opposite order to sympy
    targets.reverse()
    if (not float_params) and SymGateRegister.is_registered(op.type):
        u_mat = SymGateRegister.get_func(op.type)(op.params)
    else:
        try:
            # use internal tket unitary definition
            u_mat = ImmutableMatrix(op.get_unitary())
        except RuntimeError as e:
            # to catch tket failure to get Op unitary
            # most likely due to symbolic parameters.
            raise ValueError(
                f"{op.type} is not supported for symbolic conversion."
                " Try registering your own symbolic matrix representation"
                " with SymGateRegister.func."
            ) from e
    return symgate.UGate(targets, u_mat)


def circuit_to_symbolic_gates(circ: Circuit) -> Mul:
    """Generate a multiplication expression of sympy gates from Circuit

    :param circ: Input circuit
    :raises ValueError: If circ does not match a unitary operation.
    :return: Symbolic gate multiplication expression.
    """
    outmat = symgate.IdentityGate(0)
    nqb = circ.n_qubits
    qubit_map = {qb: nqb - 1 - i for i, qb in enumerate(circ.qubits)}
    for com in circ:
        op = com.op
        if op.type == OpType.Barrier:
            continue
        args = com.args
        try:
            targs = [qubit_map[q] for q in args]  # type: ignore
        except KeyError as e:
            raise ValueError(
                f"Gates can only act on qubits. Operation {com} not valid."
            ) from e
        gate = _op_to_sympy_gate(op, targs)

        outmat = gate * outmat

    for i in range(len(qubit_map)):
        outmat = symgate.IdentityGate(i) * outmat

    return outmat * sympy.exp(circ.phase * sympy.pi * I)


def circuit_to_symbolic_unitary(circ: Circuit) -> ImmutableMatrix:
    """Generate a symbolic unitary from Circuit.

    Unitary matches pytket default ILO BasisOrder.

    :param circ: Input circuit
    :return: Symbolic unitary.
    """
    gates = circuit_to_symbolic_gates(circ)
    nqb = circ.n_qubits
    try:
        return cast("ImmutableMatrix", represent(gates, nqubits=circ.n_qubits))
    except NotImplementedError:
        # sympy can't represent n>1 qubit unitaries very well
        # so if it fails we will just calculate columns using the statevectors
        # for all possible input basis states
        matrix_dim = 1 << nqb
        input_states = (Qubit(f"{i:0{nqb}b}") for i in range(matrix_dim))
        outmat = Matrix([])
        for col, input_state in enumerate(input_states):
            outmat = outmat.col_insert(col, represent(qapply(gates * input_state)))

        return ImmutableMatrix(outmat)


def circuit_apply_symbolic_qubit(circ: Circuit, input_qb: Expr) -> Qubit:
    """Apply circuit to an input state to calculate output symbolic state.

    :param circ: Input Circuit.
    :param input_qb: Sympy Qubit expression corresponding to a state.
    :return: Output state after circuit acts on input_qb.
    """
    gates = circuit_to_symbolic_gates(circ)

    return cast("Qubit", qapply(gates * input_qb))


def circuit_apply_symbolic_statevector(
    circ: Circuit, input_state: np.ndarray | ImmutableMatrix | None = None
) -> ImmutableMatrix:
    """Apply circuit to an optional input statevector
    to calculate output symbolic statevector.
    If no input statevector given, the all zero state is assumed.
    Statevector follows pytket default ILO BasisOrder.

    :param circ: Input Circuit.
    :param input_state: Input statevector as a column vector, defaults to None.
    :return: Symbolic state after circ acts on input_state.
    """
    if input_state:
        input_qb = matrix_to_qubit(input_state)
    else:
        input_qb = Qubit("0" * circ.n_qubits)
    return cast(
        "ImmutableMatrix",
        represent(circuit_apply_symbolic_qubit(circ, cast("Qubit", input_qb))),
    )
