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

from typing import Union
from sympy import Expr
from pytket.circuit import Circuit, OpType

Param = Union[float, "Expr"]


def approx_0_mod_2(x: Param, eps: float = 1e-10) -> bool:
    """Check if parameter is approximately 0 mod 2 up to eps precision.

    :param param: Parameter, float or sympy expression.
    :type param: Param
    :param eps: Tolerance, defaults to 1e-10
    :type eps: float, optional
    :return: Approximately 0 boolean.
    :rtype: bool
    """
    if isinstance(x, Expr) and not x.is_constant():  # type: ignore
        return False
    x = float(x)
    x %= 2
    return min(x, 2 - x) < eps


def int_half(angle: float) -> int:
    """Assume angle is approximately an even integer, and return the half

    :param angle: Float angle
    :type angle: float
    :return: Integer half of angle
    :rtype: int
    """
    #
    two_x = round(angle)
    assert not two_x % 2
    return two_x // 2


def _TK1_to_RxRy(a: Param, b: Param, c: Param) -> Circuit:
    return Circuit(1).Rx(-0.5, 0).Ry(c, 0).Rx(b, 0).Ry(a, 0).Rx(0.5, 0)


def _TK1_to_X_SX_Rz(a: Param, b: Param, c: Param) -> Circuit:
    circ = Circuit(1)
    correction_phase = 0.0

    # all phase identities use, for integer k,
    # Rx(2k) = Rz(2k) = (-1)^{k}I

    # _approx_0_mod_2 checks if parameters are constant
    # so they can be assumed to be constant
    if approx_0_mod_2(b):
        circ.Rz(a + c, 0)
        # b = 2k, if k is odd, then Rx(b) = -I
        correction_phase += int_half(float(b))

    elif approx_0_mod_2(b + 1):
        # Use Rx(2k-1) = i(-1)^{k}X
        correction_phase += -0.5 + int_half(float(b) - 1)
        if approx_0_mod_2(a - c):
            circ.X(0)
            # a - c = 2m
            # overall operation is (-1)^{m}Rx(2k -1)
            correction_phase += int_half(float(a - c))

        else:
            circ.Rz(c, 0).X(0).Rz(a, 0)

    elif approx_0_mod_2(b - 0.5) and approx_0_mod_2(a) and approx_0_mod_2(c):
        # a = 2k, b = 2m+0.5, c = 2n
        # Rz(2k)Rx(2m + 0.5)Rz(2n) = (-1)^{k+m+n}e^{-i \pi /4} SX
        circ.SX(0)
        correction_phase += (
            int_half(float(b) - 0.5) + int_half(float(a)) + int_half(float(c)) - 0.25
        )
    elif approx_0_mod_2(b - 0.5):
        # Use SX.Rz(2m-0.5).SX = (-1)^{m}e^{i \pi /4} Rz(-0.5).SX.Rz(-0.5)
        circ.Rz(c, 0).SX(0).Rz(a, 0)
        correction_phase += int_half(float(b) - 0.5) - 0.25
    elif approx_0_mod_2(b + 0.5) and approx_0_mod_2(a) and approx_0_mod_2(c):
        # a = 2k, b = 2m-0.5, c = 2n
        # Rz(2k)Rx(2m - 0.5)Rz(2n) = (-1)^{k+m+n}e^{i \pi /4} X.SX
        circ.X(0).SX(0)
        correction_phase += (
            int_half(float(b) + 0.5) + int_half(float(a)) + int_half(float(c)) + 0.25
        )
    elif approx_0_mod_2(b + 0.5):
        # Use SX.Rz(2m+0.5).SX = (-1)^{m}e^{i \pi /4} Rz(0.5).SX.Rz(0.5)
        circ.Rz(c + 1, 0).SX(0).Rz(a + 1, 0)
        correction_phase += int_half(float(b) - 1.5) - 0.25
    elif approx_0_mod_2(a - 0.5) and approx_0_mod_2(c - 0.5):
        # Rz(2k + 0.5)Rx(b)Rz(2m + 0.5) = -i(-1)^{k+m}SX.Rz(1-b).SX
        circ.SX(0).Rz(1 - b, 0).SX(0)
        correction_phase += int_half(float(a) - 0.5) + int_half(float(c) - 0.5) - 0.5
    else:
        circ.Rz(c + 0.5, 0).SX(0).Rz(b - 1, 0).SX(0).Rz(a + 0.5, 0)
        correction_phase += -0.5

    circ.add_phase(correction_phase)
    return circ


def _TK1_to_U(a: Param, b: Param, c: Param) -> Circuit:
    circ = Circuit(1)
    circ.add_gate(OpType.U3, [b, a - 0.5, c + 0.5], [0])
    circ.add_phase(-0.5 * (a + c))
    return circ
