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

from typing import Set, Callable, Dict, FrozenSet
from pytket.circuit import Circuit, OpType
import pytket._tket.circuit_library as _library
from pytket.passes import RebaseCustom, SquashCustom, SquashRzPhasedX

from ._decompositions import Param, _TK1_to_X_SX_Rz, _TK1_to_RxRy, _TK1_to_U
from .._tket.passes import BasePass


class NoAutoRebase(Exception):
    """Automatic rebase could not be found."""


_CX_CIRCS: Dict[OpType, Callable[[], "Circuit"]] = {
    OpType.CX: _library.CX,
    OpType.ZZMax: _library.CX_using_ZZMax,
    OpType.XXPhase: _library.CX_using_XXPhase_0,
    OpType.ECR: _library.CX_using_ECR,
    OpType.CZ: _library.H_CZ_H,
}


def _TK2_using_TK2(a: Param, b: Param, c: Param) -> Circuit:
    return Circuit(2).TK2(a, b, c, 0, 1)


_TK2_CIRCS: Dict[OpType, Callable[[Param, Param, Param], "Circuit"]] = {
    OpType.TK2: _TK2_using_TK2,
    OpType.ZZPhase: _library.TK2_using_ZZPhase,
    OpType.CX: _library.TK2_using_CX,
    OpType.ZZMax: _library.TK2_using_ZZMax,
}

_TK2_CIRCS_WIRE_SWAP: Dict[OpType, Callable[[Param, Param, Param], "Circuit"]] = {
    OpType.TK2: _library.TK2_using_TK2_or_swap,
    OpType.ZZPhase: _library.TK2_using_ZZPhase_and_swap,
    OpType.CX: _library.TK2_using_CX_and_swap,
    OpType.ZZMax: _library.TK2_using_ZZMax_and_swap,
}


def get_cx_decomposition(gateset: Set[OpType]) -> Circuit:
    """Return a Circuit expressing a CX in terms of a two qubit gate in the
    gateset if one is available, raise an error otherwise.

    :param gateset: Target gate set.
    :type gateset: Set[OpType]
    :raises NoAutoRebase: No suitable CX decomposition found.
    :return: Decomposition circuit.
    :rtype: Circuit
    """
    if any((matching := k) in gateset for k in _CX_CIRCS):
        return _CX_CIRCS[matching]()
    raise NoAutoRebase("No known decomposition from CX to available gateset.")


def get_tk2_decomposition(
    gateset: Set[OpType],
    allow_swaps: bool,
) -> Callable[[Param, Param, Param], "Circuit"]:
    """Return a function to construct a circuit expressing a TK2 in terms of gates in
    the given gateset, if such a function is available.

    :param gateset: target gate set
    :raises NoAutoRebase: no suitable TK2 decomposition found
    :return: function to decompose TK2 gates
    """
    if allow_swaps:
        for k, fn in _TK2_CIRCS_WIRE_SWAP.items():
            if k in gateset:
                return fn
    else:
        for k, fn in _TK2_CIRCS.items():
            if k in gateset:
                return fn

    raise NoAutoRebase("No known decomposition from TK2 to given gateset")


_TK1_circs: Dict[FrozenSet[OpType], Callable[[Param, Param, Param], "Circuit"]] = {
    frozenset({OpType.TK1}): _library.TK1_to_TK1,
    frozenset({OpType.PhasedX, OpType.Rz}): _library.TK1_to_PhasedXRz,
    frozenset({OpType.Rx, OpType.Rz}): _library.TK1_to_RzRx,
    frozenset({OpType.Ry, OpType.Rx}): _TK1_to_RxRy,
    frozenset({OpType.Rz, OpType.H}): _library.TK1_to_RzH,
    frozenset({OpType.Rz, OpType.SX, OpType.X}): _TK1_to_X_SX_Rz,
    frozenset({OpType.Rz, OpType.SX}): _TK1_to_X_SX_Rz,
    frozenset({OpType.Rz, OpType.SX}): _library.TK1_to_RzSX,
    frozenset({OpType.U3}): _TK1_to_U,
}


def get_TK1_decomposition_function(
    gateset: Set[OpType],
) -> Callable[[Param, Param, Param], "Circuit"]:
    """Return a function for generating TK1 equivalent circuits, which take the
    three TK1 parameters as arguments and return a TK1 equivalent single qubit
    circuit. If no such function is available, raise an error.

    :raises NoAutoRebase: No suitable TK1 decomposition found.
    :return: TK1 decomposition function.
    :rtype: Callable[[Param, Param, Param], "Circuit"]
    """
    subsets = [k for k in _TK1_circs if k.issubset(gateset)]
    if subsets:
        # find the largest available subset
        # as in general more available gates leads to smaller circuits
        matching = max(subsets, key=len)
        return _TK1_circs[matching]
    raise NoAutoRebase("No known decomposition from TK1 to available gateset.")


def auto_rebase_pass(gateset: Set[OpType], allow_swaps: bool = False) -> BasePass:
    """Attempt to generate a rebase pass automatically for the given target
    gateset.

    Checks if there are known existing decompositions
    to target gateset and TK1 to target gateset and uses those to construct a
    custom rebase.
    Raises an error if no known decompositions can be found, in which case try
    using RebaseCustom with your own decompositions.

    In addition to the gate types in ``gateset``, any ``Measure``, ``Reset`` and
    ``Collapse`` operations in the original circuit are retained. Conditional
    operations are also allowed. ``Phase`` gates may also be introduced.

    :param gateset: Set of supported OpTypes, target gate set.
    :type gateset: FrozenSet[OpType]
    :raises NoAutoRebase: No suitable decomposition found.
    :return: Rebase pass.
    :rtype: RebaseCustom
    """
    tk1 = get_TK1_decomposition_function(gateset)
    # if the gateset has CX but not TK2, and implicit wire swaps not allowed:
    # rebase via CX
    if OpType.CX in gateset and OpType.TK2 not in gateset and not allow_swaps:
        return RebaseCustom(gateset, _library.CX(), tk1)
    # in other cases, try to rebase via TK2 first
    try:
        return RebaseCustom(gateset, get_tk2_decomposition(gateset, allow_swaps), tk1)
    except NoAutoRebase:
        pass
    try:
        return RebaseCustom(gateset, get_cx_decomposition(gateset), tk1)
    except NoAutoRebase:
        raise NoAutoRebase(
            "No known decomposition from CX or TK2 to available gateset."
        )


def auto_squash_pass(gateset: Set[OpType]) -> BasePass:
    """Attempt to generate a squash pass automatically for the given target
    single qubit gateset.

    :param gateset: Available single qubit gateset
    :type gateset: Set[OpType]
    :return: Squash to target gateset
    :rtype: SquashCustom
    """
    if {OpType.Rz, OpType.PhasedX} <= gateset:
        return SquashRzPhasedX()

    return SquashCustom(gateset, get_TK1_decomposition_function(gateset))
