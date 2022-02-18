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

from typing import Set, Callable, Dict, FrozenSet
from pytket.circuit import Circuit, OpType  # type: ignore
from pytket._tket.circuit import _library  # type: ignore
from pytket.passes import RebaseCustom, SquashCustom  # type: ignore

from ._decompositions import Param, _TK1_to_X_SX_Rz, _TK1_to_RxRy, _TK1_to_U


class NoAutoRebase(Exception):
    """Automatic rebase could not be found."""


_CX_CIRCS: Dict[OpType, Callable[[], "Circuit"]] = {
    OpType.CX: _library._CX,
    OpType.ZZMax: _library._CX_using_ZZMax,
    OpType.XXPhase: _library._CX_using_XXPhase_0,
    OpType.ECR: _library._CX_using_ECR,
    OpType.CZ: _library._H_CZ_H,
}


def get_cx_decomposition(gateset: Set[OpType]) -> Circuit:
    """Return a Circuit expressing a CX in terms of a two qubit gate in the
    gateset if one is available, raise an error otherwise.

    :param gateset: Target gate set.
    :type gateset: Set[OpType]
    :raises NoAutoRebase: No suitable CX decomposition found.
    :return: Decomposuition circuit.
    :rtype: Circuit
    """
    if any((matching := k) in gateset for k in _CX_CIRCS):
        return _CX_CIRCS[matching]()
    raise NoAutoRebase("No known decomposition from CX to available gateset.")


_TK1_circs: Dict[FrozenSet[OpType], Callable[[Param, Param, Param], "Circuit"]] = {
    frozenset({OpType.TK1}): _library._TK1_to_TK1,
    frozenset({OpType.PhasedX, OpType.Rz}): _library._TK1_to_PhasedXRz,
    frozenset({OpType.Rx, OpType.Rz}): _library._TK1_to_RzRx,
    frozenset({OpType.Ry, OpType.Rx}): _TK1_to_RxRy,
    frozenset({OpType.Rz, OpType.H}): _library._TK1_to_RzH,
    frozenset({OpType.Rz, OpType.SX, OpType.X}): _TK1_to_X_SX_Rz,
    frozenset({OpType.Rz, OpType.SX}): _TK1_to_X_SX_Rz,
    frozenset({OpType.Rz, OpType.SX}): _library._TK1_to_RzSX,
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


def auto_rebase_pass(gateset: Set[OpType]) -> RebaseCustom:
    """Attempt to generate a rebase pass automatically for the given target
    gateset.

    Checks if there are known existing decompositions from CX
    to target gateset and TK1 to target gateset and uses those to construct a
    custom rebase.
    Raises an error if no known decompositions can be found, in which case try
    using RebaseCustom with your own decompositions.

    :param gateset: Set of supported OpTypes, target gate set.
    :type gateset: FrozenSet[OpType]
    :raises NoAutoRebase: No suitable CX or TK1 decomposition found.
    :return: Rebase pass.
    :rtype: RebaseCustom
    """
    return RebaseCustom(
        gateset, get_cx_decomposition(gateset), get_TK1_decomposition_function(gateset)
    )


def auto_squash_pass(gateset: Set[OpType]) -> SquashCustom:
    """Attempt to generate a squash pass automatically for the given target
    single qubit gateset.

    :param gateset: Available single qubit gateset
    :type gateset: Set[OpType]
    :return: Squash to target gateset
    :rtype: SquashCustom
    """
    return SquashCustom(gateset, get_TK1_decomposition_function(gateset))
