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

"""Named types for convenience"""
from typing import Union, Sequence

from sympy import Expr

from pytket._tket.circuit import Op
from pytket._tket.unit_id import UnitID, Qubit, Bit

ParamType = Union[float, Expr]
"""Type used for circuit parameters that can either
 be a floating point number or symbolic"""
PhasePolynomialDict = dict[tuple[bool, ...], ParamType]
"""Dict type used to define a phase polynomial.
A dict that maps Bitstrings(tuples) to ParamType"""
PhasePolynomialSequence = Sequence[tuple[Sequence[bool], ParamType]]
"""Sequence type used to define a phase polynomial.
 A sequence of Bitstring(sequence) - ParamType pairs"""
BitstringToOpMap = dict[tuple[bool, ...], Op]
BitstringToOpList = Sequence[tuple[Sequence[bool], Op]]
BitstringToTensoredOpMap = dict[tuple[bool, ...], Sequence[Op]]
BitstringToTensoredOpList = Sequence[tuple[Sequence[bool], Sequence[Op]]]
PermutationMap = dict[tuple[bool, ...], Sequence[bool]]
PermutationList = Sequence[tuple[Sequence[bool], Sequence[bool]]]
UnitIdType = Union[UnitID, Qubit, Bit]
UnitIdMap = dict[UnitIdType, UnitIdType]
RenameUnitsMap = UnitIdMap
