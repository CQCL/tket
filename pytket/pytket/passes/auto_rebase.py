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

import warnings
from typing import Set
from pytket.circuit import OpType
from pytket.passes import AutoRebase, AutoSquash, BasePass


class NoAutoRebase(Exception):
    """Automatic rebase could not be found."""

    def __init__(self, *args: object, **kwargs: dict) -> None:
        warnings.warn(
            "The `NoAutoRebase` exception is deprecated along with "
            "`auto_rebase_pass()` and `auto_squash_pass()`. "
            "`AutoRebase` and `AutoSquash` will raise `RuntimeError` instead.",
            DeprecationWarning,
        )
        super().__init__(*args, **kwargs)


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
    :rtype: AutoRebase

    .. deprecated:: 1.29.0
       Will be removed after pytket 1.33, please use `AutoRebase` instead.
    """
    warnings.warn(
        "The `auto_rebase_pass()` method is deprecated and will be removed after "
        "pytket 1.33, please use `AutoRebase` instead.",
        DeprecationWarning,
    )
    try:
        return AutoRebase(gateset, allow_swaps)
    except RuntimeError:
        raise NoAutoRebase(
            "No known decomposition from TK1, CX or TK2 to available gateset."
        )


def auto_squash_pass(gateset: Set[OpType]) -> BasePass:
    """Attempt to generate a squash pass automatically for the given target
    single qubit gateset.

    :param gateset: Available single qubit gateset
    :type gateset: Set[OpType]
    :return: Squash to target gateset
    :rtype: AutoSquash

    .. deprecated:: 1.29.0
        Will be removed after pytket 1.33, please use `AutoSquash` instead.
    """
    warnings.warn(
        "The `auto_squash_pass()` method is deprecated and will be removed after "
        "pytket 1.33, please use `AutoSquash` instead.",
        DeprecationWarning,
    )
    try:
        return AutoSquash(gateset)
    except RuntimeError:
        raise NoAutoRebase("No known decomposition from TK1 to available gateset.")
