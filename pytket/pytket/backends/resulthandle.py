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

"""ResultHandle class
"""

from typing import Tuple, Type, Union, Iterator, overload
from ast import literal_eval
from collections.abc import Sequence

# mypy doesn't think you can pass the tuple to Union
BasicHashType = Union[int, float, complex, str, bool, bytes]
_ResultIdTuple = Tuple[
    Union[Type[int], Type[float], Type[complex], Type[str], Type[bool], Type[bytes]],
    ...,
]


class ResultHandle(Sequence):
    """Object to store multidimensional identifiers for a circuit sent to a backend for
    execution.

    Initialisation arguments must be hashable basic types.

    Note that a `ResultHandle` may be either persistent or transient, depending on the
    backend: consult the :py:attr:`pytket.backends.Backend.persistent_handles` property
    to determine this.
    """

    def __init__(self, *args: BasicHashType):
        self._identifiers = tuple(args)

    @classmethod
    def from_str(cls, string: str) -> "ResultHandle":
        """Construct ResultHandle from string (output from str())

        :raises ValueError: If string format is invalid
        :return: Instance of ResultHandle
        :rtype: ResultHandle
        """
        try:
            evaltuple = literal_eval(string)  # will raise ValueError if failed
            if (not isinstance(evaltuple, tuple)) or (
                not all(
                    isinstance(arg, (int, float, complex, str, bool, bytes))
                    for arg in evaltuple
                )
            ):
                raise ValueError  # type check failed
            return cls(*evaltuple)
        except ValueError:
            raise ValueError("ResultHandle string format invalid.")

    def __hash__(self) -> int:
        return hash(self._identifiers)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ResultHandle):
            return NotImplemented
        return bool(self._identifiers == other._identifiers)

    def __str__(self) -> str:
        return str(self._identifiers)

    def __repr__(self) -> str:
        return "ResultHandle" + repr(self._identifiers)

    def __iter__(self) -> Iterator:
        return iter(self._identifiers)

    def __len__(self) -> int:
        return len(self._identifiers)

    @overload
    def __getitem__(self, key: int) -> BasicHashType: ...

    @overload
    def __getitem__(self, key: slice) -> Tuple[BasicHashType, ...]: ...

    def __getitem__(
        self, key: Union[int, slice]
    ) -> Union[BasicHashType, Tuple[BasicHashType, ...]]:
        # weird logic required to make mypy happy, can't just
        # return self._identifiers[key]
        if isinstance(key, slice):
            return self._identifiers[key]
        return self._identifiers[key]
