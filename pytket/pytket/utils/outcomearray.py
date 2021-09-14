# Copyright 2019-2021 Cambridge Quantum Computing
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

"""`OutcomeArray` class and associated methods."""
import operator
from functools import reduce
from typing import Counter, List, Sequence, Dict, Tuple, Any, cast

import numpy as np


class OutcomeArray(np.ndarray):
    """
    Array of measured outcomes from qubits. Derived class of `numpy.ndarray`.

    Bitwise outcomes are compressed into unsigned 8-bit integers, each
    representing up to 8 qubit measurements. Each row is a repeat measurement.

    :param width: Number of bit entries stored, less than or equal to the bit
        capacity of the array.
    :type width: int
    :param n_outcomes: Number of outcomes stored.
    :type n_outcomes: int
    """

    def __new__(cls, input_array: np.ndarray, width: int):  # type: ignore
        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        obj = np.asarray(input_array).view(cls)
        # add the new attribute to the created instance
        if len(obj.shape) != 2 or obj.dtype != np.uint8:
            raise ValueError(
                "OutcomeArray must be a two dimensional array of dtype uint8."
            )
        bitcapacity = obj.shape[-1] * 8
        if width > bitcapacity:
            raise ValueError(
                f"Width {width} is larger than maxium bitlength of "
                f"array: {bitcapacity}."
            )
        obj._width = width
        # Finally, we must return the newly created object:
        return obj

    def __array_finalize__(self, obj):  # type: ignore
        # see InfoArray.__array_finalize__ for comments
        if obj is None:
            return
        self._width: int = getattr(obj, "_width", None)

    @property
    def width(self) -> int:
        """Number of bit entries stored, less than or equal to the bit capacity of the
        array."""
        return self._width

    @property
    def n_outcomes(self) -> Any:
        """Number of outcomes stored."""
        return self.shape[0]

    def __hash__(self) -> int:
        return hash((self.tobytes(), self.width))

    def __eq__(self, other: "OutcomeArray") -> bool:  # type: ignore
        return bool(np.array_equal(self, other) and self.width == other.width)

    @classmethod
    def from_readouts(cls, readouts: Sequence[Sequence[int]]) -> "OutcomeArray":
        """Create OutcomeArray from a 2D array like object of read-out integers,
        e.g. [[1, 1, 0], [0, 1, 1]]"""
        readouts_ar = np.array(readouts, dtype=int)
        return cls(np.packbits(readouts_ar, axis=-1), readouts_ar.shape[-1])

    def to_readouts(self) -> np.ndarray:
        """Convert OutcomeArray to a 2D array of readouts, each row a separate outcome
        and each column a bit value."""
        return cast(
            np.ndarray, np.asarray(np.unpackbits(self, axis=-1))[..., : self.width]
        )

    def to_readout(self) -> np.ndarray:
        """Convert a singleton to a single readout (1D array)"""
        if self.n_outcomes > 1:
            raise ValueError(f"Not a singleton: {self.n_outcomes} readouts")
        return cast(np.ndarray, self.to_readouts()[0])

    def to_intlist(self, big_endian: bool = True) -> List[int]:
        """Express each outcome as an integer corresponding to the bit values.

        :param big_endian: whether to use big endian encoding (or little endian
            if False), defaults to True
        :type big_endian: bool, optional
        :return: List of integers, each corresponding to an outcome.
        :rtype: List[int]
        """
        if big_endian:
            array = self
        else:
            array = OutcomeArray.from_readouts(np.fliplr(self.to_readouts()))  # type: ignore
        bitcapacity = array.shape[-1] * 8
        intify = lambda bytear: reduce(
            operator.or_, (num << (8 * i) for i, num in enumerate(bytear[::-1])), 0
        ) >> (bitcapacity - array.width)
        intar = np.apply_along_axis(intify, -1, array)  # type: ignore
        return list(intar)

    @classmethod
    def from_ints(
        cls, ints: Sequence[int], width: int, big_endian: bool = True
    ) -> "OutcomeArray":
        """Create OutcomeArray from iterator of integers corresponding to outcomes
         where the bitwise representation of the integer corresponds to the readouts.

        :param ints: Iterable of outcome integers
        :type ints: Iterable[int]
        :param width: Number of qubit measurements
        :type width: int
        :param big_endian: whether to use big endian encoding (or little endian
            if False), defaults to True
        :type big_endian: bool, optional
        :return: OutcomeArray instance
        :rtype: OutcomeArray
        """
        n_ints = len(ints)
        bitstrings = (
            bin(int_val)[2:].zfill(width)[:: (-1) ** (not big_endian)]
            for int_val in ints
        )
        bitar = np.frombuffer(  # type: ignore
            "".join(bitstrings).encode("ascii"), dtype=np.uint8, count=n_ints * width
        ) - ord("0")
        bitar.resize((n_ints, width))
        return cls.from_readouts(bitar)

    def counts(self) -> Counter["OutcomeArray"]:
        """Calculate counts of outcomes in OutcomeArray

        :return: Counter of outcome, number of instances
        :rtype: Counter[OutcomeArray]
        """
        ars, count_vals = np.unique(self, axis=0, return_counts=True)  # type: ignore
        width = self.width
        ars = [OutcomeArray(x[None, :], width) for x in ars]
        return Counter(dict(zip(ars, count_vals)))

    def choose_indices(self, indices: List[int]) -> "OutcomeArray":
        """Permute ordering of bits in outcomes or choose subset of bits.
        e.g. [1, 0, 2] acting on a bitstring of length 4 swaps bit locations 0 & 1,
        leaves 2 in the same place and deletes location 3.

        :param indices: New locations for readout bits.
        :type indices: List[int]
        :return: New array corresponding to given permutation.
        :rtype: OutcomeArray
        """
        return OutcomeArray.from_readouts(self.to_readouts()[..., indices])

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON serializable dictionary representation of the OutcomeArray.

        :return: JSON serializable dictionary
        :rtype: Dict[str, Any]
        """
        return {"width": self.width, "array": self.tolist()}

    @classmethod
    def from_dict(cls, ar_dict: Dict[str, Any]) -> "OutcomeArray":
        """Create an OutcomeArray from JSON serializable dictionary (as created by
        `to_dict`).

        :param dict: Dictionary representation of OutcomeArray.
        :type indices: Dict[str, Any]
        :return: Instance of OutcomeArray
        :rtype: OutcomeArray
        """
        return OutcomeArray(
            np.array(ar_dict["array"], dtype=np.uint8), width=ar_dict["width"]
        )


def readout_counts(
    ctr: Counter[OutcomeArray],
) -> Counter[Tuple[int, ...]]:
    """Convert counts from :py:class:`OutcomeArray` types to tuples of ints."""
    return Counter({tuple(oa.to_readout()): n for oa, n in ctr.items()})
