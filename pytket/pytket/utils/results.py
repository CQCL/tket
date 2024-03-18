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

from typing import Dict, List, Tuple, Union

import numpy as np
from pytket.circuit import BasisOrder

StateTuple = Tuple[int, ...]
CountsDict = Dict[StateTuple, Union[int, float]]
KwargTypes = Union[int, float, str, None]


class BitPermuter:
    """Class for permuting the bits in an integer

    Enables inverse permuation and uses caching to speed up common uses.

    """

    def __init__(self, permutation: Tuple[int, ...]):
        """Constructor

        :param permutation: Map from current bit index (big-endian) to its new position,
            encoded as a list.
        :type permutation: Tuple[int, ...]
        :raises ValueError: Input permutation is not valid complete permutation
        of all bits
        """
        if sorted(permutation) != list(range(len(permutation))):
            raise ValueError("Permutation is not a valid complete permutation.")
        self.perm = tuple(permutation)
        self.n_bits = len(self.perm)
        self.int_maps: Tuple[Dict[int, int], Dict[int, int]] = ({}, {})

    def permute(self, val: int, inverse: bool = False) -> int:
        """Return input with bit values permuted.

        :param val: input integer
        :type val: int
        :param inverse: whether to use the inverse permutation, defaults to False
        :type inverse: bool, optional
        :return: permuted integer
        :rtype: int
        """
        perm_map, other_map = self.int_maps[:: (-1) ** inverse]
        if val in perm_map:
            return perm_map[val]

        res = 0
        for source_index, target_index in enumerate(self.perm):
            if inverse:
                target_index, source_index = source_index, target_index
            # if source bit set
            if val & (1 << (self.n_bits - 1 - source_index)):
                # set target bit
                res |= 1 << (self.n_bits - 1 - target_index)

        perm_map[val] = res
        other_map[res] = val
        return res

    def permute_all(self) -> List[int]:
        """Permute all integers within bit-width specified by permutation.

        :return: List of permuted outputs.
        :rtype: List
        """
        return list(map(self.permute, range(1 << self.n_bits)))


def counts_from_shot_table(shot_table: np.ndarray) -> Dict[Tuple[int, ...], int]:
    """Summarises a shot table into a dictionary of counts for each observed outcome.

    :param shot_table: Table of shots from a pytket backend.
    :type shot_table: np.ndarray
    :return: Dictionary mapping observed readouts to the number of times observed.
    :rtype: Dict[Tuple[int, ...], int]
    """
    shot_values, counts = np.unique(shot_table, axis=0, return_counts=True)
    return {tuple(s): c for s, c in zip(shot_values, counts)}


def probs_from_counts(
    counts: Dict[Tuple[int, ...], int]
) -> Dict[Tuple[int, ...], float]:
    """Converts raw counts of observed outcomes into the observed probability
    distribution.

    :param counts: Dictionary mapping observed readouts to the number of times observed.
    :type counts: Dict[Tuple[int, ...], int]
    :return: Probability distribution over observed readouts.
    :rtype: Dict[Tuple[int, ...], float]
    """
    total = sum(counts.values())
    return {outcome: c / total for outcome, c in counts.items()}


def _index_to_readout(
    index: int, width: int, basis: BasisOrder = BasisOrder.ilo
) -> Tuple[int, ...]:
    return tuple(
        (index >> i) & 1 for i in range(width)[:: (-1) ** (basis == BasisOrder.ilo)]
    )


def _reverse_bits_of_index(index: int, width: int) -> int:
    """Reverse bits of a readout/statevector index to change :py:class:`BasisOrder`.
    Values in tket are ILO-BE (2 means [bit0, bit1] == [1, 0]).
    Values in qiskit are DLO-BE (2 means [bit1, bit0] == [1, 0]).
    Note: Since ILO-BE (DLO-BE) is indistinguishable from DLO-LE (ILO-LE), this can also
    be seen as changing the endianness of the value.

    :param n: Value to reverse
    :type n: int
    :param width: Number of bits in bitstring
    :type width: int
    :return: Integer value of reverse bitstring
    :rtype: int
    """
    permuter = BitPermuter(tuple(range(width - 1, -1, -1)))
    return permuter.permute(index)


def _compute_probs_from_state(state: np.ndarray, min_p: float = 1e-10) -> np.ndarray:
    """
    Converts statevector to a probability vector.
    Set probabilities lower than `min_p` to 0.

    :param state: A statevector.
    :type state: np.ndarray
    :param min_p: Minimum probability to include in result
    :type min_p: float
    :return: Probability vector.
    :rtype: np.ndarray
    """
    probs = state.real**2 + state.imag**2
    probs /= sum(probs)
    ignore = probs < min_p
    probs[ignore] = 0
    probs /= sum(probs)
    return probs


def probs_from_state(
    state: np.ndarray, min_p: float = 1e-10
) -> Dict[Tuple[int, ...], float]:
    """
    Converts statevector to the probability distribution over readouts in the
    computational basis. Ignores probabilities lower than `min_p`.

    :param state: Full statevector with big-endian encoding.
    :type state: np.ndarray
    :param min_p: Minimum probability to include in result
    :type min_p: float
    :return: Probability distribution over readouts.
    :rtype: Dict[Tuple[int], float]
    """
    width = get_n_qb_from_statevector(state)
    probs = _compute_probs_from_state(state, min_p)
    return {_index_to_readout(i, width): p for i, p in enumerate(probs) if p != 0}


def int_dist_from_state(state: np.ndarray, min_p: float = 1e-10) -> Dict[int, float]:
    """
    Converts statevector to the probability distribution over
    its indices. Ignores probabilities lower than `min_p`.

    :param state: A statevector.
    :type state: np.ndarray
    :param min_p: Minimum probability to include in result
    :type min_p: float
    :return: Probability distribution over the vector's indices.
    :rtype: Dict[int, float]
    """
    probs = _compute_probs_from_state(state, min_p)
    return {i: p for i, p in enumerate(probs) if p != 0}


def get_n_qb_from_statevector(state: np.ndarray) -> int:
    """Given a statevector, returns the number of qubits described

    :param state: Statevector to inspect
    :type state: np.ndarray
    :raises ValueError: If the dimension of the statevector is not a power of 2
    :return: `n` such that `len(state) == 2 ** n`
    :rtype: int
    """
    n_qb = int(np.log2(state.shape[0]))
    if 2**n_qb != state.shape[0]:
        raise ValueError("Size is not a power of 2")
    return n_qb


def _assert_compatible_state_permutation(
    state: np.ndarray, permutation: Tuple[int, ...]
) -> None:
    """Asserts that a statevector and a permutation list both refer to the same number
    of qubits

    :param state: Statevector
    :type state: np.ndarray
    :param permutation: Permutation of qubit indices, encoded as a list.
    :type permutation: Tuple[int, ...]
    :raises ValueError: [description]
    """
    n_qb = len(permutation)
    if 2**n_qb != state.shape[0]:
        raise ValueError("Invalid permutation: length does not match number of qubits")


def permute_qubits_in_statevector(
    state: np.ndarray, permutation: Tuple[int, ...]
) -> np.ndarray:
    """Rearranges a statevector according to a permutation of the qubit indices.

    >>> # A 3-qubit state:
    >>> state = np.array([0.0, 0.0625, 0.1875, 0.25, 0.375, 0.4375, 0.5, 0.5625])
    >>> permutation = [1, 0, 2] # swap qubits 0 and 1
    >>> # Apply the permutation that swaps indices 2 (="010") and 4 (="100"), and swaps
    >>> # indices 3 (="011") and 5 (="101"):
    >>> permute_qubits_in_statevector(state, permutation)
    array([0.    , 0.0625, 0.375 , 0.4375, 0.1875, 0.25  , 0.5   , 0.5625])

    :param state: Original statevector.
    :type state: np.ndarray
    :param permutation: Map from current qubit index (big-endian) to its new position,
        encoded as a list.
    :type permutation: Tuple[int, ...]
    :return: Updated statevector.
    :rtype: np.ndarray
    """
    _assert_compatible_state_permutation(state, permutation)
    permuter = BitPermuter(permutation)
    return state[permuter.permute_all()]


def permute_basis_indexing(
    matrix: np.ndarray, permutation: Tuple[int, ...]
) -> np.ndarray:
    """Rearranges the first dimensions of an array (statevector or unitary)
     according to a permutation of the bit indices in the binary representation
     of row indices.

    :param matrix: Original unitary matrix
    :type matrix: np.ndarray
    :param permutation: Map from current qubit index (big-endian) to its new position,
        encoded as a list
    :type permutation: Tuple[int, ...]
    :return: Updated unitary matrix
    :rtype: np.ndarray
    """
    _assert_compatible_state_permutation(matrix, permutation)
    permuter = BitPermuter(permutation)

    result: np.ndarray = matrix[permuter.permute_all(), ...]
    return result


def permute_rows_cols_in_unitary(
    matrix: np.ndarray, permutation: Tuple[int, ...]
) -> np.ndarray:
    """Rearranges the rows of a unitary matrix according to a permutation of the qubit
    indices.

    :param matrix: Original unitary matrix
    :type matrix: np.ndarray
    :param permutation: Map from current qubit index (big-endian) to its new position,
        encoded as a list
    :type permutation: Tuple[int, ...]
    :return: Updated unitary matrix
    :rtype: np.ndarray
    """
    _assert_compatible_state_permutation(matrix, permutation)
    permuter = BitPermuter(permutation)
    all_perms = permuter.permute_all()
    permat: np.ndarray = matrix[:, all_perms]
    permat = permat[all_perms, :]
    return permat


def compare_statevectors(first: np.ndarray, second: np.ndarray) -> bool:
    """Check approximate equality up to global phase for statevectors.

    :param first: First statevector.
    :type first: np.ndarray
    :param second: Second statevector.
    :type second: np.ndarray
    :return: Approximate equality.
    :rtype: bool
    """
    return bool(np.isclose(np.abs(np.vdot(first, second)), 1))


def compare_unitaries(first: np.ndarray, second: np.ndarray) -> bool:
    """Check approximate equality up to global phase for unitaries.

    :param first: First unitary.
    :type first: np.ndarray
    :param second: Second unitary.
    :type second: np.ndarray
    :return: Approximate equality.
    :rtype: bool
    """
    conjug_prod = first @ second.conjugate().transpose()
    identity = np.identity(conjug_prod.shape[0], dtype=complex)
    return bool(np.allclose(conjug_prod, identity * conjug_prod[0, 0]))
