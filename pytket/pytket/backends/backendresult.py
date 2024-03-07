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

"""`BackendResult` class and associated methods."""
from typing import (
    Optional,
    Any,
    Sequence,
    Iterable,
    List,
    Tuple,
    Dict,
    Counter,
    NamedTuple,
    Collection,
    Type,
    TypeVar,
    cast,
)
import operator
from functools import reduce
import warnings

import numpy as np

from pytket.circuit import (
    BasisOrder,
    Bit,
    Circuit,
    Qubit,
    UnitID,
    _DEBUG_ZERO_REG_PREFIX,
    _DEBUG_ONE_REG_PREFIX,
)
from pytket.utils.distribution import EmpiricalDistribution, ProbabilityDistribution
from pytket.utils.results import (
    probs_from_state,
    get_n_qb_from_statevector,
    permute_basis_indexing,
    permute_rows_cols_in_unitary,
)
from pytket.utils.outcomearray import OutcomeArray, readout_counts

from .backend_exceptions import InvalidResultType


class StoredResult(NamedTuple):
    """NamedTuple with optional fields for all result types."""

    counts: Optional[Counter[OutcomeArray]] = None
    shots: Optional[OutcomeArray] = None
    state: Optional[np.ndarray] = None
    unitary: Optional[np.ndarray] = None
    density_matrix: Optional[np.ndarray] = None


class BackendResult:
    """Encapsulate generic results from pytket Backend instances.

    In the case of a real quantum device or a shots-based simulator
    a BackendResult will typically be a collection of measurements (shots and counts).

    Results can also be the output of ideal simulations of circuits.
    These can take the form of statevectors, unitary arrays or density matrices.

    :param q_bits: Sequence of qubits.
    :param c_bits: Sequence of classical bits.
    :param counts: The counts in the result.
    :param shots: The shots in the result.
    :param state: The resulting statevector (from a statevector simulation).
    :param unitary: The resulting unitary operator (from a unitary simulation).
    :param density_matrix: The resulting density matrix
        (from a density-matrix simulator).
    :param ppcirc: If provided, classical postprocessing to be applied to all measured
        results (i.e. shots and counts).
    """

    def __init__(
        self,
        *,
        q_bits: Optional[Sequence[Qubit]] = None,
        c_bits: Optional[Sequence[Bit]] = None,
        counts: Optional[Counter[OutcomeArray]] = None,
        shots: Optional[OutcomeArray] = None,
        state: Any = None,
        unitary: Any = None,
        density_matrix: Any = None,
        ppcirc: Optional[Circuit] = None,
    ):
        # deal with mutable defaults
        if q_bits is None:
            q_bits = []
        if c_bits is None:
            c_bits = []

        self._counts = counts
        self._shots = shots

        self._state = state
        self._unitary = unitary
        self._density_matrix = density_matrix

        self._ppcirc = ppcirc

        self.c_bits: Dict[Bit, int] = dict()
        self.q_bits: Dict[Qubit, int] = dict()

        def _process_unitids(
            var: Sequence[UnitID], attr: str, lent: int, uid: Type[UnitID]
        ) -> None:
            if var:
                setattr(self, attr, dict((unit, i) for i, unit in enumerate(var)))
                if lent != len(var):
                    raise ValueError(
                        (
                            f"Length of {attr} ({len(var)}) does not"
                            f" match input data dimensions ({lent})."
                        )
                    )
            else:
                setattr(self, attr, dict((uid(i), i) for i in range(lent)))  # type: ignore

        if self.contains_measured_results:
            _bitlength = 0
            if self._counts is not None:
                if shots is not None:
                    raise ValueError(
                        "Provide either counts or shots, both is not valid."
                    )
                try:
                    _bitlength = next(self._counts.elements()).width
                except StopIteration:
                    _bitlength = len(c_bits)

            if self._shots is not None:
                _bitlength = self._shots.width

            _process_unitids(c_bits, "c_bits", _bitlength, Bit)

        if self.contains_state_results:
            _n_qubits = 0
            if self._unitary is not None:
                _n_qubits = int(np.log2(self._unitary.shape[-1]))
            elif self._state is not None:
                _n_qubits = get_n_qb_from_statevector(self._state)
            elif self._density_matrix is not None:
                _n_qubits = int(np.log2(self._density_matrix.shape[-1]))

            _process_unitids(q_bits, "q_bits", _n_qubits, Qubit)

    def __repr__(self) -> str:
        return (
            "BackendResult(q_bits={s.q_bits},c_bits={s.c_bits},counts={s._counts},"
            "shots={s._shots},state={s._state},unitary={s._unitary},"
            "density_matrix={s._density_matrix})".format(s=self)
        )

    @property
    def contains_measured_results(self) -> bool:
        """Whether measured type results (shots or counts) are stored"""
        return (self._counts is not None) or (self._shots is not None)

    @property
    def contains_state_results(self) -> bool:
        """Whether state type results (state vector or unitary or density_matrix)
        are stored"""
        return (
            (self._state is not None)
            or (self._unitary is not None)
            or (self._density_matrix is not None)
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, BackendResult):
            return NotImplemented
        return (
            self.q_bits == other.q_bits
            and self.c_bits == other.c_bits
            and (
                (self._shots is None and other._shots is None)
                or cast(OutcomeArray, self._shots) == cast(OutcomeArray, other._shots)
            )
            and self._counts == other._counts
            and np.array_equal(self._state, other._state)
            and np.array_equal(self._unitary, other._unitary)
            and np.array_equal(self._density_matrix, other._density_matrix)
        )

    def get_bitlist(self) -> List[Bit]:
        """Return list of Bits in internal storage order.

        :raises AttributeError: BackendResult does not include a Bits list.
        :return: Sorted list of Bits.
        :rtype: List[Bit]
        """
        return _sort_keys_by_val(self.c_bits)

    def get_qbitlist(self) -> List[Qubit]:
        """Return list of Qubits in internal storage order.

        :raises AttributeError: BackendResult does not include a Qubits list.
        :return: Sorted list of Qubits.
        :rtype: List[Qubit]
        """

        return _sort_keys_by_val(self.q_bits)

    def _get_measured_res(
        self, bits: Sequence[Bit], ppcirc: Optional[Circuit] = None
    ) -> StoredResult:
        vals: Dict[str, Any] = {}

        if not self.contains_measured_results:
            raise InvalidResultType("shots/counts")

        if self._ppcirc is not None:
            if ppcirc is None:
                ppcirc = self._ppcirc
            else:
                raise ValueError("Postprocessing circuit already provided.")

        try:
            chosen_readouts = [self.c_bits[bit] for bit in bits]
        except KeyError:
            raise ValueError("Requested Bit not in result.")

        if self._counts is not None:
            if ppcirc is not None:
                # Modify self._counts:
                new_counts: Counter[OutcomeArray] = Counter()
                for oa, n in self._counts.items():
                    readout = oa.to_readout()
                    values = {bit: bool(readout[i]) for bit, i in self.c_bits.items()}
                    new_values = ppcirc._classical_eval(values)
                    new_oa = OutcomeArray.from_readouts(
                        [[int(new_values[bit]) for bit in bits]]
                    )
                    new_counts[new_oa] += n
            else:
                new_counts = self._counts
            vals["counts"] = reduce(
                operator.add,
                (
                    Counter({outcome.choose_indices(chosen_readouts): count})
                    for outcome, count in new_counts.items()
                ),
                Counter(),
            )
        if self._shots is not None:
            if ppcirc is not None:
                # Modify self._shots:
                readouts = self._shots.to_readouts()
                new_readouts = []
                for i in range(self._shots.n_outcomes):
                    readout = readouts[i, :]
                    values = {bit: bool(readout[i]) for bit, i in self.c_bits.items()}
                    new_values = ppcirc._classical_eval(values)
                    new_readout = [0] * self._shots.width
                    for bit, i in self.c_bits.items():
                        if new_values[bit]:
                            new_readout[i] = 1
                    new_readouts.append(new_readout)
                new_shots = OutcomeArray.from_readouts(new_readouts)
            else:
                new_shots = self._shots
            vals["shots"] = new_shots.choose_indices(chosen_readouts)

        return StoredResult(**vals)

    def _permute_statearray_qb_labels(
        self,
        array: np.ndarray,
        relabling_map: Dict[Qubit, Qubit],
    ) -> np.ndarray:
        """Permute statevector/unitary according to a relabelling of Qubits.

        :param array: The statevector or unitary
        :type array: np.ndarray
        :param relabling_map: Map from original Qubits to new.
        :type relabling_map: Dict[Qubit, Qubit]
        :return: Permuted array.
        :rtype: np.ndarray
        """
        original_labeling: Sequence["Qubit"] = self.get_qbitlist()
        n_labels = len(original_labeling)
        permutation = [0] * n_labels
        for i, orig_qb in enumerate(original_labeling):
            permutation[i] = original_labeling.index(relabling_map[orig_qb])
        if permutation == list(range(n_labels)):
            # Optimization: nothing to do; return original array.
            return array
        permuter = (
            permute_basis_indexing
            if len(array.shape) == 1
            else permute_rows_cols_in_unitary
        )
        return permuter(array, tuple(permutation))

    def _get_state_res(self, qubits: Sequence[Qubit]) -> StoredResult:
        vals: Dict[str, Any] = {}
        if not self.contains_state_results:
            raise InvalidResultType("state/unitary/density_matrix")

        if not _check_permuted_sequence(qubits, self.q_bits):
            raise ValueError(
                "For state/unitary/density_matrix results only a permutation of"
                " all qubits can be requested."
            )
        qb_mapping = {selfqb: qubits[index] for selfqb, index in self.q_bits.items()}
        if self._state is not None:
            vals["state"] = self._permute_statearray_qb_labels(self._state, qb_mapping)
        if self._unitary is not None:
            vals["unitary"] = self._permute_statearray_qb_labels(
                self._unitary, qb_mapping
            )
        if self._density_matrix is not None:
            vals["density_matrix"] = self._permute_statearray_qb_labels(
                self._density_matrix, qb_mapping
            )
        return StoredResult(**vals)

    def get_result(
        self,
        request_ids: Optional[Sequence[UnitID]] = None,
        basis: BasisOrder = BasisOrder.ilo,
        ppcirc: Optional[Circuit] = None,
    ) -> StoredResult:
        """Retrieve all results, optionally according to a specified UnitID ordering
         or subset.

        :param request_ids: Ordered set of either Qubits or Bits for which to
            retrieve results, defaults to None in which case all results are returned.
            For statevector/unitary/density_matrix results some permutation of all
            qubits must be requested.
            For measured results (shots/counts), some subset of the relevant bits must
            be requested.
        :type request_ids: Optional[Sequence[UnitID]], optional
        :param basis: Toggle between ILO (increasing lexicographic order of bit ids) and
            DLO (decreasing lexicographic order) for column ordering if request_ids is
            None. Defaults to BasisOrder.ilo.
        :param ppcirc: Classical post-processing circuit to apply to measured results
        :raises ValueError: Requested UnitIds (request_ids) contain a mixture of qubits
            and bits.
        :raises RuntimeError: Classical bits not set.
        :raises ValueError: Requested (Qu)Bit not in result.
        :raises RuntimeError: "Qubits not set."
        :raises ValueError: For state/unitary/density_matrix results only a permutation
            of all qubits can be requested.
        :return: All stored results corresponding to requested IDs.
        :rtype: StoredResult
        """
        if request_ids is None:
            if self.contains_measured_results:
                request_ids = sorted(
                    self.c_bits.keys(), reverse=(basis == BasisOrder.dlo)
                )
            elif self.contains_state_results:
                request_ids = sorted(
                    self.q_bits.keys(), reverse=(basis == BasisOrder.dlo)
                )
            else:
                raise InvalidResultType("No results stored.")

        if all(isinstance(i, Bit) for i in request_ids):
            return self._get_measured_res(request_ids, ppcirc)  # type: ignore

        if all(isinstance(i, Qubit) for i in request_ids):
            return self._get_state_res(request_ids)  # type: ignore

        raise ValueError(
            "Requested UnitIds (request_ids) contain a mixture of qubits and bits."
        )

    def get_shots(
        self,
        cbits: Optional[Sequence[Bit]] = None,
        basis: BasisOrder = BasisOrder.ilo,
        ppcirc: Optional[Circuit] = None,
    ) -> np.ndarray:
        """Return shots if available.

        :param cbits: ordered subset of Bits, returns all results by default, defaults
         to None
        :type cbits: Optional[Sequence[Bit]], optional
        :param basis: Toggle between ILO (increasing lexicographic order of bit ids) and
            DLO (decreasing lexicographic order) for column ordering if cbits is None.
            Defaults to BasisOrder.ilo.
        :param ppcirc: Classical post-processing circuit to apply to measured results
        :raises InvalidResultType: Shot results are not available
        :return: 2D array of readouts, each row a separate outcome and each column a
         bit value.
        :rtype: np.ndarray

        The order of the columns follows the order of `cbits`, if provided.
        """
        if cbits is None:
            cbits = sorted(self.c_bits.keys(), reverse=(basis == BasisOrder.dlo))
        res = self.get_result(cbits, ppcirc=ppcirc)
        if res.shots is not None:
            return res.shots.to_readouts()
        raise InvalidResultType("shots")

    def get_counts(
        self,
        cbits: Optional[Sequence[Bit]] = None,
        basis: BasisOrder = BasisOrder.ilo,
        ppcirc: Optional[Circuit] = None,
    ) -> Counter[Tuple[int, ...]]:
        """Return counts of outcomes if available.

        :param cbits: ordered subset of Bits, returns all results by default, defaults
         to None
        :type cbits: Optional[Sequence[Bit]], optional
        :param basis: Toggle between ILO (increasing lexicographic order of bit ids) and
            DLO (decreasing lexicographic order) for column ordering if cbits is None.
            Defaults to BasisOrder.ilo.
        :param ppcirc: Classical post-processing circuit to apply to measured results
        :raises InvalidResultType: Counts are not available
        :return: Counts of outcomes
        :rtype: Counter[Tuple(int)]
        """
        if cbits is None:
            cbits = sorted(self.c_bits.keys(), reverse=(basis == BasisOrder.dlo))
        res = self.get_result(cbits, ppcirc=ppcirc)
        if res.counts is not None:
            return readout_counts(res.counts)
        if res.shots is not None:
            return readout_counts(res.shots.counts())
        raise InvalidResultType("counts")

    def get_state(
        self,
        qbits: Optional[Sequence[Qubit]] = None,
        basis: BasisOrder = BasisOrder.ilo,
    ) -> np.ndarray:
        """Return statevector if available.

        :param qbits: permutation of Qubits, defaults to None
        :type qbits: Optional[Sequence[Qubit]], optional
        :param basis: Toggle between ILO (increasing lexicographic order of qubit ids)
            and DLO (decreasing lexicographic order) for column ordering if qbits is
            None. Defaults to BasisOrder.ilo.
        :raises InvalidResultType: Statevector not available
        :return: Statevector, (complex 1-D numpy array)
        :rtype: np.ndarray
        """
        if qbits is None:
            qbits = sorted(self.q_bits.keys(), reverse=(basis == BasisOrder.dlo))
        res = self.get_result(qbits)
        if res.state is not None:
            return res.state
        if res.unitary is not None:
            state: np.ndarray = res.unitary[:, 0]
            return state
        raise InvalidResultType("state")

    def get_unitary(
        self,
        qbits: Optional[Sequence[Qubit]] = None,
        basis: BasisOrder = BasisOrder.ilo,
    ) -> np.ndarray:
        """Return unitary if available.

        :param qbits: permutation of Qubits, defaults to None
        :type qbits: Optional[Sequence[Qubit]], optional
        :param basis: Toggle between ILO (increasing lexicographic order of qubit ids)
            and DLO (decreasing lexicographic order) for column ordering if qbits is
            None. Defaults to BasisOrder.ilo.
        :raises InvalidResultType: Statevector not available
        :return: Unitary, (complex 2-D numpy array)
        :rtype: np.ndarray
        """
        if qbits is None:
            qbits = sorted(self.q_bits.keys(), reverse=(basis == BasisOrder.dlo))
        res = self.get_result(qbits)
        if res.unitary is not None:
            return res.unitary
        raise InvalidResultType("unitary")

    def get_density_matrix(
        self,
        qbits: Optional[Sequence[Qubit]] = None,
        basis: BasisOrder = BasisOrder.ilo,
    ) -> np.ndarray:
        """Return density_matrix if available.

        :param qbits: permutation of Qubits, defaults to None
        :type qbits: Optional[Sequence[Qubit]], optional
        :param basis: Toggle between ILO (increasing lexicographic order of qubit ids)
            and DLO (decreasing lexicographic order) for column ordering if qbits is
            None. Defaults to BasisOrder.ilo.
        :raises InvalidResultType: Statevector not available
        :return: density_matrix, (complex 2-D numpy array)
        :rtype: np.ndarray
        """
        if qbits is None:
            qbits = sorted(self.q_bits.keys(), reverse=(basis == BasisOrder.dlo))
        res = self.get_result(qbits)
        if res.density_matrix is not None:
            return res.density_matrix
        raise InvalidResultType("density_matrix")

    def get_distribution(
        self, units: Optional[Sequence[UnitID]] = None
    ) -> Dict[Tuple[int, ...], float]:
        """Calculate an exact or approximate probability distribution over outcomes.

        If the exact statevector is known, the exact probability distribution is
        returned. Otherwise, if measured results are available the distribution
        is estimated from these results.

        This method is deprecated. Please use :py:meth:`get_empirical_distribution` or
        :py:meth:`get_probability_distribution` instead.

        :param units: Optionally provide the Qubits or Bits
            to marginalise the distribution over, defaults to None
        :type units: Optional[Sequence[UnitID]], optional
        :return: A distribution as a map from bitstring to probability.
        :rtype: Dict[Tuple[int, ...], float]
        """
        warnings.warn(
            "The `BackendResult.get_distribution()` method is deprecated: "
            "please use `get_empirical_distribution()` or "
            "`get_probability_distribution()` instead.",
            DeprecationWarning,
        )
        try:
            state = self.get_state(units)  # type: ignore
            return probs_from_state(state)
        except InvalidResultType:
            counts = self.get_counts(units)  # type: ignore
            total = sum(counts.values())
            dist = {outcome: count / total for outcome, count in counts.items()}
            return dist

    def get_empirical_distribution(
        self, bits: Optional[Sequence[Bit]] = None
    ) -> EmpiricalDistribution[Tuple[int, ...]]:
        """Convert to a :py:class:`pytket.utils.distribution.EmpiricalDistribution`
        where the observations are sequences of 0s and 1s.

        :param bits: Optionally provide the :py:class:`Bit` s over which to
            marginalize the distribution.
        :return: A distribution where the observations are sequences of 0s and 1s.
        """
        if not self.contains_measured_results:
            raise InvalidResultType(
                "Empirical distribution only available for measured result types."
            )
        return EmpiricalDistribution(self.get_counts(bits))

    def get_probability_distribution(
        self, qubits: Optional[Sequence[Qubit]] = None, min_p: float = 0.0
    ) -> ProbabilityDistribution[Tuple[int, ...]]:
        """Convert to a :py:class:`pytket.utils.distribution.ProbabilityDistribution`
        where the possible outcomes are sequences of 0s and 1s.

        :param qubits: Optionally provide the :py:class:`Qubit` s over which to
            marginalize the distribution.
        :param min_p: Optional probability below which to ignore values (for
            example to avoid spurious values due to rounding errors in
            statevector computations). Default 0.
        :return: A distribution where the possible outcomes are tuples of 0s and 1s.
        """
        if not self.contains_state_results:
            raise InvalidResultType(
                "Probability distribution only available for statevector result types."
            )
        state = self.get_state(qubits)
        return ProbabilityDistribution(probs_from_state(state), min_p=min_p)

    def get_debug_info(self) -> Dict[str, float]:
        """Calculate the success rate of each assertion averaged across shots.

        Each assertion in pytket is decomposed into a sequence of transformations
        and measurements. An assertion is successful if and only if all its associated
        measurements yield the correct results.

        :return: The debug results as a map from assertion to average success rate.
        :rtype: Dict[str, float]
        """
        _tket_debug_zero_prefix = _DEBUG_ZERO_REG_PREFIX + "_"
        _tket_debug_one_prefix = _DEBUG_ONE_REG_PREFIX + "_"
        debug_bit_dict: Dict[str, Dict[str, Any]] = {}
        for bit in self.c_bits:
            if bit.reg_name.startswith(_tket_debug_zero_prefix):
                expectation = 0
                assertion_name = bit.reg_name.split(_tket_debug_zero_prefix, 1)[1]
            elif bit.reg_name.startswith(_tket_debug_one_prefix):
                expectation = 1
                assertion_name = bit.reg_name.split(_tket_debug_one_prefix, 1)[1]
            else:
                continue
            if assertion_name not in debug_bit_dict:
                debug_bit_dict[assertion_name] = {"bits": [], "expectations": []}
            debug_bit_dict[assertion_name]["bits"].append(bit)
            debug_bit_dict[assertion_name]["expectations"].append(expectation)

        debug_result_dict: Dict[str, float] = {}
        for assertion_name, bits_info in debug_bit_dict.items():
            counts = self.get_counts(bits_info["bits"])
            debug_result_dict[assertion_name] = counts[
                tuple(bits_info["expectations"])
            ] / sum(counts.values())
        return debug_result_dict

    def to_dict(self) -> Dict[str, Any]:
        """Generate a dictionary serialized representation of BackendResult,
         suitable for writing to JSON.

        :return: JSON serializable dictionary.
        :rtype: Dict[str, Any]
        """
        outdict: Dict[str, Any] = dict()
        outdict["qubits"] = [q.to_list() for q in self.get_qbitlist()]
        outdict["bits"] = [c.to_list() for c in self.get_bitlist()]
        if self._shots is not None:
            outdict["shots"] = self._shots.to_dict()
        if self._counts is not None:
            outdict["counts"] = [
                {"outcome": oc.to_dict(), "count": count}
                for oc, count in self._counts.items()
            ]

        if self._state is not None:
            outdict["state"] = _complex_ar_to_dict(self._state)
        if self._unitary is not None:
            outdict["unitary"] = _complex_ar_to_dict(self._unitary)
        if self._density_matrix is not None:
            outdict["density_matrix"] = _complex_ar_to_dict(self._density_matrix)

        return outdict

    @classmethod
    def from_dict(cls, res_dict: Dict[str, Any]) -> "BackendResult":
        """Construct BackendResult object from JSON serializable dictionary
         representation, as generated by BackendResult.to_dict.

        :return: Instance of BackendResult constructed from dictionary.
        :rtype: BackendResult
        """
        init_dict = dict.fromkeys(
            (
                "q_bits",
                "c_bits",
                "shots",
                "counts",
                "state",
                "unitary",
                "density_matrix",
            )
        )

        if "qubits" in res_dict:
            init_dict["q_bits"] = [Qubit.from_list(tup) for tup in res_dict["qubits"]]
        if "bits" in res_dict:
            init_dict["c_bits"] = [Bit.from_list(tup) for tup in res_dict["bits"]]
        if "shots" in res_dict:
            init_dict["shots"] = OutcomeArray.from_dict(res_dict["shots"])
        if "counts" in res_dict:
            init_dict["counts"] = Counter(
                {
                    OutcomeArray.from_dict(elem["outcome"]): elem["count"]
                    for elem in res_dict["counts"]
                }
            )
        if "state" in res_dict:
            init_dict["state"] = _complex_ar_from_dict(res_dict["state"])
        if "unitary" in res_dict:
            init_dict["unitary"] = _complex_ar_from_dict(res_dict["unitary"])
        if "density_matrix" in res_dict:
            init_dict["density_matrix"] = _complex_ar_from_dict(
                res_dict["density_matrix"]
            )

        return BackendResult(**init_dict)


T = TypeVar("T")


def _sort_keys_by_val(dic: Dict[T, int]) -> List[T]:
    if not dic:
        return []
    vals, _ = zip(*sorted(dic.items(), key=lambda x: x[1]))
    return list(cast(Iterable[T], vals))


def _check_permuted_sequence(first: Collection[Any], second: Collection[Any]) -> bool:
    return len(first) == len(second) and set(first) == set(second)


def _complex_ar_to_dict(ar: np.ndarray) -> Dict[str, List]:
    """Dictionary of real, imaginary parts of complex array, each in list form."""
    return {"real": ar.real.tolist(), "imag": ar.imag.tolist()}


def _complex_ar_from_dict(dic: Dict[str, List]) -> np.ndarray:
    """Construct complex array from dictionary of real and imaginary parts"""

    out = np.array(dic["real"], dtype=complex)
    out.imag = np.array(dic["imag"], dtype=float)
    return out
