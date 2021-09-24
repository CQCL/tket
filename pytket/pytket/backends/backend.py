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

""" Abstract base class for all Backend encapsulations."""
import warnings
from abc import ABC, abstractmethod
from typing import (
    Dict,
    Iterable,
    List,
    Optional,
    Sequence,
    Tuple,
    Union,
    Any,
    cast,
    overload,
)
from importlib import import_module
from types import ModuleType

import numpy as np
from typing_extensions import Literal

from pytket.circuit import BasisOrder, Bit, Circuit, OpType  # type: ignore
from pytket.passes import BasePass  # type: ignore
from pytket.predicates import Predicate  # type: ignore
from pytket.utils.outcomearray import OutcomeArray
from pytket.utils.results import KwargTypes, counts_from_shot_table

from .backend_exceptions import (
    CircuitNotValidError,
    InvalidResultType,
    CircuitNotRunError,
)
from .backendinfo import BackendInfo
from .backendresult import BackendResult
from .resulthandle import ResultHandle, _ResultIdTuple
from .status import CircuitStatus

ResultCache = Dict[str, Any]


class Backend(ABC):
    """
    This abstract class defines the structure of a backend as something that
    can run quantum circuits and produce output as at least one of shots,
    counts, state, or unitary
    """

    _supports_shots = False
    _supports_counts = False
    _supports_state = False
    _supports_unitary = False
    _supports_density_matrix = False
    _supports_expectation = False
    _expectation_allows_nonhermitian = False
    _supports_contextual_optimisation = False
    _persistent_handles = False

    def __init__(self) -> None:
        self._cache: Dict[ResultHandle, ResultCache] = {}

    @staticmethod
    def empty_result(circuit: Circuit, n_shots: int) -> BackendResult:
        n_bits = len(circuit.bits)
        empty_readouts = [[0] * n_bits for _ in range(n_shots)]
        shots = OutcomeArray.from_readouts(empty_readouts)
        c_bits = [Bit(index) for index in range(n_bits)]
        return BackendResult(shots=shots, c_bits=c_bits)

    @property
    @abstractmethod
    def required_predicates(self) -> List[Predicate]:
        """
        The minimum set of predicates that a circuit must satisfy before it can
        be successfully run on this backend.

        :return: Required predicates.
        :rtype: List[Predicate]
        """
        ...

    def valid_circuit(self, circuit: Circuit) -> bool:
        """
        Checks that the circuit satisfies all of required_predicates.

        :param circuit: The circuit to check.
        :type circuit: Circuit
        :return: Whether or not all of required_predicates are satisfied.
        :rtype: bool
        """
        return all(pred.verify(circuit) for pred in self.required_predicates)

    def _check_all_circuits(
        self, circuits: Iterable[Circuit], nomeasure_warn: Optional[bool] = None
    ) -> bool:

        if nomeasure_warn is None:
            nomeasure_warn = not (
                self._supports_state
                or self._supports_unitary
                or self._supports_density_matrix
                or self._supports_expectation
            )
        for i, circ in enumerate(circuits):
            errors = (
                CircuitNotValidError(i, repr(pred))
                for pred in self.required_predicates
                if not pred.verify(circ)
            )
            for error in errors:
                raise error
            if nomeasure_warn:
                if circ.n_gates_of_type(OpType.Measure) < 1:
                    warnings.warn(
                        f"Circuit with index {i} in submitted does not contain a "
                        "measure operation."
                    )
        return True

    @abstractmethod
    def default_compilation_pass(self, optimisation_level: int = 1) -> BasePass:
        """
        A suggested compilation pass that will guarantee the resulting circuit
        will be suitable to run on this backend with as few preconditions as
        possible.

        :param optimisation_level: The level of optimisation to perform during
            compilation. Level 0 just solves the device constraints without
            optimising. Level 1 additionally performs some light optimisations.
            Level 2 adds more intensive optimisations that can increase compilation
            time for large circuits. Defaults to 1.
        :type optimisation_level: int, optional
        :return: Compilation pass guaranteeing required predicates.
        :rtype: BasePass
        """
        ...

    def compile_circuit(self, circuit: Circuit, optimisation_level: int = 1) -> None:
        """Apply the default_compilation_pass to a circuit in place.

        As well as applying a degree of optimisation (controlled by the
        `optimisation_level` parameter), this method tries to ensure that the circuit
        can be run on the backend (i.e. successfully passed to
        :py:meth:`process_circuits`), for example by rebasing to the supported gate set,
        or routing to match the connectivity of the device. However, this is not always
        possible, for example if the circuit contains classical operations that are not
        supported by the backend. You may use :py:meth:`valid_circuit` to check whether
        the circuit meets the backend's requirements after compilation. This validity
        check is included in :py:meth:`process_circuits` by default, before any circuits
        are submitted to the backend.

        If the validity check fails, you can obtain more information about the failure
        by iterating through the predicates in the `required_predicates` property of the
        backend, and running the :py:meth:`verify` method on each in turn with your
        circuit.

        :param circuit: The circuit to compile.
        :type circuit: Circuit
        :param optimisation_level: The level of optimisation to perform during
            compilation. Level 0 just solves the device constraints without
            optimising. Level 1 additionally performs some light optimisations.
            Level 2 adds more intensive optimisations that can increase compilation
            time for large circuits. Defaults to 1.
        :type optimisation_level: int, optional
        """
        warnings.warn(
            "compile_circuit is deprecated and will be removed in a future pytket.",
            DeprecationWarning,
        )
        self.default_compilation_pass(optimisation_level).apply(circuit)

    def get_compiled_circuit(
        self, circuit: Circuit, optimisation_level: int = 1
    ) -> Circuit:
        """
        Return a single circuit compiled with default_compilation_pass. See
        :py:meth:`Backend.get_compiled_circuits`.
        """
        return_circuit = circuit.copy()
        self.default_compilation_pass(optimisation_level).apply(return_circuit)
        return return_circuit

    def get_compiled_circuits(
        self, circuits: Sequence[Circuit], optimisation_level: int = 1
    ) -> List[Circuit]:
        """Compile a sequence of circuits with default_compilation_pass
        and return the list of compiled circuits (does not act in place).

        As well as applying a degree of optimisation (controlled by the
        `optimisation_level` parameter), this method tries to ensure that the circuits
        can be run on the backend (i.e. successfully passed to
        :py:meth:`process_circuits`), for example by rebasing to the supported gate set,
        or routing to match the connectivity of the device. However, this is not always
        possible, for example if the circuit contains classical operations that are not
        supported by the backend. You may use :py:meth:`valid_circuit` to check whether
        the circuit meets the backend's requirements after compilation. This validity
        check is included in :py:meth:`process_circuits` by default, before any circuits
        are submitted to the backend.

        If the validity check fails, you can obtain more information about the failure
        by iterating through the predicates in the `required_predicates` property of the
        backend, and running the :py:meth:`verify` method on each in turn with your
        circuit.

        :param circuits: The circuits to compile.
        :type circuit: Sequence[Circuit]
        :param optimisation_level: The level of optimisation to perform during
            compilation. Level 0 just solves the device constraints without
            optimising. Level 1 additionally performs some light optimisations.
            Level 2 adds more intensive optimisations that can increase compilation
            time for large circuits. Defaults to 1.
        :type optimisation_level: int, optional
        :return: Compiled circuits.
        :rtype: List[Circuit]
        """
        return [self.get_compiled_circuit(c, optimisation_level) for c in circuits]

    @property
    @abstractmethod
    def _result_id_type(self) -> _ResultIdTuple:
        """Identifier type signature for ResultHandle for this backend.

        :return: Type signature (tuple of hashable types)
        :rtype: _ResultIdTuple
        """
        ...

    def _check_handle_type(self, reshandle: ResultHandle) -> None:
        """Check a result handle is valid for this backend, raises TypeError if not.

        :param reshandle: Handle to check
        :type reshandle: ResultHandle
        :raises TypeError: Types of handle identifiers don't match those of backend.
        """
        if (len(reshandle) != len(self._result_id_type)) or not all(
            isinstance(idval, ty) for idval, ty in zip(reshandle, self._result_id_type)
        ):
            raise TypeError(
                "{0!r} does not match expected identifier types {1}".format(
                    reshandle, self._result_id_type
                )
            )

    @property
    def persistent_handles(self) -> bool:
        """
        Whether the backend produces `ResultHandle` objects that can be reused with
        other instances of the backend class.
        """
        return self._persistent_handles

    def process_circuit(
        self,
        circuit: Circuit,
        n_shots: Optional[int] = None,
        valid_check: bool = True,
        **kwargs: KwargTypes,
    ) -> ResultHandle:
        """
        Submit a single circuit to the backend for running. See
        :py:meth:`Backend.process_circuits`.
        """

        return self.process_circuits(
            [circuit], n_shots=n_shots, valid_check=valid_check, **kwargs
        )[0]

    @abstractmethod
    def process_circuits(
        self,
        circuits: Sequence[Circuit],
        n_shots: Optional[Union[int, Sequence[int]]] = None,
        valid_check: bool = True,
        **kwargs: KwargTypes,
    ) -> List[ResultHandle]:

        """
        Submit circuits to the backend for running. The results will be stored
        in the backend's result cache to be retrieved by the corresponding
        get_<data> method.

        If the `postprocess` keyword argument is set to True, and the backend supports
        the feature (see  :py:meth:`supports_contextual_optimisation`), then contextual
        optimisatioons are applied before running the circuit and retrieved results will
        have any necessary classical postprocessing applied. This is not enabled by
        default.

        Use keyword arguments to specify parameters to be used in submitting circuits
        See specific Backend derived class for available parameters, from the following
        list:

        * `seed`: RNG seed for simulators
        * `postprocess`: if True, apply contextual optimisations

        :param circuits: Circuits to process on the backend.
        :type circuits: Sequence[Circuit]
        :param n_shots: Number of shots to run per circuit. Optionally, this can be
            a list of shots specifying the number of shots for each circuit separately.
            None is to be used for state/unitary simulators. Defaults to None.
        :type n_shots: Optional[Union[int, Iterable[int]], optional
        :param valid_check: Explicitly check that all circuits satisfy all required
            predicates to run on the backend. Defaults to True
        :type valid_check: bool, optional
        :return: Handles to results for each input circuit, as an interable in
            the same order as the circuits.
        :rtype: List[ResultHandle]
        """
        ...

    @overload
    @staticmethod
    def _get_n_shots_as_list(
        n_shots: Union[None, int, Sequence[Optional[int]]],
        n_circuits: int,
        optional: Literal[False],
    ) -> List[int]:
        ...

    @overload
    @staticmethod
    def _get_n_shots_as_list(
        n_shots: Union[None, int, Sequence[Optional[int]]],
        n_circuits: int,
        optional: Literal[True],
        set_zero: Literal[True],
    ) -> List[int]:
        ...

    @overload
    @staticmethod
    def _get_n_shots_as_list(
        n_shots: Union[None, int, Sequence[Optional[int]]],
        n_circuits: int,
        optional: bool = True,
        set_zero: bool = False,
    ) -> Union[List[Optional[int]], List[int]]:
        ...

    @staticmethod
    def _get_n_shots_as_list(
        n_shots: Union[None, int, Sequence[Optional[int]]],
        n_circuits: int,
        optional: bool = True,
        set_zero: bool = False,
    ) -> Union[List[Optional[int]], List[int]]:
        """
        Convert any admissible n_shots value into List[Optional[int]] format.

        This validates the n_shots argument for process_circuits. If a single
        value is passed, this value is broadcast to the number of circuits.
        Additional boolean flags control how the argument is validated.
        Raises an exception if n_shots is in an invalid format.

        :param n_shots: The argument to be validated.
        :type n_shots: Union[None, int, Sequence[Optional[int]]]
        :param n_circuits: Length of the converted argument returned.
        :type n_circuits: int
        :param optional: Whether n_shots can be None (default: True).
        :type optional: bool
        :param set_zero: Whether None values should be set to 0 (default: False).
        :type set_zero: bool
        :return: a list of length `n_circuits`, the converted argument
        """

        n_shots_list: List[Optional[int]] = []

        def validate_n_shots(n: Optional[int]) -> bool:
            return optional or (n is not None and n > 0)

        if set_zero and not optional:
            ValueError("set_zero cannot be true when optional is false")

        if hasattr(n_shots, "__iter__"):
            assert not isinstance(n_shots, int)
            assert n_shots is not None

            if not all(map(validate_n_shots, n_shots)):
                raise ValueError(
                    "n_shots values are required for all circuits for this backend"
                )
            n_shots_list = list(n_shots)
        else:
            assert n_shots is None or isinstance(n_shots, int)

            if not validate_n_shots(n_shots):
                raise ValueError("Parameter n_shots is required for this backend")
            # convert n_shots to a list
            n_shots_list = [n_shots] * n_circuits

        if len(n_shots_list) != n_circuits:
            raise ValueError("The length of n_shots and circuits must match")

        if set_zero:
            # replace None with 0
            n_shots_list = list(map(lambda n: n or 0, n_shots_list))

        return n_shots_list

    @abstractmethod
    def circuit_status(self, handle: ResultHandle) -> CircuitStatus:
        """
        Return a CircuitStatus reporting the status of the circuit execution
        corresponding to the ResultHandle
        """
        ...

    def empty_cache(self) -> None:
        """Manually empty the result cache on the backend."""
        self._cache = {}

    def pop_result(self, handle: ResultHandle) -> Optional[ResultCache]:
        """Remove cache entry corresponding to handle from the cache and return.

        :param handle: ResultHandle object
        :type handle: ResultHandle
        :return: Cache entry corresponding to handle, if it was present
        :rtype: Optional[ResultCache]
        """
        return self._cache.pop(handle, None)

    @property
    def characterisation(self) -> Optional[dict]:
        """Retrieve the characterisation targeted by the backend if it exists.

        :return: The characterisation that this backend targets if it exists. The
            characterisation object contains device-specific information such as gate
            error rates.
        :rtype: Optional[dict]
        """
        raise NotImplementedError(
            "Backend does not support retrieving characterisation."
        )

    @property
    def backend_info(self) -> Optional[BackendInfo]:
        """Retrieve all Backend properties in a BackendInfo object, including
        device architecture, supported gate set, gate errors and other hardware-specific
        information.

        :return: The BackendInfo describing this backend if it exists.
        :rtype: Optional[BackendInfo]
        """
        raise NotImplementedError("Backend does not provide any device properties.")

    def get_result(self, handle: ResultHandle, **kwargs: KwargTypes) -> BackendResult:
        """Return a BackendResult corresponding to the handle.

        Use keyword arguments to specify parameters to be used in retrieving results.
        See specific Backend derived class for available parameters, from the following
        list:

        * `timeout`: maximum time to wait for remote job to finish
        * `wait`: polling interval between remote calls to check job status

        :param handle: handle to results
        :type handle: ResultHandle
        :return: Results corresponding to handle.
        :rtype: BackendResult
        """
        self._check_handle_type(handle)
        if handle in self._cache and "result" in self._cache[handle]:
            return cast(BackendResult, self._cache[handle]["result"])
        raise CircuitNotRunError(handle)

    def get_results(
        self, handles: Iterable[ResultHandle], **kwargs: KwargTypes
    ) -> List[BackendResult]:
        """Return results corresponding to handles.

        :param handles: Iterable of handles
        :return: List of results

        Keyword arguments are as for `get_result`, and apply to all jobs.
        """
        return [self.get_result(handle, **kwargs) for handle in handles]

    def _process(
        self, circuit: Circuit, **kwargs: KwargTypes
    ) -> Tuple[BackendResult, ResultHandle]:
        handle = self.process_circuit(circuit, **kwargs)  # type: ignore
        result = self.get_result(handle, **kwargs)
        self.pop_result(handle)
        return result, handle

    def cancel(self, handle: ResultHandle) -> None:
        """
        Cancel a job.

        :param handle: handle to job
        :type handle: ResultHandle
        :raises NotImplementedError: If backend does not support job cancellation
        """
        raise NotImplementedError("Backend does not support job cancellation.")

    @property
    def supports_shots(self) -> bool:
        """Does this backend support shot result retrieval via `get_shots`."""
        return self._supports_shots

    def get_shots(
        self,
        circuit: Circuit,
        n_shots: Optional[int] = None,
        basis: BasisOrder = BasisOrder.ilo,
        valid_check: bool = True,
        **kwargs: KwargTypes,
    ) -> np.ndarray:
        """Obtain the table of shots from an experiment. Accepts a
        :py:class:`~pytket.circuit.Circuit` to be run and immediately returned. This
        will fail if the circuit does not match the device's requirements.

        If the `postprocess` keyword argument is set to True, and the backend supports
        the feature (see  :py:meth:`supports_contextual_optimisation`), then contextual
        optimisations are applied before running the circuit and retrieved results will
        have any necessary classical postprocessing applied. This is not enabled by
        default.

        :param circuit: The circuit to run
        :type circuit: Circuit
        :param n_shots: Number of shots to generate from the circuit. Defaults to None
        :type n_shots: Optional[int], optional
        :param basis: Toggle between ILO (increasing lexicographic order of bit ids) and
            DLO (decreasing lexicographic order) for column ordering. Defaults to
            BasisOrder.ilo.
        :type basis: BasisOrder, optional
        :param valid_check: Explicitly check that the circuit satisfies all of the
            required predicates before running. Defaults to True
        :type valid_check: bool, optional
        :raises NotImplementedError: If backend implementation does not support shot
            table retrieval
        :return: Table of shot results. Each row is a single shot, with columns ordered
            by classical bit order (according to `basis`). Entries are 0 or 1
            corresponding to qubit basis states.
        :rtype: np.ndarray
        """
        result, _ = self._process(
            circuit, n_shots=n_shots, valid_check=valid_check, **kwargs
        )
        c_bits = (
            sorted(result.c_bits.keys(), reverse=(basis is not BasisOrder.ilo))
            if result.c_bits
            else None
        )
        return result.get_shots(c_bits)

    @property
    def supports_counts(self) -> bool:
        """Does this backend support counts result retrieval via `get_counts`."""
        return self._supports_counts

    def get_counts(
        self,
        circuit: Circuit,
        n_shots: Optional[int] = None,
        basis: BasisOrder = BasisOrder.ilo,
        valid_check: bool = True,
        **kwargs: KwargTypes,
    ) -> Dict[Tuple[int, ...], int]:
        """Obtain a summary of results, accumulating the shots for each result from an
        experiment. Accepts a
        :py:class:`~pytket.circuit.Circuit` to be run and immediately returned. This
        will fail if the circuit does not match the device's requirements.

        If the `postprocess` keyword argument is set to True, and the backend supports
        the feature (see  :py:meth:`supports_contextual_optimisation`), then contextual
        optimisatioons are applied before running the circuit and retrieved results will
        have any necessary classical postprocessing applied. This is not enabled by
        default.

        :param circuit: The circuit to run
        :type circuit: Circuit
        :param n_shots: Number of shots to generate from the circuit. Defaults to None
        :type n_shots: Optional[int], optional
        :param basis: Toggle between ILO (increasing lexicographic order of bit ids) and
            DLO (decreasing lexicographic order) for column ordering. Defaults to
            BasisOrder.ilo.
        :type basis: BasisOrder, optional
        :param valid_check: Explicitly check that the circuit satisfies all of the
            required predicates before running. Defaults to True
        :type valid_check: bool, optional
        :raises NotImplementedError: If backend implementation does not support counts
            retrieval
        :return: Dictionary mapping observed readouts to the number of times observed.
        :rtype: Dict[Tuple[int, ...], int]
        """

        result, _ = self._process(
            circuit, n_shots=n_shots, valid_check=valid_check, **kwargs
        )
        c_bits = (
            sorted(result.c_bits.keys(), reverse=(basis is not BasisOrder.ilo))
            if result.c_bits
            else None
        )
        try:
            return result.get_counts(c_bits)
        except InvalidResultType:
            shots = self.get_shots(
                circuit, n_shots=n_shots, basis=basis, valid_check=valid_check, **kwargs
            )
            return counts_from_shot_table(shots)

    @property
    def supports_state(self) -> bool:
        """Does this backend support statevector retrieval via `get_state`."""
        return self._supports_state

    def get_state(
        self,
        circuit: Circuit,
        basis: BasisOrder = BasisOrder.ilo,
        valid_check: bool = True,
    ) -> np.ndarray:
        """Obtain a statevector from a simulation. Accepts a
        :py:class:`~pytket.circuit.Circuit` to be run and immediately returned. This
        will fail if the circuit does not match the simulator's requirements.

        :param circuit: The circuit to run
        :type circuit: Circuit
        :param basis: Toggle between ILO-BE (increasing lexicographic order of bit ids,
            big-endian) and DLO-BE (decreasing lexicographic order, big-endian) for
            ordering the coefficients. Defaults to BasisOrder.ilo.
        :type basis: BasisOrder, optional
        :param valid_check: Explicitly check that the circuit satisfies all of the
            required predicates before running. Defaults to True
        :type valid_check: bool, optional
        :raises NotImplementedError: If backend implementation does not support
            statevector retrieval
        :return: A big-endian statevector for the circuit in the encoding given by
            `basis`; e.g. :math:`[a_{00}, a_{01}, a_{10}, a_{11}]` where :math:`a_{01}`
            is the amplitude of the :math:`\\left|01\\right>` state (in ILO, this means
            qubit q[0] is in state :math:`\\left|0\\right>` and q[1] is in state
            :math:`\\left|1\\right>`, and the reverse in DLO)
        :rtype: np.ndarray
        """

        result, _ = self._process(circuit, valid_check=valid_check)
        q_bits = (
            sorted(result.q_bits.keys(), reverse=(basis is not BasisOrder.ilo))
            if result.q_bits
            else None
        )
        return result.get_state(q_bits)

    @property
    def supports_unitary(self) -> bool:
        """Does this backend support unitary retrieval via `get_unitary`."""
        return self._supports_unitary

    @property
    def supports_density_matrix(self) -> bool:
        """Does this backend support density matrix retrieval via
        `get_density_matrix`."""
        return self._supports_density_matrix

    @property
    def supports_expectation(self) -> bool:
        """Does this backend support expectation value calculation for operators."""
        return self._supports_expectation

    @property
    def expectation_allows_nonhermitian(self) -> bool:
        """If expectations are supported, is the operator allowed to be non-Hermitan?"""
        return self._expectation_allows_nonhermitian

    @property
    def supports_contextual_optimisation(self) -> bool:
        """Does this backend support contextual optimisation?

        See :py:meth:`process_circuits`."""
        return self._supports_contextual_optimisation

    def _get_extension_module(self) -> Optional[ModuleType]:
        """Return the extension module of the backend if it belongs to a
        pytket-extension package.

        :return: The extension module of the backend if it belongs to a pytket-extension
            package.
        :rtype: Optional[ModuleType]
        """
        mod_parts = self.__class__.__module__.split(".")[:3]
        if not (mod_parts[0] == "pytket" and mod_parts[1] == "extensions"):
            return None
        return import_module(".".join(mod_parts))

    @property
    def __extension_name__(self) -> Optional[str]:
        """Retrieve the extension name of the backend if it belongs to a
        pytket-extension package.

        :return: The extension name of the backend if it belongs to a pytket-extension
            package.
        :rtype: Optional[str]
        """
        try:
            return self._get_extension_module().__extension_name__  # type: ignore
        except AttributeError:
            return None

    @property
    def __extension_version__(self) -> Optional[str]:
        """Retrieve the extension version of the backend if it belongs to a
        pytket-extension package.

        :return: The extension version of the backend if it belongs to a
            pytket-extension package.
        :rtype: Optional[str]
        """
        try:
            return self._get_extension_module().__extension_version__  # type: ignore
        except AttributeError:
            return None
