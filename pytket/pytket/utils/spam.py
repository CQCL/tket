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

import itertools
from functools import lru_cache

from math import ceil, log2
from collections import OrderedDict
from typing import Dict, Iterable, List, Tuple, Counter, cast, Optional, Callable, Union
import numpy as np
from pytket.circuit import Circuit, Qubit, Bit, Node, CircBox, OpType
from pytket.backends import Backend
from pytket.passes import DecomposeBoxes, FlattenRegisters
from pytket.backends.backendresult import BackendResult
from pytket.utils.outcomearray import OutcomeArray
from pytket.utils.results import CountsDict, StateTuple


ParallelMeasures = List[Dict[Qubit, Bit]]


def compress_counts(
    counts: Dict[StateTuple, float], tol: float = 1e-6, round_to_int: bool = False
) -> CountsDict:
    """Filter counts to remove states that have a count value (which can be a
    floating-point number) below a tolerance, and optionally round to an
    integer.

    :param counts: Input counts
    :type counts: Dict[StateTuple, float]
    :param tol: Value below which counts are pruned. Defaults to 1e-6.
    :type tol: float, optional
    :param round_to_int: Whether to round each count to an integer. Defaults to False.
    :type round_to_int: bool, optional

    :return: Filtered counts
    :rtype: CountsDict
    """
    valprocess: Callable[[float], Union[int, float]] = (
        lambda x: int(round(x)) if round_to_int else x
    )
    processed_pairs = (
        (key, valprocess(val)) for key, val in counts.items() if val > tol
    )
    return {key: val for key, val in processed_pairs if val > 0}


@lru_cache(maxsize=128)
def binary_to_int(bintuple: Tuple[int]) -> int:
    """Convert a binary tuple to corresponding integer, with most significant bit as
    the first element of tuple.

    :param bintuple: Binary tuple
    :type bintuple: Tuple[int]

    :return:
    Integer :rtype: int
    """
    integer = 0
    for index, bitset in enumerate(reversed(bintuple)):
        if bitset:
            integer |= 1 << index
    return integer


@lru_cache(maxsize=128)
def int_to_binary(val: int, dim: int) -> Tuple[int, ...]:
    """Convert an integer to corresponding binary tuple, with most significant bit as
     the first element of tuple.

    :param val: input integer
    :type val: int
    :param dim: Bit width
    :type dim: int

    :return: Binary tuple of width dim
    :rtype: Tuple[int, ...]
    """
    return tuple(map(int, format(val, "0{}b".format(dim))))


#########################################
### _compute_dot and helper functions ###
###
### With thanks to
### https://math.stackexchange.com/a/3423910
### and especially
### https://gist.github.com/ahwillia/f65bc70cb30206d4eadec857b98c4065
### on which this code is based.
def _unfold(tens: np.ndarray, mode: int, dims: List[int]) -> np.ndarray:
    """Unfolds tensor into matrix.

    :param tens: Tensor with shape equivalent to dimensions
    :type tens: np.ndarray
    :param mode: Specifies axis move to front of matrix in unfolding of tensor
    :type mode: int
    :param dims: Gives shape of tensor passed
    :type dims: List[int]

    :return: Matrix with shape (dims[mode], prod(dims[/mode]))
    :rtype: np.ndarray
    """
    if mode == 0:
        return tens.reshape(dims[0], -1)
    else:
        return np.moveaxis(tens, mode, 0).reshape(dims[mode], -1)


def _refold(vec: np.ndarray, mode: int, dims: List[int]) -> np.ndarray:
    """Refolds vector into tensor.

    :param vec: Tensor with length equivalent to the product of dimensions given in
        dims
    :type vec: np.ndarray
    :param mode: Axis tensor was unfolded along
    :type mode: int
    :param dims: Shape of tensor
    :type dims: List[int]

    :return: Tensor folded from vector with shape equivalent to given dimensions
    :rtype: np.ndarray
    """
    if mode == 0:
        return vec.reshape(dims)
    else:
        # Reshape and then move dims[mode] back to its
        # appropriate spot (undoing the `unfold` operation).
        tens = vec.reshape([dims[mode]] + [d for m, d in enumerate(dims) if m != mode])
        return np.moveaxis(tens, 0, mode)


def _compute_dot(submatrices: Iterable[np.ndarray], vector: np.ndarray) -> np.ndarray:
    """Multiplies the kronecker product of the given submatrices with given vector.

    :param submatrices: Submatrices multiplied
    :type submatrices: Iterable[np.ndarray]
    :param vector: Vector multplied
    :type vector: np.ndarray

    :return: Kronecker product of arguments
    :rtype: np.ndarray
    """
    dims = [A.shape[0] for A in submatrices]
    vt = vector.reshape(dims)
    for i, A in enumerate(submatrices):
        vt = _refold(A @ _unfold(vt, i, dims), i, dims)
    return vt.ravel()


def _bayesian_iteration(
    submatrices: Iterable[np.ndarray],
    measurements: np.ndarray,
    t: np.ndarray,
    epsilon: float,
) -> np.ndarray:
    """Transforms T corresponds to a Bayesian iteration, used to modfiy
    measurements.

    :param submatrices: submatrices to be inverted and applied to measurements.
    :type submatrices: Iterable[np.ndarray]
    :param measurements: Probability distribution over set of states to be amended.
    :type measurements: np.ndarray
    :param t: Some transform to act on measurements.
    :type t: np.ndarray
    :param epsilon: A stabilization parameter to define an affine transformation for
        application to submatrices, eliminating zero probabilities.
    :type epsilon: float

    :return: Transformed distribution vector.
    :rtype: np.ndarray
    """
    # Transform t according to the Bayesian iteration
    # The parameter epsilon is a stabilization parameter which defines an affine
    # transformation to apply to the submatrices to eliminate zero probabilities. This
    # transformation preserves the property that all columns sum to 1
    if epsilon == 0:
        # avoid copying if we don't need to
        As = submatrices
    else:
        As = [
            epsilon / submatrix.shape[0] + (1 - epsilon) * submatrix
            for submatrix in submatrices
        ]
    z = _compute_dot(As, t)
    if np.isclose(z, 0).any():
        raise ZeroDivisionError
    return cast(
        np.ndarray, t * _compute_dot([A.transpose() for A in As], measurements / z)
    )


def _bayesian_iterative_correct(
    submatrices: Iterable[np.ndarray],
    measurements: np.ndarray,
    tol: float = 1e-5,
    max_it: Optional[int] = None,
) -> np.ndarray:
    """Finds new states to represent application of inversion of submatrices on
    measurements. Converges when update states within tol range of previously
    tested states.

    :param submatrices: Matrices comprising the pure noise characterisation.
    :type submatrices: Iterable[np.ndarray]
    :param input_vector: Vector corresponding to some counts distribution.
    :type input_vector: np.ndarray
    :param tol: tolerance of closeness of found results
    :type tol: float
    :param max_it: Maximum number of inversions attempted to correct results.
    :type max_it: int
    """
    # based on method found in https://arxiv.org/abs/1910.00129

    vector_size = measurements.size
    # uniform initial
    true_states = np.full(vector_size, 1 / vector_size)
    prev_true = true_states.copy()
    converged = False
    count = 0
    epsilon: float = 0  # stabilization parameter, adjusted dynamically
    while not converged:
        if max_it:
            if count >= max_it:
                break
            count += 1
        try:
            true_states = _bayesian_iteration(
                submatrices, measurements, true_states, epsilon
            )
            converged = np.allclose(true_states, prev_true, atol=tol)
            prev_true = true_states.copy()
        except ZeroDivisionError:
            # Shift the stabilization parameter up a bit (always < 0.5).
            epsilon = 0.99 * epsilon + 0.01 * 0.5

    return true_states


def reduce_matrix(indices_to_remove: List[int], matrix: np.ndarray) -> np.ndarray:
    """Removes indices from indices_to_remove from binary associated to indexing of
    matrix, producing a new transition matrix. To do so, it assigns all transition
    probabilities as the given state in the remaining indices binary, with the removed
    binary in state 0. This is an assumption on the noise made because it is likely
    that unmeasured qubits will be in that state.

    :param indices_to_remove: Binary index of state matrix is mapping to be removed.
    :type indices_to_remove: List[int]
    :param matrix: Transition matrix where indices correspond to some binary state.
    :type matrix: np.ndarray

    :return: Transition matrix with removed entries.
    :rtype: np.ndarray
    """

    new_n_qubits = int(log2(matrix.shape[0])) - len(indices_to_remove)
    if new_n_qubits == 0:
        return np.array([])
    bin_map = dict()
    mat_dim = 1 << new_n_qubits
    for index in range(mat_dim):
        # get current binary
        bina = list(int_to_binary(index, new_n_qubits))
        # add 0's to fetch old binary to set values from
        for i in sorted(indices_to_remove):
            bina.insert(i, 0)
        # get index of values
        bin_map[index] = binary_to_int(tuple(bina))

    new_mat = np.zeros((mat_dim,) * 2, dtype=float)
    for i in range(len(new_mat)):
        old_row_index = bin_map[i]
        for j in range(len(new_mat)):
            old_col_index = bin_map[j]
            new_mat[i, j] = matrix[old_row_index, old_col_index]
    return new_mat


def reduce_matrices(
    entries_to_remove: List[Tuple[int, int]], matrices: List[np.ndarray]
) -> List[np.ndarray]:
    """Removes some dimensions from some matrices.

    :param entries_to_remove: Via indexing, details dimensions to be removed.
    :type entries_to_remove: List[Tuple[int, int]]
    :param matrices: All matrices to have dimensions removed.
    :type matrices: List[np.ndarray]

    :return: Matrices with some dimensions removed.
    :rtype: List[np.ndarray]
    """
    organise: Dict[int, List] = dict({k: [] for k in range(len(matrices))})
    for unused in entries_to_remove:
        # unused[0] is index in matrices
        # unused[1] is qubit index in matrix
        organise[unused[0]].append(unused[1])
    output_matrices = [reduce_matrix(organise[m], matrices[m]) for m in organise]
    normalised_mats = [
        mat / np.sum(mat, axis=0) for mat in [x for x in output_matrices if len(x) != 0]
    ]
    return normalised_mats


class SpamCorrecter:
    """A class for generating "state preparation and measurement" (SPAM) calibration
    experiments for ``pytket`` backends, and correcting counts generated from them.

    Supports saving calibrated state to a dictionary format, and restoring from the
    dictionary.
    """

    def __init__(
        self, qubit_subsets: List[List[Node]], backend: Optional[Backend] = None
    ):
        """Construct a new `SpamCorrecter`.

        :param qubit_subsets: A list of lists of correlated Nodes of an `Architecture`.
            Qubits within the same list are assumed to only have SPAM errors correlated
            with each other. Thus to allow SPAM errors between all qubits you should
            provide a single list.
        :type qubit_subsets: List[List[Node]]
        :param backend: Backend on which the experiments are intended to be run
            (optional). If provided, the qubits in `qubit_subsets` must be nodes in the
            backend's associated `Architecture`. If not provided, it is assumed that the
            experiment will be run on an `Architecture`with the nodes in
            `qubit_subsets`, and furthermore that the intended architecture natively
            supports X gates.

        :raises ValueError: There are repeats in the `qubit_subsets` specification.
        """
        self.correlations = qubit_subsets

        self.all_qbs = [qb for subset in qubit_subsets for qb in subset]

        def to_tuple(inp: list[Node]) -> tuple:
            return tuple(inp)

        self.subsets_matrix_map = OrderedDict.fromkeys(
            sorted(map(to_tuple, self.correlations), key=len, reverse=True)
        )
        # ordered from largest to smallest via OrderedDict & sorted
        self.subset_dimensions = [len(subset) for subset in self.subsets_matrix_map]

        if len(self.all_qbs) != len(set(self.all_qbs)):
            raise ValueError("Qubit subsets are not mutually disjoint.")

        xcirc = Circuit(1).X(0)
        if backend is not None:
            if backend.backend_info is None:
                raise ValueError("No architecture associated with backend.")
            nodes = backend.backend_info.nodes
            if not all(node in nodes for node in self.all_qbs):
                raise ValueError("Nodes do not all belong to architecture.")
            backend.default_compilation_pass().apply(xcirc)
            FlattenRegisters().apply(xcirc)

        self.xbox = CircBox(xcirc)

    def calibration_circuits(self) -> List[Circuit]:
        """Generate calibration circuits according to the specified correlations.

        :return: A list of calibration circuits to be run on the machine. The circuits
            should be processed without compilation. Results from these circuits must
            be given back to this class (via the `calculate_matrices` method) in the
            same order.
        :rtype: List[Circuit]
        """

        major_state_dimensions = self.subset_dimensions[0]
        n_circuits = 1 << major_state_dimensions
        # output
        self.prepared_circuits = []
        self.state_infos = []

        # set up base circuit for appending xbox to
        base_circuit = Circuit()
        c_reg = []
        for index, qb in enumerate(self.all_qbs):
            base_circuit.add_qubit(qb)
            c_bit = Bit(index)
            c_reg.append(c_bit)
            base_circuit.add_bit(c_bit)

        # generate state circuits for given correlations
        for major_state_index in range(n_circuits):
            state_circuit = base_circuit.copy()
            # get bit string corresponding to basis state of biggest subset of qubits
            major_state = int_to_binary(major_state_index, major_state_dimensions)
            new_state_dicts = {}
            # parallelise circuits, run uncorrelated subsets
            # characterisation in parallel
            for dim, qubits in zip(self.subset_dimensions, self.subsets_matrix_map):
                # add state to prepared states
                new_state_dicts[qubits] = major_state[:dim]
                # find only qubits that are expected to be in 1 state,
                # add xbox to given qubits
                for flipped_qb in itertools.compress(qubits, major_state[:dim]):
                    state_circuit.add_circbox(self.xbox, [flipped_qb])
            # Decompose boxes, add barriers to preserve circuit, add measures
            DecomposeBoxes().apply(state_circuit)
            for qb, cb in zip(self.all_qbs, c_reg):
                state_circuit.Measure(qb, cb)

            # add to returned types
            self.prepared_circuits.append(state_circuit)
            self.state_infos.append((new_state_dicts, state_circuit.qubit_to_bit_map))
        return self.prepared_circuits

    def calculate_matrices(self, results_list: List[BackendResult]) -> None:
        """Calculate the calibration matrices from the results of running calibration
        circuits.

        :param results_list: List of results from Backend. Must be in the same order as
             the corresponding circuits generated by `calibration_circuits`.
        :type counts_list: List[BackendResult]

        :raises RuntimeError: Calibration circuits have not been generated yet.
        """
        if not self.state_infos:
            raise RuntimeError(
                "Ensure calibration states/circuits have been calculated first."
            )

        counter = 0
        self.node_index_dict: Dict[Node, Tuple[int, int]] = dict()

        for qbs, dim in zip(self.subsets_matrix_map, self.subset_dimensions):
            # for a subset with n qubits, create a 2^n by 2^n matrix
            self.subsets_matrix_map[qbs] = np.zeros((1 << dim,) * 2, dtype=float)
            for i in range(len(qbs)):
                qb = qbs[i]
                self.node_index_dict[qb] = (counter, i)
            counter += 1

        for result, state_info in zip(results_list, self.state_infos):
            state_dict = state_info[0]
            qb_bit_map = state_info[1]
            for qb_sub in self.subsets_matrix_map:
                # bits of counts to consider
                bits = [qb_bit_map[q] for q in qb_sub]
                counts_dict = result.get_counts(cbits=bits)
                for measured_state, count in counts_dict.items():
                    # intended state
                    prepared_state_index = binary_to_int(state_dict[qb_sub])
                    # produced state
                    measured_state_index = binary_to_int(measured_state)
                    # update characterisation matrix
                    M = self.subsets_matrix_map[qb_sub]
                    assert type(M) is np.ndarray
                    M[measured_state_index, prepared_state_index] += count

        # normalise everything
        self.characterisation_matrices = [
            mat / np.sum(cast(np.ndarray, mat), axis=0)
            for mat in self.subsets_matrix_map.values()
        ]

    def get_parallel_measure(self, circuit: Circuit) -> ParallelMeasures:
        """For a given circuit, produces and returns a ParallelMeasures object required
         for correcting counts results.

        :param circuit: Circuit with some Measure operations.
        :type circuit: Circuit

        :return: A list of dictionaries mapping Qubit to Bit where each separate
            dictionary details some set of Measurement operations run in parallel.
        :rtype: ParallelMeasures
        """
        parallel_measure = [circuit.qubit_to_bit_map]
        # implies mid-circuit measurements, or that at least missing
        # bits need to be checked for Measure operation
        if len(parallel_measure[0]) != len(circuit.bits):
            used_bits = set(parallel_measure[0].values())
            for mc in circuit.commands_of_type(OpType.Measure):
                bit = mc.bits[0]
                if bit not in used_bits:
                    # mid-circuit measure, add as a separate parallel measure
                    parallel_measure.append({mc.qubits[0]: bit})
        return parallel_measure

    def correct_counts(
        self,
        result: BackendResult,
        parallel_measures: ParallelMeasures,
        method: str = "bayesian",
        options: Optional[Dict] = None,
    ) -> BackendResult:
        """Modifies count distribution for result, such that the inversion of the pure
        noise map represented by characterisation matrices is applied to it.

        :param result: BackendResult object to be negated by pure noise object.
        :type result: BackendResult
        :param parallel_measures: Used to permute corresponding BackendResult object so
             counts order matches noise characterisation and to amend characterisation
             matrices to correct the right bits. SpamCorrecter.get_parallel_measure
             returns the required object for a given circuit.
        :type parallel_measures: ParallelMeasures

        :raises ValueError: Measured qubit in result not characterised.

        :return: A new result object with counts modified to reflect SPAM correction.
        :rtype: BackendResult
        """
        # the correction process assumes that when passed a list of matrices
        #  and a distribution to correct, that the j rows of matrix i
        # corrects for the i, i+1,...i+j states in the passed distribution
        # given information of which bits are measured on which qubits from
        # parallel_measures, the following first produces matrices such that
        # this condition is true

        char_bits_order = []
        correction_matrices = []

        for mapping in parallel_measures:
            # reduce_matrices removes given qubits corresponding entries from
            # characterisation matrices
            unused_qbs = set(self.all_qbs.copy())
            for q in mapping:
                # no q duplicates as mapping is dict from qubit to bit
                if q not in unused_qbs:
                    raise ValueError(
                        "Measured qubit {} is not characterised by "
                        "SpamCorrecter".format(q)
                    )
                unused_qbs.remove(q)  # type:ignore[arg-type]
                char_bits_order.append(mapping[q])
            correction_matrices.extend(
                reduce_matrices(
                    [self.node_index_dict[q] for q in unused_qbs],
                    self.characterisation_matrices,
                )
            )

        # get counts object for returning later
        counts = result.get_counts(cbits=char_bits_order)
        in_vec = np.zeros(1 << len(char_bits_order), dtype=float)
        # turn from counts to probability distribution
        for state, count in counts.items():
            in_vec[binary_to_int(state)] = count
        Ncounts = np.sum(in_vec)
        in_vec_norm = in_vec / Ncounts

        # with counts and characterisation matrices orders matching,
        # correct distribution
        if method == "invert":
            try:
                subinverts = [
                    np.linalg.inv(submatrix) for submatrix in correction_matrices
                ]
            except np.linalg.LinAlgError:
                raise ValueError(
                    "Unable to invert calibration matrix: please re-run "
                    "calibration experiments or use an alternative correction method."
                )
            # assumes that order of rows in flattened subinverts equals
            # order of bits in input vector
            outvec = _compute_dot(subinverts, in_vec_norm)
            # The entries of v will always sum to 1, but they may not all
            #  be in the range [0,1]. In order to make them genuine
            # probabilities (and thus generate meaningful counts),
            # we adjust them by setting all negative values to 0 and scaling
            #  the remainder.
            outvec[outvec < 0] = 0
            outvec /= sum(outvec)

        elif method == "bayesian":
            if options is None:
                options = {}
            tol_val = options.get("tol", 1 / Ncounts)
            maxit = options.get("maxiter", None)
            outvec = _bayesian_iterative_correct(
                correction_matrices, in_vec_norm, tol=tol_val, max_it=maxit
            )

        else:
            valid_methods = ("invert", "bayesian")
            raise ValueError("Method must be one of: ", *valid_methods)

        outvec *= Ncounts

        # counter object with binary from distribution
        corrected_counts = {
            int_to_binary(index, len(char_bits_order)): Bcount
            for index, Bcount in enumerate(outvec)
        }
        counter = Counter(
            {
                OutcomeArray.from_readouts([key]): ceil(val)
                for key, val in corrected_counts.items()
            }
        )
        # produce and return BackendResult object
        return BackendResult(counts=counter, c_bits=char_bits_order)

    def to_dict(self) -> Dict:
        """Get calibration information as a dictionary.

        :return: Dictionary output
        :rtype: Dict
        """
        correlations = []
        for subset in self.correlations:
            correlations.append([(uid.reg_name, uid.index) for uid in subset])

        node_index_hashable = [
            ((uid.reg_name, uid.index), self.node_index_dict[uid])
            for uid in self.node_index_dict
        ]
        char_matrices = [m.tolist() for m in self.characterisation_matrices]
        self_dict = {
            "correlations": correlations,
            "node_index_dict": node_index_hashable,
            "characterisation_matrices": char_matrices,
        }
        return self_dict

    @classmethod
    def from_dict(class_obj, d: Dict) -> "SpamCorrecter":
        """Build a `SpamCorrecter` instance from a dictionary in the format returned
        by `to_dict`.

        :return: Dictionary of calibration information.
        :rtype: SpamCorrecter
        """
        new_inst = class_obj(
            [
                [Node(*pair)]
                for subset_tuple in d["correlations"]
                for pair in subset_tuple
            ]
        )
        new_inst.node_index_dict = dict(
            [
                (Node(*pair[0]), (int(pair[1][0]), int(pair[1][1])))
                for pair in d["node_index_dict"]
            ]
        )
        new_inst.characterisation_matrices = [
            np.array(m) for m in d["characterisation_matrices"]
        ]
        return new_inst
