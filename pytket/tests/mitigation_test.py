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

import tempfile
import json

from pytket.utils.spam import SpamCorrecter, compress_counts
from pytket.circuit import Node, Circuit, Qubit  # type: ignore
from pytket.mapping import MappingManager, LexiLabellingMethod, LexiRouteRoutingMethod  # type: ignore
from pytket.architecture import Architecture  # type: ignore
from pytket.placement import place_with_map  # type: ignore
from pytket.passes import DelayMeasures  # type: ignore
from typing import List, Dict, Counter, Tuple
from pytket.utils.outcomearray import OutcomeArray
from math import ceil
from pytket.backends.backendresult import BackendResult


def test_spam_integration() -> None:
    SEED = 120
    # test data, to avoid using backend
    calib_results: List[Dict[Tuple[int, ...], int]] = [
        {
            (0, 0, 0): 742,
            (0, 0, 1): 84,
            (0, 1, 0): 76,
            (0, 1, 1): 8,
            (1, 0, 0): 72,
            (1, 0, 1): 10,
            (1, 1, 0): 7,
            (1, 1, 1): 1,
        },
        {
            (0, 0, 0): 183,
            (0, 0, 1): 23,
            (0, 1, 0): 614,
            (0, 1, 1): 70,
            (1, 0, 0): 13,
            (1, 0, 1): 5,
            (1, 1, 0): 85,
            (1, 1, 1): 7,
        },
        {
            (0, 0, 0): 53,
            (0, 0, 1): 161,
            (0, 1, 0): 6,
            (0, 1, 1): 20,
            (1, 0, 0): 166,
            (1, 0, 1): 517,
            (1, 1, 0): 18,
            (1, 1, 1): 59,
        },
        {
            (0, 0, 0): 13,
            (0, 0, 1): 47,
            (0, 1, 0): 48,
            (0, 1, 1): 131,
            (1, 0, 0): 49,
            (1, 0, 1): 152,
            (1, 1, 0): 136,
            (1, 1, 1): 424,
        },
    ]
    bellres_counts: Dict[Tuple[int, ...], int] = {
        (0, 0, 0): 406,
        (0, 0, 1): 111,
        (0, 1, 0): 38,
        (0, 1, 1): 13,
        (1, 0, 0): 136,
        (1, 0, 1): 251,
        (1, 1, 0): 14,
        (1, 1, 1): 31,
    }
    bellres_counter = Counter(
        {
            OutcomeArray.from_readouts([key]): ceil(val)
            for key, val in bellres_counts.items()
        }
    )
    bellres = BackendResult(counts=bellres_counter)
    qbs = [Node("qx", i) for i in range(4)]
    arc = Architecture([[qbs[i], qbs[i + 1]] for i in range(3)])

    subs = [[qbs[2], qbs[0]], [qbs[1]]]
    spam = SpamCorrecter(subs)
    calib_circs = spam.calibration_circuits()

    assert len(calib_circs) == 4

    calib_counters = [
        Counter({OutcomeArray.from_readouts([key]): val for key, val in r.items()})
        for r in calib_results
    ]
    calib_brs = [BackendResult(counts=c) for c in calib_counters]
    spam.calculate_matrices(calib_brs)
    assert spam.characterisation_matrices[0].shape == (4, 4)
    assert spam.characterisation_matrices[1].shape == (2, 2)

    bellcc = Circuit(3, 3).H(0).CX(0, 2).measure_all()
    qmap = {Qubit(0): qbs[1], Qubit(1): qbs[2], Qubit(2): qbs[0]}
    place_with_map(bellcc, qmap)
    mm = MappingManager(arc)
    rbell = bellcc.copy()
    mm.route_circuit(rbell, [LexiLabellingMethod(), LexiRouteRoutingMethod()])

    def check_correction(
        counts0: Dict[Tuple[int, ...], int], counts1: Dict[Tuple[int, ...], int]
    ) -> bool:
        if (
            counts0[(0, 0, 0)] > counts1[(0, 0, 0)]
            and counts0[(1, 0, 1)] > counts1[(1, 0, 1)]
        ):
            return True
        return False

    rbell_parallel_measures = spam.get_parallel_measure(rbell)
    default_correct = spam.correct_counts(bellres, rbell_parallel_measures).get_counts()
    assert check_correction(default_correct, bellres_counts)

    bellcc_mid = Circuit()
    qbs = bellcc_mid.add_q_register("qx", 3)
    bits = bellcc_mid.add_c_register("c", 3)
    bellcc_mid.H(qbs[0]).CX(qbs[0], qbs[2]).Measure(qbs[0], bits[0]).SWAP(
        qbs[0], qbs[1]
    ).Measure(qbs[0], bits[1]).Measure(qbs[2], bits[2])

    bell_cc_parallel_measures = spam.get_parallel_measure(bellcc_mid)
    default_correct_mid = spam.correct_counts(
        bellres, bell_cc_parallel_measures
    ).get_counts()
    assert check_correction(default_correct_mid, bellres_counts)

    with tempfile.TemporaryFile(mode="w+t") as tmpjson:
        json.dump(spam.to_dict(), tmpjson)
        tmpjson.seek(0)
        newspam = SpamCorrecter.from_dict(json.load(tmpjson))

    new_default = newspam.correct_counts(bellres, rbell_parallel_measures).get_counts()
    assert default_correct == new_default

    assert check_correction(
        spam.correct_counts(
            bellres, rbell_parallel_measures, method="invert"
        ).get_counts(),
        bellres_counts,
    )
    assert check_correction(
        spam.correct_counts(
            bellres, rbell_parallel_measures, method="bayesian"
        ).get_counts(),
        bellres_counts,
    )
    assert check_correction(
        spam.correct_counts(
            bellres,
            rbell_parallel_measures,
            method="bayesian",
            options={"tol": 1e-8, "maxiter": 10},
        ).get_counts(),
        bellres_counts,
    )
    assert check_correction(
        spam.correct_counts(
            bellres, rbell_parallel_measures, method="bayesian", options={"tol": 1e-2}
        ).get_counts(),
        bellres_counts,
    )

    assert check_correction(
        spam.correct_counts(
            bellres, bell_cc_parallel_measures, method="invert"
        ).get_counts(),
        bellres_counts,
    )
    assert check_correction(
        spam.correct_counts(
            bellres, bell_cc_parallel_measures, method="bayesian"
        ).get_counts(),
        bellres_counts,
    )
    assert check_correction(
        spam.correct_counts(
            bellres,
            bell_cc_parallel_measures,
            method="bayesian",
            options={"tol": 1e-8, "maxiter": 10},
        ).get_counts(),
        bellres_counts,
    )
    assert check_correction(
        spam.correct_counts(
            bellres, bell_cc_parallel_measures, method="bayesian", options={"tol": 1e-2}
        ).get_counts(),
        bellres_counts,
    )


def test_spam_routing() -> None:
    # test spam with a pre routed circuit, using final mapped qubits to perform
    # the calibration
    raw_res = {
        (0, 0, 0, 0): 352,
        (0, 0, 0, 1): 38,
        (0, 0, 1, 0): 37,
        (0, 0, 1, 1): 15,
        (0, 1, 0, 0): 36,
        (0, 1, 0, 1): 17,
        (0, 1, 1, 0): 19,
        (0, 1, 1, 1): 58,
        (1, 0, 0, 0): 43,
        (1, 0, 0, 1): 25,
        (1, 0, 1, 0): 25,
        (1, 0, 1, 1): 41,
        (1, 1, 0, 0): 22,
        (1, 1, 0, 1): 60,
        (1, 1, 1, 0): 45,
        (1, 1, 1, 1): 167,
    }

    calib_results: List[Dict[Tuple[int, ...], int]] = [
        {
            (0, 0, 0, 0): 659,
            (0, 0, 0, 1): 65,
            (0, 0, 1, 0): 68,
            (0, 0, 1, 1): 10,
            (0, 1, 0, 0): 70,
            (0, 1, 0, 1): 12,
            (0, 1, 1, 0): 6,
            (1, 0, 0, 0): 86,
            (1, 0, 0, 1): 7,
            (1, 0, 1, 0): 7,
            (1, 0, 1, 1): 2,
            (1, 1, 0, 0): 6,
            (1, 1, 1, 0): 1,
            (1, 1, 1, 1): 1,
        },
        {
            (0, 0, 0, 0): 170,
            (0, 0, 0, 1): 555,
            (0, 0, 1, 0): 19,
            (0, 0, 1, 1): 71,
            (0, 1, 0, 0): 18,
            (0, 1, 0, 1): 50,
            (0, 1, 1, 0): 2,
            (0, 1, 1, 1): 8,
            (1, 0, 0, 0): 20,
            (1, 0, 0, 1): 72,
            (1, 0, 1, 1): 6,
            (1, 1, 0, 0): 2,
            (1, 1, 0, 1): 6,
            (1, 1, 1, 1): 1,
        },
        {
            (0, 0, 0, 0): 173,
            (0, 0, 0, 1): 17,
            (0, 0, 1, 0): 538,
            (0, 0, 1, 1): 75,
            (0, 1, 0, 0): 21,
            (0, 1, 0, 1): 3,
            (0, 1, 1, 0): 69,
            (0, 1, 1, 1): 12,
            (1, 0, 0, 0): 17,
            (1, 0, 0, 1): 4,
            (1, 0, 1, 0): 59,
            (1, 0, 1, 1): 7,
            (1, 1, 0, 0): 1,
            (1, 1, 1, 0): 4,
        },
        {
            (0, 0, 0, 0): 42,
            (0, 0, 0, 1): 150,
            (0, 0, 1, 0): 145,
            (0, 0, 1, 1): 444,
            (0, 1, 0, 0): 10,
            (0, 1, 0, 1): 19,
            (0, 1, 1, 0): 16,
            (0, 1, 1, 1): 73,
            (1, 0, 0, 0): 4,
            (1, 0, 0, 1): 17,
            (1, 0, 1, 0): 12,
            (1, 0, 1, 1): 63,
            (1, 1, 0, 1): 1,
            (1, 1, 1, 1): 4,
        },
        {
            (0, 0, 0, 0): 168,
            (0, 0, 0, 1): 26,
            (0, 0, 1, 0): 22,
            (0, 0, 1, 1): 6,
            (0, 1, 0, 0): 555,
            (0, 1, 0, 1): 65,
            (0, 1, 1, 0): 59,
            (0, 1, 1, 1): 4,
            (1, 0, 0, 0): 24,
            (1, 0, 0, 1): 2,
            (1, 0, 1, 0): 2,
            (1, 0, 1, 1): 1,
            (1, 1, 0, 0): 50,
            (1, 1, 0, 1): 11,
            (1, 1, 1, 0): 5,
        },
        {
            (0, 0, 0, 0): 60,
            (0, 0, 0, 1): 164,
            (0, 0, 1, 0): 7,
            (0, 0, 1, 1): 15,
            (0, 1, 0, 0): 158,
            (0, 1, 0, 1): 442,
            (0, 1, 1, 0): 14,
            (0, 1, 1, 1): 39,
            (1, 0, 0, 0): 7,
            (1, 0, 0, 1): 19,
            (1, 0, 1, 1): 2,
            (1, 1, 0, 0): 18,
            (1, 1, 0, 1): 49,
            (1, 1, 1, 0): 1,
            (1, 1, 1, 1): 5,
        },
        {
            (0, 0, 0, 0): 53,
            (0, 0, 0, 1): 3,
            (0, 0, 1, 0): 160,
            (0, 0, 1, 1): 17,
            (0, 1, 0, 0): 155,
            (0, 1, 0, 1): 16,
            (0, 1, 1, 0): 449,
            (0, 1, 1, 1): 55,
            (1, 0, 0, 0): 3,
            (1, 0, 1, 0): 11,
            (1, 0, 1, 1): 4,
            (1, 1, 0, 0): 17,
            (1, 1, 0, 1): 2,
            (1, 1, 1, 0): 50,
            (1, 1, 1, 1): 5,
        },
        {
            (0, 0, 0, 0): 13,
            (0, 0, 0, 1): 38,
            (0, 0, 1, 0): 34,
            (0, 0, 1, 1): 146,
            (0, 1, 0, 0): 49,
            (0, 1, 0, 1): 140,
            (0, 1, 1, 0): 120,
            (0, 1, 1, 1): 365,
            (1, 0, 0, 0): 1,
            (1, 0, 0, 1): 1,
            (1, 0, 1, 0): 7,
            (1, 0, 1, 1): 17,
            (1, 1, 0, 0): 3,
            (1, 1, 0, 1): 13,
            (1, 1, 1, 0): 13,
            (1, 1, 1, 1): 40,
        },
        {
            (0, 0, 0, 0): 175,
            (0, 0, 0, 1): 21,
            (0, 0, 1, 0): 23,
            (0, 0, 1, 1): 4,
            (0, 1, 0, 0): 26,
            (0, 1, 1, 0): 1,
            (1, 0, 0, 0): 548,
            (1, 0, 0, 1): 65,
            (1, 0, 1, 0): 68,
            (1, 0, 1, 1): 5,
            (1, 1, 0, 0): 53,
            (1, 1, 0, 1): 4,
            (1, 1, 1, 0): 6,
            (1, 1, 1, 1): 1,
        },
        {
            (0, 0, 0, 0): 52,
            (0, 0, 0, 1): 154,
            (0, 0, 1, 0): 6,
            (0, 0, 1, 1): 17,
            (0, 1, 0, 0): 8,
            (0, 1, 0, 1): 11,
            (0, 1, 1, 0): 1,
            (0, 1, 1, 1): 1,
            (1, 0, 0, 0): 169,
            (1, 0, 0, 1): 444,
            (1, 0, 1, 0): 16,
            (1, 0, 1, 1): 56,
            (1, 1, 0, 0): 15,
            (1, 1, 0, 1): 39,
            (1, 1, 1, 0): 3,
            (1, 1, 1, 1): 8,
        },
        {
            (0, 0, 0, 0): 52,
            (0, 0, 0, 1): 1,
            (0, 0, 1, 0): 153,
            (0, 0, 1, 1): 13,
            (0, 1, 0, 0): 8,
            (0, 1, 0, 1): 2,
            (0, 1, 1, 0): 14,
            (0, 1, 1, 1): 1,
            (1, 0, 0, 0): 160,
            (1, 0, 0, 1): 17,
            (1, 0, 1, 0): 446,
            (1, 0, 1, 1): 49,
            (1, 1, 0, 0): 12,
            (1, 1, 0, 1): 1,
            (1, 1, 1, 0): 66,
            (1, 1, 1, 1): 5,
        },
        {
            (0, 0, 0, 0): 13,
            (0, 0, 0, 1): 59,
            (0, 0, 1, 0): 42,
            (0, 0, 1, 1): 120,
            (0, 1, 0, 0): 2,
            (0, 1, 0, 1): 2,
            (0, 1, 1, 0): 3,
            (0, 1, 1, 1): 9,
            (1, 0, 0, 0): 44,
            (1, 0, 0, 1): 122,
            (1, 0, 1, 0): 126,
            (1, 0, 1, 1): 373,
            (1, 1, 0, 1): 13,
            (1, 1, 1, 0): 15,
            (1, 1, 1, 1): 57,
        },
        {
            (0, 0, 0, 0): 52,
            (0, 0, 0, 1): 7,
            (0, 0, 1, 0): 6,
            (0, 1, 0, 0): 149,
            (0, 1, 0, 1): 14,
            (0, 1, 1, 0): 14,
            (0, 1, 1, 1): 2,
            (1, 0, 0, 0): 163,
            (1, 0, 0, 1): 15,
            (1, 0, 1, 0): 25,
            (1, 0, 1, 1): 2,
            (1, 1, 0, 0): 459,
            (1, 1, 0, 1): 40,
            (1, 1, 1, 0): 46,
            (1, 1, 1, 1): 6,
        },
        {
            (0, 0, 0, 0): 16,
            (0, 0, 0, 1): 39,
            (0, 0, 1, 0): 3,
            (0, 0, 1, 1): 5,
            (0, 1, 0, 0): 38,
            (0, 1, 0, 1): 115,
            (0, 1, 1, 0): 6,
            (0, 1, 1, 1): 16,
            (1, 0, 0, 0): 53,
            (1, 0, 0, 1): 120,
            (1, 0, 1, 0): 3,
            (1, 0, 1, 1): 13,
            (1, 1, 0, 0): 108,
            (1, 1, 0, 1): 399,
            (1, 1, 1, 0): 19,
            (1, 1, 1, 1): 47,
        },
        {
            (0, 0, 0, 0): 14,
            (0, 0, 0, 1): 3,
            (0, 0, 1, 0): 42,
            (0, 0, 1, 1): 2,
            (0, 1, 0, 0): 33,
            (0, 1, 0, 1): 1,
            (0, 1, 1, 0): 137,
            (0, 1, 1, 1): 9,
            (1, 0, 0, 0): 47,
            (1, 0, 0, 1): 7,
            (1, 0, 1, 0): 114,
            (1, 0, 1, 1): 7,
            (1, 1, 0, 0): 131,
            (1, 1, 0, 1): 22,
            (1, 1, 1, 0): 397,
            (1, 1, 1, 1): 34,
        },
        {
            (0, 0, 0, 0): 1,
            (0, 0, 0, 1): 16,
            (0, 0, 1, 0): 17,
            (0, 0, 1, 1): 22,
            (0, 1, 0, 0): 3,
            (0, 1, 0, 1): 25,
            (0, 1, 1, 0): 31,
            (0, 1, 1, 1): 103,
            (1, 0, 0, 0): 20,
            (1, 0, 0, 1): 31,
            (1, 0, 1, 0): 49,
            (1, 0, 1, 1): 85,
            (1, 1, 0, 0): 41,
            (1, 1, 0, 1): 106,
            (1, 1, 1, 0): 110,
            (1, 1, 1, 1): 340,
        },
    ]
    qbs = [Node("qx", i) for i in range(9)]
    arc = Architecture([[qbs[i], qbs[i + 1]] for i in range(8)] + [[qbs[0], qbs[4]]])

    testc = Circuit(4, 4).H(0).CX(0, 3).CX(1, 2).CX(0, 1).CX(3, 2).measure_all()
    routed = testc.copy()
    mm = MappingManager(arc)
    mm.route_circuit(routed, [LexiLabellingMethod(), LexiRouteRoutingMethod()])
    DelayMeasures().apply(routed)
    readout = routed.qubit_readout

    spam = SpamCorrecter([list(readout.keys())])
    calib_circs = spam.calibration_circuits()

    assert len(calib_circs) == 16

    calib_counters = [
        Counter({OutcomeArray.from_readouts([key]): val for key, val in r.items()})
        for r in calib_results
    ]
    calib_brs = [BackendResult(counts=c) for c in calib_counters]
    spam.calculate_matrices(calib_brs)

    def check_correction(counts):  # type: ignore
        return (
            counts[(0, 0, 0, 0)] > raw_res[(0, 0, 0, 0)]
            and counts[(1, 1, 1, 1)] > raw_res[(1, 1, 1, 1)]
        )

    raw_counter = Counter(
        {OutcomeArray.from_readouts([key]): ceil(val) for key, val in raw_res.items()}
    )
    raw_br = BackendResult(counts=raw_counter)
    corrected_counts = spam.correct_counts(raw_br, spam.get_parallel_measure(routed))
    assert check_correction(compress_counts(corrected_counts.get_counts()))  # type: ignore


if __name__ == "__main__":
    test_spam_integration()
    test_spam_routing()
