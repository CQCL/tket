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

from pytket import logging
from pytket.circuit import Circuit, Qubit  # type: ignore
from pytket.pauli import Pauli, QubitPauliString  # type: ignore
from pytket.partition import (  # type: ignore
    PauliPartitionStrat,
    MeasurementBitMap,
    MeasurementSetup,
    measurement_reduction,
)

import pytest  # type: ignore
import platform
from typing import Any


def test_empty_setup() -> None:
    ms = MeasurementSetup()
    assert ms.verify()


def test_parity_flip() -> None:
    ms = MeasurementSetup()
    circ = Circuit(1, 1)
    circ.X(0)
    circ.Measure(0, 0)
    ms.add_measurement_circuit(circ)

    tensor = dict()
    tensor[Qubit(0)] = Pauli.Z

    mbm = MeasurementBitMap(0, [0], True)
    string = QubitPauliString(tensor)
    ms.add_result_for_term(string, mbm)
    assert ms.verify()
    ms_str = ms.__repr__()
    assert "|| (Zq[0]) ||" in ms_str


def test_reduction() -> None:
    str1 = QubitPauliString({Qubit(0): Pauli.I, Qubit(1): Pauli.Z, Qubit(2): Pauli.X})
    str2 = QubitPauliString({Qubit(0): Pauli.Y, Qubit(1): Pauli.Z, Qubit(2): Pauli.X})
    str3 = QubitPauliString({Qubit(0): Pauli.X, Qubit(1): Pauli.Y, Qubit(2): Pauli.Y})

    strings = [str1, str2, str3]
    strats = [PauliPartitionStrat.NonConflictingSets, PauliPartitionStrat.CommutingSets]
    for s in strats:
        measurements = measurement_reduction(strings, s)
        assert len(measurements.measurement_circs) == 2
        assert measurements.verify()


# readouterr() doesn't seem to work correctly on Windows, so skip. TODO investigate.
@pytest.mark.skipif(platform.system() == "Windows", reason="issues with readouterr()")
def test_error_logging(capfd: Any) -> None:
    ms = MeasurementSetup()
    circ = Circuit(2, 2)
    circ.X(0)
    circ.Measure(0, 0)
    circ.V(1)
    circ.Measure(1, 1)
    ms.add_measurement_circuit(circ)
    zi = QubitPauliString()
    zi[Qubit(0)] = Pauli.Z
    mbm = MeasurementBitMap(0, [0], False)
    ms.add_result_for_term(zi, mbm)
    logging.set_level(logging.level.err)  # type: ignore
    assert not ms.verify()
    out = capfd.readouterr().out
    assert "[error]" in out
    logging.set_level(logging.level.critical)  # type: ignore
    assert not ms.verify()
    out = capfd.readouterr().out
    assert not out


def test_serialization() -> None:
    mbm_dict = {"circ_index": 0, "bits": [0], "invert": True}
    mbm = MeasurementBitMap(0, [0], True)
    j_mbm = mbm.to_dict()
    assert j_mbm == mbm_dict
    assert mbm_dict == MeasurementBitMap.from_dict(mbm_dict).to_dict()

    ms = MeasurementSetup()
    circ = Circuit(2, 2)
    circ.X(0)
    circ.Measure(0, 0)
    circ.V(1)
    circ.Measure(1, 1)
    circ2 = Circuit(2, 2)
    circ2.Measure(0, 0)
    circ2.Measure(1, 1)
    ms.add_measurement_circuit(circ)
    ms.add_measurement_circuit(circ2)
    zi = QubitPauliString()
    zi[Qubit(0)] = Pauli.Z
    iz = QubitPauliString()
    iz[Qubit(1)] = Pauli.Z
    xx = QubitPauliString()
    xx[Qubit(0)] = Pauli.X
    xx[Qubit(1)] = Pauli.X
    mbm = MeasurementBitMap(0, [0], False)
    mbm2 = MeasurementBitMap(1, [0], False)
    mbm3 = MeasurementBitMap(1, [0, 1], False)

    ms.add_result_for_term(zi, mbm)
    ms.add_result_for_term(iz, mbm)
    ms.add_result_for_term(zi, mbm2)
    ms.add_result_for_term(xx, mbm3)

    j_ms = ms.to_dict()
    assert len(j_ms["circs"]) == 2
    assert j_ms["circs"][0] == circ.to_dict()
    assert j_ms["circs"][1] == circ2.to_dict()
    assert len(j_ms["result_map"]) == 3
    assert j_ms["result_map"][0] == [iz.to_list(), [mbm.to_dict()]]
    assert j_ms["result_map"][1] == [xx.to_list(), [mbm3.to_dict()]]
    assert j_ms["result_map"][2] == [zi.to_list(), [mbm.to_dict(), mbm2.to_dict()]]
    assert MeasurementSetup.from_dict(j_ms).to_dict() == j_ms


if __name__ == "__main__":
    test_empty_setup()
    test_parity_flip()
    test_reduction()
    test_serialization()
    # test_error_logging()
