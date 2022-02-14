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

import pytest  # type: ignore


params_circuit_4 = [(4, gates) for gates in range(0, 1100, 100)]


@pytest.mark.parametrize("circuit,circuit_4", params_circuit_4, indirect=True)
def test_timing_circuit(benchmark, circuit_4, optype) -> None:  # type: ignore
    """Time depth_by_type X-gates on a basic 4-qubit circuit with single X-gate type"""
    benchmark.pedantic(
        circuit_4.depth_by_type, args=[optype.X], iterations=1000, rounds=4
    )


params_circuit_random = [
    (qubit, layer) for qubit, layer in zip(range(0, 110, 10), range(0, 1100, 100))
]


@pytest.mark.parametrize(
    "circuit, circuit_random", params_circuit_random, indirect=True  # type: ignore
)
def test_timing_random_circuit(benchmark, circuit_random, optype) -> None:
    """Time depth_by_type 2-qubit gates on a randomly generated circuit"""
    benchmark.pedantic(
        time_circuit.depth_by_type, args=[optype.CX], iterations=1000, rounds=4  # type: ignore
    )
