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

import math
from pytest import fixture  # type: ignore
import random
import time
from typing import Any

from pytket.circuit import (  # type: ignore
    Circuit,
    OpType,
)


@fixture
def circuit(request: Any) -> Circuit:
    return Circuit(request.param)


@fixture
def circuit_4(request: Any, circuit: Circuit) -> Circuit:
    for nb_gates in range(request.param):
        circuit.X(0)
    return circuit


# TODO: Unify these two functions into one.
# Attemped to add them as methods of instance circuit using MethodType
# but circuit is a pybind11 object
def add_randomly_single_qubit_gates(
    nb_single_qubit_gates: int, single_qubit_indices: set, circuit: Circuit
) -> Circuit:
    for single_qubit_gate in range(nb_single_qubit_gates):
        # Randomly pick one qubit by index from the full list
        qubit_index = random.sample(single_qubit_indices, 1)[0]
        # Remove it from the list
        single_qubit_indices.remove(qubit_index)
        # Choose randomly between H- and S-gates
        if random.randrange(2) == 0:
            circuit.H(qubit_index)
        else:
            circuit.S(qubit_index)
    return circuit


def add_randomly_double_qubit_gates(
    nb_double_qubit_gates: int, double_qubit_indices: set, circuit: Circuit
) -> Circuit:
    for double_qubit_gate in range(nb_double_qubit_gates):
        # Randomly pick two qubit by index from the full list
        qubit_index_0 = random.sample(double_qubit_indices, 1)[0]
        # Remove it from the list
        double_qubit_indices.remove(qubit_index_0)
        qubit_index_1 = random.sample(double_qubit_indices, 1)[0]
        double_qubit_indices.remove(qubit_index_1)
        circuit.CX(qubit_index_0, qubit_index_1)
    return circuit


@fixture
def circuit_random(request: Any, circuit: Circuit) -> Circuit:
    single_qubit_gate_density = 0.5  # Proportion of single vs. double qubit gates
    double_qubit_gate_density = 0.5
    nb_single_qubit_gates = math.floor(single_qubit_gate_density * circuit.n_qubits)
    nb_double_qubit_gates = math.floor(
        double_qubit_gate_density * circuit.n_qubits * 0.5
    )
    for layer in range(request.param):
        # Recreate these sets every iteration as they get depopulated
        single_qubit_indices = set(i for i in range(circuit.n_qubits))
        double_qubit_indices = single_qubit_indices
        circuit = add_randomly_single_qubit_gates(
            nb_single_qubit_gates, single_qubit_indices, circuit
        )
        circuit = add_randomly_double_qubit_gates(
            nb_double_qubit_gates, double_qubit_indices, circuit
        )

    # Save the circuit to be run in tket benchmarks to a timestamped file
    timestr = time.strftime("%Y-%m-%d_%H:%M:%S")
    filename = (
        "circuit_random_nb_qubits="
        + str(circuit.n_qubits)
        + "_nb_layers="
        + str(request.param)
        + "_"
        + timestr
        + ".tkc"
    )
    circuit._to_tkc_file(filename)
    return circuit


@fixture
def optype() -> OpType:
    return OpType
