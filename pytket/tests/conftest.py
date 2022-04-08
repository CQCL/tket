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

import os
import math
from pytest import fixture, yield_fixture  # type: ignore
import random
import time
from typing import Any, Generator

from pytket.circuit import (  # type: ignore
    Circuit,
    OpType,
)
from pytket.qir.qir import circuit_to_qir, ExtendedModule, QUANTINUUM_GATES


@fixture
def file_name() -> str:
    return "SimpleCircuit.ll"


@fixture
def ext_module_quantinuum_gateset() -> ExtendedModule:
    em = ExtendedModule(
        name="Simple module for Quantinuum gateset.",
        num_qubits=2,
        num_results=1,
        gateset=QUANTINUUM_GATES,
    )
    em.module.builder.call(em.h, [em.module.qubits[0]])  # type: ignore
    em.module.builder.call(em.x, [em.module.qubits[1]])  # type: ignore
    em.module.builder.call(em.y, [em.module.qubits[0]])  # type: ignore
    em.module.builder.call(em.z, [em.module.qubits[1]])  # type: ignore
    em.module.builder.call(em.rx, [0.0, em.module.qubits[1]])  # type: ignore
    em.module.builder.call(em.ry, [1.0, em.module.qubits[0]])  # type: ignore
    em.module.builder.call(em.rz, [2.0, em.module.qubits[1]])  # type: ignore
    em.module.builder.call(em.phx, [1.0, 2.0, em.module.qubits[1]])  # type: ignore
    em.module.builder.call(
        em.cnot, [em.module.qubits[0], em.module.qubits[1]]  # type: ignore
    )
    em.module.builder.call(
        em.zzmax, [em.module.qubits[1], em.module.qubits[0]]  # type: ignore
    )
    em.module.builder.call(
        em.zzph, [1.0, em.module.qubits[0], em.module.qubits[1]]  # type: ignore
    )
    em.module.builder.call(
        em.mz, [em.module.qubits[0], em.module.results[0]]  # type: ignore
    )
    return em


@fixture
def circuit_quantinuum_gateset(file_name: str) -> Generator:
    c = Circuit(2, 2)
    c.H(0)
    c.X(1)
    c.Y(0)
    c.Z(1)
    c.CX(0, 1)
    c.add_gate(OpType.ZZMax, [1, 0])
    c.Rx(0.0, 1)
    c.Ry(1.0, 0)
    c.Rz(2.0, 1)
    c.add_gate(OpType.PhasedX, [1.5, 2.5], [1])
    c.add_gate(OpType.ZZPhase, [1.0], [0, 1])
    c.Measure(1, 1)
    circuit_to_qir(c, file_name, QUANTINUUM_GATES)
    yield
    os.remove(file_name)


@fixture
def circuit_pyqir_gateset(file_name: str) -> Generator:
    c = Circuit(2, 2)
    c.H(0)
    c.X(1)
    c.Y(0)
    c.Z(1)
    c.S(0)
    c.Sdg(1)
    c.T(0)
    c.Tdg(1)
    c.add_gate(OpType.Reset, [0])
    c.CX(0, 1)
    c.CZ(1, 0)
    c.Rx(0.0, 1)
    c.Ry(1.0, 0)
    c.Rz(2.0, 1)
    c.Measure(1, 1)
    circuit_to_qir(c, file_name)
    yield
    os.remove(file_name)


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
