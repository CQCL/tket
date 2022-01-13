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

project(tket-simulation)

add_library(tket-simulation
    BitOperations.cpp
    CircuitSimulator.cpp
    DecomposeCircuit.cpp
    GateNode.cpp
    GateNodesBuffer.cpp
    PauliExpBoxUnitaryCalculator.cpp)

target_include_directories(tket-simulation
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${TKET_SIMULATION_INCLUDE_DIR}
    ${TKET_SIMULATION_INCLUDE_DIR}/Simulation
    ${TKET_UTILS_INCLUDE_DIR}
    ${TKET_OPTYPE_INCLUDE_DIR}
    ${TKET_OPS_INCLUDE_DIR}
    ${TKET_GATE_INCLUDE_DIR}
    ${TKET_CIRCUIT_INCLUDE_DIR})

target_link_libraries(tket-simulation
    PRIVATE
    tket-utils
    tket-optype
    tket-gate
    tket-circuit
    ${CONAN_LIBS})
