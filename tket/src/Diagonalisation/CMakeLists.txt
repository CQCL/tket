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

project(tket-Diagonalisation)

add_library(tket-Diagonalisation
    DiagUtils.cpp
    Diagonalisation.cpp
    PauliPartition.cpp)

target_include_directories(tket-Diagonalisation
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${TKET_Diagonalisation_INCLUDE_DIR}
    ${TKET_Diagonalisation_INCLUDE_DIR}/Diagonalisation
    ${TKET_Utils_INCLUDE_DIR}
    ${TKET_OpType_INCLUDE_DIR}
    ${TKET_Clifford_INCLUDE_DIR}
    ${TKET_Ops_INCLUDE_DIR}
    ${TKET_Gate_INCLUDE_DIR}
    ${TKET_Graphs_INCLUDE_DIR}
    ${TKET_PauliGraph_INCLUDE_DIR}
    ${TKET_Circuit_INCLUDE_DIR})

target_link_libraries(tket-Diagonalisation
    PRIVATE
    tket-Utils
    tket-OpType
    tket-Ops
    tket-Gate
    tket-Graphs
    tket-PauliGraph
    tket-Circuit
    ${CONAN_LIBS})
