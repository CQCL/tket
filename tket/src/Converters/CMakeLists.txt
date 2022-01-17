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

project(tket-converters)

add_library(tket-converters
    CliffTableauConverters.cpp
    PauliGadget.cpp
    PauliGraphConverters.cpp
    Gauss.cpp
    PhasePoly.cpp
    UnitaryTableauBox.cpp
    UnitaryTableauConverters.cpp)

target_include_directories(tket-converters
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${TKET_Converters_INCLUDE_DIR}
    ${TKET_Converters_INCLUDE_DIR}/Converters
    ${TKET_Utils_INCLUDE_DIR}
    ${TKET_OpType_INCLUDE_DIR}
    ${TKET_Clifford_INCLUDE_DIR}
    ${TKET_Ops_INCLUDE_DIR}
    ${TKET_Gate_INCLUDE_DIR}
    ${TKET_PauliGraph_INCLUDE_DIR}
    ${TKET_Circuit_INCLUDE_DIR}
    ${TKET_Diagonalisation_INCLUDE_DIR})

target_link_libraries(tket-converters
    PRIVATE
    tket-utils
    tket-optype
    tket-clifford
    tket-ops
    tket-gate
    tket-pauligraph
    tket-circuit
    tket-diagonalisation
    ${CONAN_LIBS})
