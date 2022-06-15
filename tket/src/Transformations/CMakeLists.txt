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

project(tket-${COMP})

if (NOT ${COMP} STREQUAL "Transformations")
    message(FATAL_ERROR "Unexpected component name.")
endif()

add_library(tket-${COMP}
    Combinator.cpp
    Rebase.cpp
    BasicOptimisation.cpp
    PauliOptimisation.cpp
    CliffordOptimisation.cpp
    CliffordReductionPass.cpp
    OptimisationPass.cpp
    PhaseOptimisation.cpp
    ControlledGates.cpp
    Decomposition.cpp
    Replacement.cpp
    MeasurePass.cpp
    ContextualReduction.cpp
    ThreeQubitSquash.cpp
    SingleQubitSquash.cpp
    PhasedXFrontier.cpp
    PQPSquash.cpp
    StandardSquash.cpp)

list(APPEND DEPS_${COMP}
    Architecture
    Characterisation
    Circuit
    Clifford
    Converters
    Gate
    Graphs
    Ops
    OpType
    PauliGraph
    Utils
    ZX)

foreach(DEP ${DEPS_${COMP}})
    target_include_directories(
        tket-${COMP} PRIVATE ${TKET_${DEP}_INCLUDE_DIR})
    target_link_libraries(
        tket-${COMP} PRIVATE tket-${DEP})
endforeach()

target_include_directories(tket-${COMP}
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${TKET_${COMP}_INCLUDE_DIR}
    ${TKET_${COMP}_INCLUDE_DIR}/${COMP})

target_link_libraries(tket-${COMP} PRIVATE ${CONAN_LIBS_SYMENGINE})
