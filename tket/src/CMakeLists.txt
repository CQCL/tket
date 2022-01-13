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

cmake_minimum_required(VERSION 3.16.0)
project(tket)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_file(CONANBUILDINFO_FILE conanbuildinfo.cmake HINTS ${CMAKE_BINARY_DIR})
include(${CONANBUILDINFO_FILE})
conan_basic_setup()

IF (WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /WX /EHsc")
ELSE()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -Wunreachable-code -Wunused")
ENDIF()
if(CMAKE_CXX_COMPILER_ID MATCHES "(Apple)?Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Winconsistent-missing-override -Wloop-analysis")
endif()

set(PROFILE_COVERAGE no CACHE BOOL "Build library with profiling for test coverage")
IF (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    IF (PROFILE_COVERAGE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -fprofile-arcs -ftest-coverage")
        # Bug in gcc 10: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=95353
        IF (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-stringop-overflow")
        ENDIF()
    ENDIF()
ENDIF()

IF (WIN32)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS yes)
ELSEIF(APPLE)
    # set correct install_name
    set(CMAKE_INSTALL_NAME_DIR "@loader_path")
    set(CMAKE_BUILD_WITH_INSTALL_NAME_DIR ON)
ENDIF()

set(TKET_UTILS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Utils/include)
set(TKET_ZX_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ZX/include)
set(TKET_OPTYPE_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/OpType/include)
set(TKET_CLIFFORD_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Clifford/include)
set(TKET_OPS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Ops/include)
set(TKET_GRAPHS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Graphs/include)
set(TKET_GATE_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Gate/include)
set(TKET_PAULIGRAPH_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/PauliGraph/include)
set(TKET_CIRCUIT_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Circuit/include)
set(TKET_ARCHITECTURE_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Architecture/include)
set(TKET_SIMULATION_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Simulation/include)
set(TKET_DIAGONALISATION_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Diagonalisation/include)
set(TKET_PROGRAM_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Program/include)
set(TKET_CHARACTERISATION_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Characterisation/include)
set(TKET_CONVERTERS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Converters/include)
set(TKET_ROUTING_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Routing/include)
set(TKET_MEASUREMENT_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/MeasurementSetup/include)
set(TKET_TRANSFORM_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Transformations/include)
set(TKET_AAS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ArchAwareSynth/include)
set(TKET_PREDICATES_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Predicates/include)

add_subdirectory(Utils)
add_subdirectory(ZX)
add_subdirectory(OpType)
add_subdirectory(Clifford)
add_subdirectory(Ops)
add_subdirectory(Graphs)
add_subdirectory(Gate)
add_subdirectory(PauliGraph)
add_subdirectory(Circuit)
add_subdirectory(Architecture)
add_subdirectory(Simulation)
add_subdirectory(Diagonalisation)
add_subdirectory(Program)
add_subdirectory(Characterisation)
add_subdirectory(Converters)
add_subdirectory(Routing)
add_subdirectory(MeasurementSetup)
add_subdirectory(Transformations)
add_subdirectory(ArchAwareSynth)
add_subdirectory(Predicates)
