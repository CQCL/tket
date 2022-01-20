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

cmake_minimum_required(VERSION 3.16.0)
project(pytket)

set (CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
endif()

find_file(CONANBUILDINFO_FILE conanbuildinfo.cmake HINTS ${CMAKE_BINARY_DIR})
include(${CONANBUILDINFO_FILE})
conan_basic_setup(NO_OUTPUT_DIRS)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_BINARY_DIR}")
find_package(pybind11 REQUIRED)

set(Boost_NO_BOOST_CMAKE ON)

if (WIN32)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /WX")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror")
endif()

pybind11_add_module(circuit
    binders/circuit/main.cpp
    binders/circuit/unitid.cpp
    binders/circuit/boxes.cpp
    binders/circuit/classical.cpp
    binders/circuit/Circuit/main.cpp
    binders/circuit/Circuit/add_op.cpp
    binders/circuit/Circuit/add_classical_op.cpp)
target_include_directories(circuit PRIVATE binders/include)
target_link_libraries(circuit PRIVATE
    symengine
    tket-Circuit
    tket-Converters
    tket-Gate
    tket-Ops
    tket-OpType
    tket-Routing
    tket-Simulation
    tket-Utils)
if (WIN32)
    # For boost::uuid:
    target_link_libraries(circuit PRIVATE bcrypt)
endif()

pybind11_add_module(routing binders/routing.cpp)
target_include_directories(routing PRIVATE binders/include)
target_link_libraries(routing PRIVATE
    symengine
    tket-Architecture
    tket-Characterisation
    tket-Circuit
    tket-Gate
    tket-Graphs
    tket-Ops
    tket-OpType
    tket-Routing
    tket-Utils)
if (WIN32)
    # For boost::uuid:
    target_link_libraries(routing PRIVATE bcrypt)
endif()

pybind11_add_module(transform binders/transform.cpp)
target_include_directories(transform PRIVATE binders/include)
target_link_libraries(transform PRIVATE
    symengine
    tket-Architecture
    tket-Characterisation
    tket-Circuit
    tket-Clifford
    tket-Converters
    tket-Gate
    tket-Graphs
    tket-Ops
    tket-OpType
    tket-PauliGraph
    tket-Transformations
    tket-Utils)
if (WIN32)
    # For boost::uuid:
    target_link_libraries(transform PRIVATE bcrypt)
endif()

pybind11_add_module(predicates binders/predicates.cpp)
target_include_directories(predicates PRIVATE binders/include)
target_link_libraries(predicates PRIVATE
    symengine
    tket-ArchAwareSynth
    tket-Architecture
    tket-Characterisation
    tket-Circuit
    tket-Clifford
    tket-Converters
    tket-Gate
    tket-Graphs
    tket-Ops
    tket-OpType
    tket-PauliGraph
    tket-Predicates
    tket-Routing
    tket-Transformations
    tket-Utils)
if (WIN32)
    # For boost::uuid:
    target_link_libraries(predicates PRIVATE bcrypt)
endif()

pybind11_add_module(passes binders/passes.cpp)
target_include_directories(passes PRIVATE binders/include)
target_link_libraries(passes PRIVATE
    symengine
    tket-ArchAwareSynth
    tket-Architecture
    tket-Characterisation
    tket-Circuit
    tket-Clifford
    tket-Converters
    tket-Gate
    tket-Graphs
    tket-Ops
    tket-OpType
    tket-PauliGraph
    tket-Predicates
    tket-Routing
    tket-Transformations
    tket-Utils)
if (WIN32)
    # For boost::uuid:
    target_link_libraries(passes PRIVATE bcrypt)
endif()

pybind11_add_module(program binders/program.cpp)
target_include_directories(program PRIVATE binders/include)
target_link_libraries(program PRIVATE
    symengine
    tket-Circuit
    tket-Gate
    tket-Ops
    tket-OpType
    tket-Program
    tket-Utils)
if (WIN32)
    # For boost::uuid:
    target_link_libraries(program PRIVATE bcrypt)
endif()

pybind11_add_module(partition binders/partition.cpp)
target_include_directories(partition PRIVATE binders/include)
target_link_libraries(partition PRIVATE
    symengine
    tket-Circuit
    tket-Clifford
    tket-Converters
    tket-Diagonalisation
    tket-Gate
    tket-MeasurementSetup
    tket-Ops
    tket-OpType
    tket-PauliGraph
    tket-Utils)
if (WIN32)
    # For boost::uuid:
    target_link_libraries(partition PRIVATE bcrypt)
endif()

pybind11_add_module(pauli binders/pauli.cpp)
target_include_directories(pauli PRIVATE binders/include)
target_link_libraries(pauli PRIVATE
    symengine
    tket-Utils)

pybind11_add_module(logging binders/logging.cpp)
target_include_directories(logging PRIVATE binders/include)
target_link_libraries(logging PRIVATE
    symengine
    tket-Utils)

pybind11_add_module(utils_serialization binders/utils_serialization.cpp)
target_include_directories(utils_serialization PRIVATE binders/include)
target_link_libraries(utils_serialization PRIVATE
    symengine
    tket-Utils)

pybind11_add_module(tailoring binders/tailoring.cpp)
target_include_directories(tailoring PRIVATE binders/include)
target_link_libraries(tailoring PRIVATE
    symengine
    tket-Architecture
    tket-Characterisation
    tket-Circuit
    tket-Clifford
    tket-Converters
    tket-Diagonalisation
    tket-Gate
    tket-Graphs
    tket-Ops
    tket-OpType
    tket-PauliGraph
    tket-Utils)
if (WIN32)
    # For boost::uuid:
    target_link_libraries(tailoring PRIVATE bcrypt)
endif()

pybind11_add_module(tableau binders/tableau.cpp)
target_include_directories(tableau PRIVATE binders/include)
target_link_libraries(tableau PRIVATE
    symengine
    tket-Circuit
    tket-Clifford
    tket-Converters
    tket-Diagonalisation
    tket-Gate
    tket-Ops
    tket-OpType
    tket-PauliGraph
    tket-Utils)
if (WIN32)
    # For boost::uuid:
    target_link_libraries(tableau PRIVATE bcrypt)
endif()

pybind11_add_module(zx
    binders/zx/diagram.cpp
    binders/zx/rewrite.cpp)
target_include_directories(zx PRIVATE binders/include)
target_link_libraries(zx PRIVATE
    symengine
    tket-Utils
    tket-ZX)
