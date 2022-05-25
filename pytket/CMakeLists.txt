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

list(APPEND TKET_EXTRA_LIBS
    ${CONAN_LIBS_FMT} ${CONAN_LIBS_SPDLOG} ${CONAN_LIBS_SYMENGINE})

pybind11_add_module(circuit
    binders/circuit/main.cpp
    binders/circuit/unitid.cpp
    binders/circuit/boxes.cpp
    binders/circuit/classical.cpp
    binders/circuit/library.cpp
    binders/circuit/Circuit/main.cpp
    binders/circuit/Circuit/add_op.cpp
    binders/circuit/Circuit/add_classical_op.cpp)
target_include_directories(circuit PRIVATE binders/include)
target_link_libraries(circuit PRIVATE
    tket-Circuit
    tket-Converters
    tket-Gate
    tket-Ops
    tket-OpType
    tket-Mapping
    tket-TokenSwapping
    tket-Simulation
    tket-Utils)
target_link_libraries(circuit PRIVATE ${TKET_EXTRA_LIBS})
if (WIN32)
    # For boost::uuid:
    target_link_libraries(circuit PRIVATE bcrypt)
endif()

pybind11_add_module(mapping binders/mapping.cpp)
target_include_directories(mapping PRIVATE binders/include)
target_link_libraries(mapping PRIVATE
    tket-ArchAwareSynth
    tket-Architecture
    tket-Characterisation
    tket-Circuit
    tket-Gate
    tket-Graphs
    tket-Mapping
    tket-Ops
    tket-OpType
    tket-PauliGraph
    tket-Placement    
    tket-Transformations
    tket-TokenSwapping
    tket-Utils)

pybind11_add_module(transform binders/transform.cpp)
target_include_directories(transform PRIVATE binders/include)
target_link_libraries(transform PRIVATE
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
target_link_libraries(transform PRIVATE ${TKET_EXTRA_LIBS})
if (WIN32)
    # For boost::uuid:
    target_link_libraries(transform PRIVATE bcrypt)
endif()

pybind11_add_module(predicates binders/predicates.cpp)
target_include_directories(predicates PRIVATE binders/include)
target_link_libraries(predicates PRIVATE
    tket-ArchAwareSynth
    tket-Architecture
    tket-Characterisation
    tket-Circuit
    tket-Clifford
    tket-Converters
    tket-Gate
    tket-Graphs
    tket-Mapping
    tket-Ops
    tket-OpType
    tket-PauliGraph
    tket-Predicates
    tket-TokenSwapping
    tket-Transformations
    tket-Utils)
target_link_libraries(predicates PRIVATE ${TKET_EXTRA_LIBS})
if (WIN32)
    # For boost::uuid:
    target_link_libraries(predicates PRIVATE bcrypt)
endif()

pybind11_add_module(passes binders/passes.cpp)
target_include_directories(passes PRIVATE binders/include)
target_link_libraries(passes PRIVATE
    tket-ArchAwareSynth
    tket-Architecture
    tket-Characterisation
    tket-Circuit
    tket-Clifford
    tket-Converters
    tket-Gate
    tket-Graphs
    tket-Mapping
    tket-Ops
    tket-OpType
    tket-PauliGraph
    tket-Predicates
    tket-TokenSwapping
    tket-Transformations
    tket-Utils)
target_link_libraries(passes PRIVATE ${TKET_EXTRA_LIBS})
if (WIN32)
    # For boost::uuid:
    target_link_libraries(passes PRIVATE bcrypt)
endif()


pybind11_add_module(architecture binders/architecture.cpp)
target_include_directories(architecture PRIVATE binders/include)
target_link_libraries(architecture PRIVATE
    tket-Architecture
    tket-Graphs
    tket-Utils)
target_link_libraries(architecture PRIVATE ${TKET_EXTRA_LIBS})


pybind11_add_module(placement binders/placement.cpp)
target_include_directories(placement PRIVATE binders/include)
target_link_libraries(placement PRIVATE
    tket-Placement
    tket-Architecture
    tket-Circuit
    tket-Utils)
target_link_libraries(placement PRIVATE ${TKET_EXTRA_LIBS})


pybind11_add_module(partition binders/partition.cpp)
target_include_directories(partition PRIVATE binders/include)
target_link_libraries(partition PRIVATE
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
target_link_libraries(partition PRIVATE ${TKET_EXTRA_LIBS})
if (WIN32)
    # For boost::uuid:
    target_link_libraries(partition PRIVATE bcrypt)
endif()

pybind11_add_module(pauli binders/pauli.cpp)
target_include_directories(pauli PRIVATE binders/include)
target_link_libraries(pauli PRIVATE
    tket-Utils)
target_link_libraries(pauli PRIVATE ${TKET_EXTRA_LIBS})

pybind11_add_module(logging binders/logging.cpp)
target_include_directories(logging PRIVATE binders/include)
target_link_libraries(logging PRIVATE
    tket-Utils)
target_link_libraries(logging PRIVATE ${TKET_EXTRA_LIBS})

pybind11_add_module(utils_serialization binders/utils_serialization.cpp)
target_include_directories(utils_serialization PRIVATE binders/include)
target_link_libraries(utils_serialization PRIVATE
    tket-Utils)
target_link_libraries(utils_serialization PRIVATE ${TKET_EXTRA_LIBS})

pybind11_add_module(tailoring binders/tailoring.cpp)
target_include_directories(tailoring PRIVATE binders/include)
target_link_libraries(tailoring PRIVATE
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
target_link_libraries(tailoring PRIVATE ${TKET_EXTRA_LIBS})
if (WIN32)
    # For boost::uuid:
    target_link_libraries(tailoring PRIVATE bcrypt)
endif()

pybind11_add_module(tableau binders/tableau.cpp)
target_include_directories(tableau PRIVATE binders/include)
target_link_libraries(tableau PRIVATE
    tket-Circuit
    tket-Clifford
    tket-Converters
    tket-Diagonalisation
    tket-Gate
    tket-Ops
    tket-OpType
    tket-PauliGraph
    tket-Utils)
target_link_libraries(tableau PRIVATE ${TKET_EXTRA_LIBS})
if (WIN32)
    # For boost::uuid:
    target_link_libraries(tableau PRIVATE bcrypt)
endif()

pybind11_add_module(zx
    binders/zx/diagram.cpp
    binders/zx/rewrite.cpp)
target_include_directories(zx PRIVATE binders/include)
target_link_libraries(zx PRIVATE
    tket-Circuit
    tket-Converters
    tket-Ops
    tket-OpType
    tket-Utils
    tket-ZX)
target_link_libraries(zx PRIVATE ${TKET_EXTRA_LIBS})
