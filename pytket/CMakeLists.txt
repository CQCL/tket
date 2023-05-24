# Copyright 2019-2023 Cambridge Quantum Computing
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

cmake_minimum_required(VERSION 3.23)
project(pytket CXX)
add_definitions(-DPYBIND11_DETAILED_ERROR_MESSAGES)

find_package(tket CONFIG REQUIRED)
find_package(tklog CONFIG REQUIRED)
find_package(tkassert CONFIG REQUIRED)
find_package(symengine CONFIG REQUIRED)
find_package(gmp CONFIG REQUIRED)
find_package(pybind11_json CONFIG REQUIRED)
find_package(pybind11 CONFIG REQUIRED COMPONENTS headers)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_BINARY_DIR}")

set(Boost_NO_BOOST_CMAKE ON)

if (WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /WX")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -Wunreachable-code -Wunused")
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "(Apple)?Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
    # remove -Wno-deprecated-declarations once https://github.com/boostorg/boost/issues/688 is resolved
endif()
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # Suppress warnings coming from headers in pybind11_json/0.2.13:
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-redundant-move")
endif()

if (UNIX)
    # Allow binder libraries to load other shared libraries from same directory.
    set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
    set(CMAKE_INSTALL_RPATH "\${ORIGIN}")
endif()

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
    tket::tket pybind11::headers pybind11_json::pybind11_json)
if (WIN32)
    target_link_libraries(circuit PRIVATE
        symengine gmp::gmpxx gmp::gmp tkassert::tkassert)
    # For boost::uuid:
    target_link_libraries(circuit PRIVATE bcrypt)
endif()

pybind11_add_module(mapping binders/mapping.cpp)
target_include_directories(mapping PRIVATE binders/include)
target_link_libraries(mapping PRIVATE
    tket::tket pybind11::headers)

pybind11_add_module(transform binders/transform.cpp)
target_include_directories(transform PRIVATE binders/include)
target_link_libraries(transform PRIVATE
    tket::tket pybind11::headers)
if (WIN32)
    # For boost::uuid:
    target_link_libraries(transform PRIVATE bcrypt)
endif()

pybind11_add_module(predicates binders/predicates.cpp)
target_include_directories(predicates PRIVATE binders/include)
target_link_libraries(predicates PRIVATE
    tket::tket pybind11::headers pybind11_json::pybind11_json)
if (WIN32)
    # For boost::uuid:
    target_link_libraries(predicates PRIVATE bcrypt)
endif()

pybind11_add_module(passes binders/passes.cpp)
target_include_directories(passes PRIVATE binders/include)
target_link_libraries(passes PRIVATE
    tket::tket pybind11::headers pybind11_json::pybind11_json)
if (WIN32)
    target_link_libraries(passes PRIVATE
        symengine gmp::gmpxx gmp::gmp)
    # For boost::uuid:
    target_link_libraries(passes PRIVATE bcrypt)
endif()

pybind11_add_module(architecture binders/architecture.cpp)
target_include_directories(architecture PRIVATE binders/include)
target_link_libraries(architecture PRIVATE
    tket::tket pybind11::headers pybind11_json::pybind11_json)

pybind11_add_module(placement binders/placement.cpp)
target_include_directories(placement PRIVATE binders/include)
target_link_libraries(placement PRIVATE
    tket::tket pybind11::headers pybind11_json::pybind11_json)

pybind11_add_module(partition binders/partition.cpp)
target_include_directories(partition PRIVATE binders/include)
target_link_libraries(partition PRIVATE
    tket::tket pybind11_json::pybind11_json)
if (WIN32)
    # For boost::uuid:
    target_link_libraries(partition PRIVATE bcrypt)
endif()

pybind11_add_module(pauli binders/pauli.cpp)
target_include_directories(pauli PRIVATE binders/include)
target_link_libraries(pauli PRIVATE
    tket::tket pybind11::headers pybind11_json::pybind11_json)

pybind11_add_module(logging binders/logging.cpp)
target_include_directories(logging PRIVATE binders/include)
target_link_libraries(logging PRIVATE
    tket::tket pybind11::headers)

pybind11_add_module(utils_serialization binders/utils_serialization.cpp)
target_include_directories(utils_serialization PRIVATE binders/include)
target_link_libraries(utils_serialization PRIVATE
    tket::tket pybind11::headers pybind11_json::pybind11_json)

pybind11_add_module(tailoring binders/tailoring.cpp)
target_include_directories(tailoring PRIVATE binders/include)
target_link_libraries(tailoring PRIVATE
    tket::tket pybind11::headers)
if (WIN32)
    target_link_libraries(tailoring PRIVATE
        symengine gmp::gmpxx gmp::gmp)
    # For boost::uuid:
    target_link_libraries(tailoring PRIVATE bcrypt)
endif()

pybind11_add_module(tableau binders/tableau.cpp)
target_include_directories(tableau PRIVATE binders/include)
target_link_libraries(tableau PRIVATE
    tket::tket pybind11::headers)
if (WIN32)
    # For boost::uuid:
    target_link_libraries(tableau PRIVATE bcrypt)
endif()

pybind11_add_module(zx
    binders/zx/diagram.cpp
    binders/zx/rewrite.cpp)
target_include_directories(zx PRIVATE binders/include)
target_link_libraries(zx PRIVATE
    tket::tket pybind11::headers)
if (WIN32)
    target_link_libraries(zx PRIVATE
        symengine gmp::gmpxx gmp::gmp)
endif()

install(TARGETS
    logging
    utils_serialization
    circuit
    passes
    predicates
    partition
    pauli
    mapping
    transform
    tailoring
    tableau
    zx
    placement
    architecture
    LIBRARY DESTINATION lib)
