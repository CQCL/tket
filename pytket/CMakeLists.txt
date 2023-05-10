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

find_package(Boost CONFIG REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)
find_package(tket CONFIG REQUIRED)
find_package(tklog CONFIG REQUIRED)
find_package(tkassert CONFIG REQUIRED)
find_package(tkrng CONFIG REQUIRED)
find_package(tktokenswap CONFIG REQUIRED)
find_package(tkwsm CONFIG REQUIRED)
find_package(tkassert CONFIG REQUIRED)
find_package(gmp CONFIG)
if (NOT gmp_FOUND)
    find_package(PkgConfig REQUIRED)
    pkg_search_module(gmp REQUIRED IMPORTED_TARGET gmp)
endif()
find_package(SymEngine CONFIG REQUIRED)
find_package(pybind11 CONFIG REQUIRED)
find_package(pybind11_json CONFIG REQUIRED)

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

if (NOT TARGET pybind11_json::pybind11_json)
    add_library(pybind11_json::pybind11_json ALIAS pybind11_json)
endif()
if (NOT TARGET gmp::gmp)
    add_library(gmp::gmp ALIAS PkgConfig::gmp)
endif()
if (NOT TARGET symengine::symengine)
    add_library(symengine::symengine ALIAS symengine)
endif()

list(APPEND lib_deps
    tket::tket
    tklog::tklog
    tkassert::tkassert
    pybind11::headers
    pybind11_json::pybind11_json
    Eigen3::Eigen
    gmp::gmp
    symengine::symengine)
if (WIN32)
    list(APPEND lib_deps bcrypt) # For boost::uuid
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
target_link_libraries(circuit PRIVATE ${lib_deps})

pybind11_add_module(mapping binders/mapping.cpp)
target_include_directories(mapping PRIVATE binders/include)
target_link_libraries(mapping PRIVATE ${lib_deps})

pybind11_add_module(transform binders/transform.cpp)
target_include_directories(transform PRIVATE binders/include)
target_link_libraries(transform PRIVATE ${lib_deps})

pybind11_add_module(predicates binders/predicates.cpp)
target_include_directories(predicates PRIVATE binders/include)
target_link_libraries(predicates PRIVATE ${lib_deps})

pybind11_add_module(passes binders/passes.cpp)
target_include_directories(passes PRIVATE binders/include)
target_link_libraries(passes PRIVATE ${lib_deps})

pybind11_add_module(architecture binders/architecture.cpp)
target_include_directories(architecture PRIVATE binders/include)
target_link_libraries(architecture PRIVATE ${lib_deps})

pybind11_add_module(placement binders/placement.cpp)
target_include_directories(placement PRIVATE binders/include)
target_link_libraries(placement PRIVATE ${lib_deps})

pybind11_add_module(partition binders/partition.cpp)
target_include_directories(partition PRIVATE binders/include)
target_link_libraries(partition PRIVATE ${lib_deps})

pybind11_add_module(pauli binders/pauli.cpp)
target_include_directories(pauli PRIVATE binders/include)
target_link_libraries(pauli PRIVATE ${lib_deps})

pybind11_add_module(logging binders/logging.cpp)
target_include_directories(logging PRIVATE binders/include)
target_link_libraries(logging PRIVATE ${lib_deps})

pybind11_add_module(utils_serialization binders/utils_serialization.cpp)
target_include_directories(utils_serialization PRIVATE binders/include)
target_link_libraries(utils_serialization PRIVATE ${lib_deps})

pybind11_add_module(tailoring binders/tailoring.cpp)
target_include_directories(tailoring PRIVATE binders/include)
target_link_libraries(tailoring PRIVATE ${lib_deps})

pybind11_add_module(tableau binders/tableau.cpp)
target_include_directories(tableau PRIVATE binders/include)
target_link_libraries(tableau PRIVATE ${lib_deps})

pybind11_add_module(zx
    binders/zx/diagram.cpp
    binders/zx/rewrite.cpp)
target_include_directories(zx PRIVATE binders/include)
target_link_libraries(zx PRIVATE ${lib_deps})

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
