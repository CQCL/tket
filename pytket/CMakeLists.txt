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

set (CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

cmake_minimum_required(VERSION 3.16.0)

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
endif()

project(pytket)

find_file(CONANBUILDINFO_FILE conanbuildinfo.cmake HINTS ${CMAKE_BINARY_DIR})
include(${CONANBUILDINFO_FILE})
conan_basic_setup(NO_OUTPUT_DIRS)

set(Boost_NO_BOOST_CMAKE ON)

if (WIN32)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /WX")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror")
endif()

function(build_module module)
    pybind11_add_module(${module} NO_EXTRAS ${ARGN})
    target_include_directories(${module} PRIVATE binders/include)
    target_link_libraries(${module} PRIVATE ${CONAN_LIBS})
endfunction(build_module)

build_module(circuit
    binders/circuit/main.cpp
    binders/circuit/unitid.cpp
    binders/circuit/boxes.cpp
    binders/circuit/classical.cpp
    binders/circuit/Circuit/main.cpp
    binders/circuit/Circuit/add_op.cpp
    binders/circuit/Circuit/add_classical_op.cpp)
build_module(routing binders/routing.cpp)
build_module(transform binders/transform.cpp)
build_module(predicates binders/predicates.cpp)
build_module(passes binders/passes.cpp)
build_module(program binders/program.cpp)
build_module(partition binders/partition.cpp)
build_module(pauli binders/pauli.cpp)
build_module(logging binders/logging.cpp)
build_module(utils_serialization binders/utils_serialization.cpp)
build_module(tailoring binders/tailoring.cpp)
build_module(zx
    binders/zx/diagram.cpp
    binders/zx/rewrite.cpp)
