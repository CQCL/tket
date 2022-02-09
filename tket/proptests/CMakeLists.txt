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
project(tket-proptests)

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
    # /Zc:__cplusplus flag needed for C++20 support: see
    # https://github.com/emil-e/rapidcheck/issues/266
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /Zc:__cplusplus")
ENDIF()

add_executable(proptest
    ComparisonFunctions.cpp
    proptest.cpp
)

target_link_libraries(proptest PRIVATE
    rapidcheck
    tket-Architecture
    tket-Circuit
    tket-Gate
    tket-Ops
    tket-OpType
    tket-Predicates
    tket-TokenSwapping
    tket-Mapping
    tket-Simulation
    tket-Transformations
    tket-Utils)

target_link_libraries(proptest PRIVATE
    ${CONAN_LIBS_FMT} ${CONAN_LIBS_SPDLOG} ${CONAN_LIBS_SYMENGINE})
