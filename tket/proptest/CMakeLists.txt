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
project(proptest-tket CXX)

find_package(Boost CONFIG REQUIRED)
find_package(gmp CONFIG)
if (NOT gmp_FOUND)
    find_package(PkgConfig REQUIRED)
    pkg_search_module(gmp REQUIRED IMPORTED_TARGET gmp)
endif()
find_package(SymEngine CONFIG REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)
find_package(tklog CONFIG REQUIRED)
find_package(tkassert CONFIG REQUIRED)
find_package(tkrng CONFIG REQUIRED)
find_package(tktokenswap CONFIG REQUIRED)
find_package(tkwsm CONFIG REQUIRED)
find_package(tket CONFIG REQUIRED)
find_package(rapidcheck CONFIG REQUIRED)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

IF (WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /WX /EHsc")
ELSE()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -Wunreachable-code -Wunused")
ENDIF()

if(CMAKE_CXX_COMPILER_ID MATCHES "(Apple)?Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
    # remove -Wno-deprecated-declarations once https://github.com/boostorg/boost/issues/688 is resolved
endif()

add_executable(proptest-tket
    src/ComparisonFunctions.cpp
    src/proptest.cpp
)

if (NOT TARGET gmp::gmp)
    add_library(gmp::gmp ALIAS PkgConfig::gmp)
endif()
if (NOT TARGET symengine::symengine)
    add_library(symengine::symengine ALIAS symengine)
endif()
if (NOT TARGET rapidcheck::rapidcheck)
    add_library(rapidcheck::rapidcheck ALIAS rapidcheck)
endif()

target_link_libraries(proptest-tket PRIVATE tket::tket)
target_link_libraries(proptest-tket PRIVATE Boost::headers)
target_link_libraries(proptest-tket PRIVATE gmp::gmp)
target_link_libraries(proptest-tket PRIVATE symengine::symengine)
target_link_libraries(proptest-tket PRIVATE Eigen3::Eigen)
target_link_libraries(proptest-tket PRIVATE rapidcheck::rapidcheck)

install(TARGETS proptest-tket DESTINATION "."
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        )
