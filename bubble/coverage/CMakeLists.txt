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

# ONLY FOR
# TKET TEST COVERAGE

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

cmake_minimum_required(VERSION 3.16.0)

project(tket-test-coverage)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(CodeCoverage.cmake)

set(BUBBLE_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

include(${CMAKE_CURRENT_SOURCE_DIR}/src/tketCMake.cmake)

include(${CMAKE_CURRENT_SOURCE_DIR}/src/tketsrcfiles.cmake)

find_file(CONANBUILDINFO_FILE conanbuildinfo.cmake HINTS ${CMAKE_BINARY_DIR})
include(${CONANBUILDINFO_FILE})
conan_basic_setup()

# BOOST_ALLOW_DEPRECATED_HEADERS: Silence messages about deprecated headers during build.
set(CMAKE_CXX_FLAGS "-DBOOST_ALLOW_DEPRECATED_HEADERS")

set(Boost_USE_STATIC_LIBS ON)

add_definitions(-DALL_LOGS)

set(BUBBLE_TESTS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/tests)

include(tests/tkettestsfiles.cmake)

set(TEST_EXE test_tket)

add_executable(${TEST_EXE} ${TKET_SOURCES} ${TEST_SOURCES})

target_link_libraries(${TEST_EXE} ${CONAN_LIBS})

target_include_directories(${TEST_EXE} PRIVATE ${BUBBLE_SRC_DIR})

append_coverage_compiler_flags()

setup_target_for_coverage_gcovr_html(
    NAME coverage
    EXECUTABLE test_tket
    BASE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
