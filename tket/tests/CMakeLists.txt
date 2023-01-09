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

cmake_minimum_required(VERSION 3.16.0)
project(tket-tests)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_file(CONANBUILDINFO_FILE_TKET_TESTS conanbuildinfo.cmake HINTS ${CMAKE_CURRENT_BINARY_DIR})
include(${CONANBUILDINFO_FILE_TKET_TESTS})
conan_basic_setup()

IF (WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /WX /EHsc")
ELSE()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -Wunreachable-code -Wunused")
ENDIF()

if(CMAKE_CXX_COMPILER_ID MATCHES "(Apple)?Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
    # remove -Wno-deprecated-declarations once https://github.com/boostorg/boost/issues/688 is resolved
endif()

IF (WIN32)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS yes)
ENDIF()

add_definitions(-DALL_LOGS)

set(TKET_TESTS_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_definitions(-DTEST_CIRCUITS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/test_circuits")

include(tkettestutilsfiles.cmake)
include(tkettestsfiles.cmake)

add_executable(test_tket ${TESTUTILS_SOURCES} ${TEST_SOURCES})

target_link_libraries(test_tket PRIVATE ${CONAN_LIBS})
