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

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

cmake_minimum_required(VERSION 3.16.0)

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(tket-tests)

find_file(CONANBUILDINFO_FILE conanbuildinfo.cmake HINTS ${CMAKE_BINARY_DIR})
include(${CONANBUILDINFO_FILE})
conan_basic_setup()

IF (WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /WX /wd4068 /wd4267 /wd4805 /Wv:18 /EHsc /bigobj")
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS yes)
ENDIF()

add_definitions(-DALL_LOGS)


set(BUBBLE_TESTS_DIR ${CMAKE_CURRENT_SOURCE_DIR})

include(tkettestsfiles.cmake)

set(TEST_EXE test_tket)

add_executable(${TEST_EXE} ${TEST_SOURCES})
target_link_libraries(${TEST_EXE} ${CONAN_LIBS})