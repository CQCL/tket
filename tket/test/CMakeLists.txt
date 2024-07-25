# Copyright 2019-2024 Cambridge Quantum Computing
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
project(test-tket CXX)

find_package(tket CONFIG REQUIRED)
find_package(Catch2 CONFIG REQUIRED)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

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

FILE(GLOB_RECURSE JSON_FILES src/*.json)
FILE(GLOB_RECURSE TEST_SOURCES src/test_*.cpp)
SET(TEST_UTILS
    src/tests_main.cpp
    src/testutil.cpp
    src/CircuitsForTesting.cpp
    src/Graphs/EdgeSequence.cpp
    src/Graphs/EdgeSequenceColouringParameters.cpp
    src/Graphs/GraphTestingRoutines.cpp
    src/Graphs/RandomGraphGeneration.cpp
    src/Graphs/RandomPlanarGraphs.cpp
    src/TokenSwapping/Data/FixedCompleteSolutions.cpp
    src/TokenSwapping/Data/FixedSwapSequences.cpp
    src/TokenSwapping/TestUtils/ArchitectureEdgesReimplementation.cpp
    src/TokenSwapping/TestUtils/BestTsaTester.cpp
    src/TokenSwapping/TestUtils/DebugFunctions.cpp
    src/TokenSwapping/TestUtils/DecodedProblemData.cpp
    src/TokenSwapping/TestUtils/FullTsaTesting.cpp
    src/TokenSwapping/TestUtils/GetRandomSet.cpp
    src/TokenSwapping/TestUtils/PartialTsaTesting.cpp
    src/TokenSwapping/TestUtils/ProblemGeneration.cpp
    src/TokenSwapping/TestUtils/TestStatsStructs.cpp
    src/Gate/GatesData.cpp
    src/Simulation/ComparisonFunctions.cpp)

add_executable(test-tket ${TEST_UTILS} ${TEST_SOURCES})
target_link_libraries(test-tket PRIVATE tket::tket Catch2::Catch2WithMain)

set(WITH_COVERAGE no CACHE BOOL "Link library with profiling for test coverage")
IF (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    IF (WITH_COVERAGE)
        target_link_options(test-tket PUBLIC --coverage)
    ENDIF()
ENDIF()

install(FILES ${JSON_FILES} DESTINATION bin)

install(TARGETS test-tket DESTINATION "."
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        )
