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

project(test-tktokenswap)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
endif()

find_file(CONANBUILDINFO_FILE conanbuildinfo.cmake HINTS ${CMAKE_BINARY_DIR})
include(${CONANBUILDINFO_FILE})
conan_basic_setup()

add_definitions(-DALL_LOGS)

add_executable(test-tktokenswap
    TableLookup/test_CanonicalRelabelling.cpp
    TableLookup/test_ExactMappingLookup.cpp
    TableLookup/test_FilteredSwapSequences.cpp
    TableLookup/test_SwapSequenceReductions.cpp
    TableLookup/test_SwapSequenceTable.cpp
    TestUtils/test_DebugFunctions.cpp
    TSAUtils/test_SwapFunctions.cpp
    test_SwapList.cpp
    test_SwapListOptimiser.cpp
    test_VectorListHybrid.cpp
    test_VectorListHybridSkeleton.cpp
    # test utilities:
    Data/FixedCompleteSolutions.cpp
    Data/FixedSwapSequences.cpp
    TableLookup/NeighboursFromEdges.cpp
    TableLookup/PermutationTestUtils.cpp
    TableLookup/SwapSequenceReductionTester.cpp
    TestUtils/ArchitectureEdgesReimplementation.cpp
    TestUtils/DebugFunctions.cpp
    TestUtils/DecodedProblemData.cpp
    TestUtils/GetRandomSet.cpp
    TestUtils/ProblemGeneration.cpp
    TestUtils/TestStatsStructs.cpp
)

target_include_directories(test-tktokenswap PRIVATE )

target_link_libraries(test-tktokenswap PRIVATE
    tklog tkassert tkrng tktokenswap ${CONAN_LIBS_CATCH2})

set(WITH_COVERAGE no CACHE BOOL "Link library with profiling for test coverage")
IF (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    IF (WITH_COVERAGE)
        target_link_options(test-tktokenswap PUBLIC --coverage)
    ENDIF()
ENDIF()
