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
project(test-tktokenswap CXX)

find_package(tktokenswap CONFIG REQUIRED)
find_package(tkrng CONFIG REQUIRED)
find_package(Catch2 CONFIG REQUIRED)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
endif()

add_executable(test-tktokenswap
    src/TableLookup/test_CanonicalRelabelling.cpp
    src/TableLookup/test_ExactMappingLookup.cpp
    src/TableLookup/test_FilteredSwapSequences.cpp
    src/TableLookup/test_SwapSequenceReductions.cpp
    src/TableLookup/test_SwapSequenceTable.cpp
    src/TestUtils/test_DebugFunctions.cpp
    src/TSAUtils/test_SwapFunctions.cpp
    src/test_SwapList.cpp
    src/test_SwapListOptimiser.cpp
    src/test_VectorListHybrid.cpp
    src/test_VectorListHybridSkeleton.cpp
    src/Data/FixedCompleteSolutions.cpp
    src/Data/FixedSwapSequences.cpp
    src/TableLookup/NeighboursFromEdges.cpp
    src/TableLookup/PermutationTestUtils.cpp
    src/TableLookup/SwapSequenceReductionTester.cpp
    src/TestUtils/ArchitectureEdgesReimplementation.cpp
    src/TestUtils/DebugFunctions.cpp
    src/TestUtils/DecodedProblemData.cpp
    src/TestUtils/GetRandomSet.cpp
    src/TestUtils/ProblemGeneration.cpp
    src/TestUtils/TestStatsStructs.cpp)

target_link_libraries(test-tktokenswap PRIVATE tktokenswap::tktokenswap)
target_link_libraries(test-tktokenswap PRIVATE tkrng::tkrng)
target_link_libraries(test-tktokenswap PRIVATE Catch2::Catch2WithMain)

set(WITH_COVERAGE no CACHE BOOL "Link library with profiling for test coverage")
IF (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    IF (WITH_COVERAGE)
        target_link_options(test-tktokenswap PUBLIC --coverage)
    ENDIF()
ENDIF()

install(TARGETS test-tktokenswap DESTINATION "."
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        )
