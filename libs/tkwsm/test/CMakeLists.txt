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
project(test-tkwsm CXX)

find_package(tkwsm CONFIG REQUIRED)
find_package(tkassert CONFIG REQUIRED)
find_package(tkrng CONFIG REQUIRED)
find_package(Catch2 CONFIG REQUIRED)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
endif()

add_executable(test-tkwsm
    src/Common/test_BitFunctions.cpp
    src/Common/test_DyadicFraction.cpp
    src/Common/test_GeneralUtils.cpp
    src/Common/test_LogicalStack.cpp
    src/EndToEndWrappers/test_SolutionWSM.cpp
    src/GraphTheoretic/test_FilterUtils.cpp
    src/GraphTheoretic/test_GeneralStructs.cpp
    src/GraphTheoretic/test_NeighboursData.cpp
    src/InitPlacement/test_InitialPlacementProblems.cpp
    src/InitPlacement/test_MonteCarloCompleteTargetSolution.cpp
    src/InitPlacement/test_PrunedTargetEdges.cpp
    src/InitPlacement/test_WeightedBinaryTree.cpp
    src/InitPlacement/test_WeightedSquareGrid.cpp
    src/InitPlacement/PlacementCostModelInterface.cpp
    src/InitPlacement/TestUtilsIQP.cpp
    src/InitPlacement/TestWeightedGraphData.cpp
    src/InitPlacement/WeightedBinaryTree.cpp
    src/InitPlacement/WeightedSquareGrid.cpp
    src/Searching/test_NodesRawData.cpp
    src/Searching/test_NodesRawDataTraversals.cpp
    src/SolvingProblems/test_CubicLattice.cpp
    src/SolvingProblems/test_FixedSmallGraphs.cpp
    src/SolvingProblems/test_RandomGraphs.cpp
    src/SolvingProblems/test_SnakeIntoSquareGrid.cpp
    src/SolvingProblems/test_SpecialProblems.cpp
    src/SolvingProblems/test_SquareGrids.cpp
    src/SolvingProblems/test_UnweightedProblems.cpp
    src/SolvingProblems/test_UnweightedSelfEmbeddings.cpp
    src/TestUtils/test_SquareGridGeneration.cpp
    src/TestUtils/CheckedSolution.cpp
    src/TestUtils/FixedArchitectures.cpp
    src/TestUtils/GraphGeneration.cpp
    src/TestUtils/ProblemGeneration.cpp
    src/TestUtils/ResumedSolutionChecker.cpp
    src/TestUtils/SquareGridGeneration.cpp
    src/TestUtils/TestSettings.cpp)

target_link_libraries(test-tkwsm PRIVATE tkwsm::tkwsm)
target_link_libraries(test-tkwsm PRIVATE tkassert::tkassert)
target_link_libraries(test-tkwsm PRIVATE tkrng::tkrng)
target_link_libraries(test-tkwsm PRIVATE Catch2::Catch2WithMain)

set(WITH_COVERAGE no CACHE BOOL "Link library with profiling for test coverage")
IF (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    IF (WITH_COVERAGE)
        target_link_options(test-tkwsm PUBLIC --coverage)
    ENDIF()
ENDIF()

install(TARGETS test-tkwsm DESTINATION "."
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        )
