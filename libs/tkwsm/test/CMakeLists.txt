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

project(test-tkwsm)

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

add_executable(test-tkwsm
    Common/test_BitFunctions.cpp
    Common/test_DyadicFraction.cpp
    Common/test_GeneralUtils.cpp
    Common/test_LogicalStack.cpp
    EndToEndWrappers/test_SolutionWSM.cpp
    GraphTheoretic/test_FilterUtils.cpp
    GraphTheoretic/test_GeneralStructs.cpp
    GraphTheoretic/test_NeighboursData.cpp
    InitPlacement/test_InitialPlacementProblems.cpp
    InitPlacement/test_MonteCarloCompleteTargetSolution.cpp
    InitPlacement/test_PrunedTargetEdges.cpp
    InitPlacement/test_WeightedBinaryTree.cpp
    InitPlacement/test_WeightedSquareGrid.cpp
    InitPlacement/PlacementCostModelInterface.cpp
    InitPlacement/TestUtilsIQP.cpp
    InitPlacement/TestWeightedGraphData.cpp
    InitPlacement/WeightedBinaryTree.cpp
    InitPlacement/WeightedSquareGrid.cpp
    Searching/test_NodesRawData.cpp
    Searching/test_NodesRawDataTraversals.cpp
    SolvingProblems/test_CubicLattice.cpp
    SolvingProblems/test_FixedSmallGraphs.cpp
    SolvingProblems/test_RandomGraphs.cpp
    SolvingProblems/test_SnakeIntoSquareGrid.cpp
    SolvingProblems/test_SpecialProblems.cpp
    SolvingProblems/test_SquareGrids.cpp
    SolvingProblems/test_UnweightedProblems.cpp
    SolvingProblems/test_UnweightedSelfEmbeddings.cpp
    TestUtils/test_SquareGridGeneration.cpp
    TestUtils/CheckedSolution.cpp
    TestUtils/FixedArchitectures.cpp
    TestUtils/GraphGeneration.cpp
    TestUtils/ProblemGeneration.cpp
    TestUtils/ResumedSolutionChecker.cpp
    TestUtils/SquareGridGeneration.cpp
    TestUtils/TestSettings.cpp
)

target_include_directories(test-tkwsm PRIVATE )

target_link_libraries(test-tkwsm PRIVATE ${CONAN_LIBS})

set(WITH_COVERAGE no CACHE BOOL "Link library with profiling for test coverage")
IF (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    IF (WITH_COVERAGE)
        target_link_options(test-tkwsm PUBLIC --coverage)
    ENDIF()
ENDIF()
