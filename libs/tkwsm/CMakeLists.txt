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
project(tkwsm CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
endif()

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

if(WIN32)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS yes)
endif()

if(APPLE)
    # set correct install_name
    set(CMAKE_INSTALL_NAME_DIR "@loader_path")
    set(CMAKE_BUILD_WITH_INSTALL_NAME_DIR ON)
endif()

add_library(tkwsm
  Common/BitFunctions.cpp
  Common/DyadicFraction.cpp
  EndToEndWrappers/MainSolver.cpp
  EndToEndWrappers/MainSolverParameters.cpp
  EndToEndWrappers/PreSearchComponents.cpp
  EndToEndWrappers/SolutionWSM.cpp
  GraphTheoretic/DerivedGraphs.cpp
  GraphTheoretic/DerivedGraphsCalculator.cpp
  GraphTheoretic/DomainInitialiser.cpp
  GraphTheoretic/FilterUtils.cpp
  GraphTheoretic/GeneralStructs.cpp
  GraphTheoretic/NearNeighboursData.cpp
  GraphTheoretic/NeighboursData.cpp
  GraphTheoretic/VertexRelabelling.cpp
  InitPlacement/EndToEndIQP.cpp
  InitPlacement/FastRandomBits.cpp
  InitPlacement/InputStructs.cpp
  InitPlacement/MonteCarloCompleteTargetSolution.cpp
  InitPlacement/MonteCarloManager.cpp
  InitPlacement/PrunedTargetEdges.cpp
  InitPlacement/SolutionJumper.cpp
  InitPlacement/UtilsIQP.cpp
  Reducing/DerivedGraphsReducer.cpp
  Reducing/DistancesReducer.cpp
  Reducing/HallSetReduction.cpp
  Reducing/ReducerWrapper.cpp
  Searching/DomainsAccessor.cpp
  Searching/NodeListTraversal.cpp
  Searching/NodesRawData.cpp
  Searching/SearchBranch.cpp
  Searching/ValueOrdering.cpp
  Searching/VariableOrdering.cpp
  Searching/WeightCalculator.cpp
  WeightPruning/WeightChecker.cpp
  WeightPruning/WeightNogoodDetector.cpp
  WeightPruning/WeightNogoodDetectorManager.cpp
  )
target_include_directories(tkwsm PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(tkwsm PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include/tkwsm)
target_include_directories(tkwsm INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(tkwsm PRIVATE
    ${CONAN_LIBS_TKLOG} ${CONAN_LIBS_TKASSERT} ${CONAN_LIBS_TKRNG})

if(MSVC)
  target_compile_options(tkwsm PRIVATE /W4 /WX  /wd4267)
else()
  target_compile_options(tkwsm PRIVATE -Wall -Wextra -Wpedantic -Werror)
endif()

set(PROFILE_COVERAGE no CACHE BOOL "Build library with profiling for test coverage")
IF (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    IF (PROFILE_COVERAGE)
        target_compile_options(tkwsm PRIVATE -g --coverage)
        target_link_options(tkwsm PUBLIC --coverage)
        # Bug in gcc 10: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=95353
        IF (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11)
            target_compile_options(tkwsm PRIVATE -Wno-stringop-overflow)
        ENDIF()
    ENDIF()
ENDIF()
