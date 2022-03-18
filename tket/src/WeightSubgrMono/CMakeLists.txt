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

project(tket-${COMP})

if (NOT ${COMP} STREQUAL "WeightSubgrMono")
    message(FATAL_ERROR "Unexpected component name.")
endif()

add_library(tket-${COMP}
  Common/BitFunctions.cpp
  Common/DyadicFraction.cpp
  Common/RNG.cpp
  DomainInitialising/DistanceCounts.cpp
  DomainInitialising/DomainInitialiser.cpp
  EndToEndWrappers/CheckedWeightBounds.cpp
  EndToEndWrappers/MainSolver.cpp
  EndToEndWrappers/MainSolverData.cpp
  EndToEndWrappers/MainSolverParameters.cpp
  GraphTheoretic/CloseNeighboursCalculator.cpp
  GraphTheoretic/DerivedGraph.cpp
  GraphTheoretic/DerivedGraphs.cpp
  GraphTheoretic/DerivedGraphsCalculator.cpp
  GraphTheoretic/DerivedGraphsFilter.cpp
  GraphTheoretic/DerivedGraphsStorage.cpp
  GraphTheoretic/DerivedGraphsUpdater.cpp
  GraphTheoretic/FilterUtils.cpp
  GraphTheoretic/GeneralStructs.cpp
  GraphTheoretic/NeighboursData.cpp
  GraphTheoretic/TriangleCounts.cpp
  Reducing/AllDiffPropagator.cpp
  Reducing/CloseVerticesFilter.cpp
  Reducing/DerivedGraphsReducer.cpp
  Reducing/HallSetReducer.cpp
  Searching/CompleteTargetManager.cpp
  Searching/EnrichedNode.cpp
  Searching/FixedData.cpp
  Searching/SearchBranch.cpp
  Searching/SearchNode.cpp
  Searching/SearchNodeWrapper.cpp
  Searching/SharedData.cpp
  Searching/SolutionStorage.cpp
  Searching/SolutionWSM.cpp
  Searching/ValueOrdering.cpp
  Searching/VariableOrdering.cpp
  Searching/WeightUpdater.cpp
  WeightPruning/WeightNogoodDetector.cpp
  WeightPruning/WeightNogoodDetectorManager.cpp
  )

list(APPEND DEPS_${COMP}
    Utils)

foreach(DEP ${DEPS_${COMP}})
    target_include_directories(
        tket-${COMP} PRIVATE ${TKET_${DEP}_INCLUDE_DIR})
    target_link_libraries(
        tket-${COMP} PRIVATE tket-${DEP})
endforeach()

target_include_directories(tket-${COMP}
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${TKET_${COMP}_INCLUDE_DIR}
    ${TKET_${COMP}_INCLUDE_DIR}/${COMP})
