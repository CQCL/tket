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

set(TKET_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

project(tket-compiler)

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_file(CONANBUILDINFO_FILE conanbuildinfo.cmake HINTS ${CMAKE_BINARY_DIR})
include(${CONANBUILDINFO_FILE})
conan_basic_setup()

set(TKET "tket")

IF (WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /WX /EHsc")
ELSE()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -Wunreachable-code -Wunused")
ENDIF()
if(CMAKE_CXX_COMPILER_ID MATCHES "(Apple)?Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Winconsistent-missing-override -Wloop-analysis")
endif()

set(PROFILE_COVERAGE no CACHE BOOL "Build library with profiling for test coverage")
IF (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    IF (PROFILE_COVERAGE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -fprofile-arcs -ftest-coverage")
        # Bug in gcc 10: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=95353
        IF (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-stringop-overflow")
        ENDIF()
    ENDIF()
ENDIF()

IF (WIN32)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS yes)
ELSEIF(APPLE)
    # set correct install_name
    set(CMAKE_INSTALL_NAME_DIR "@loader_path")
    set(CMAKE_BUILD_WITH_INSTALL_NAME_DIR ON)
ENDIF()

set(TKET_ARCHITECTURE_DIR ${TKET_SRC_DIR}/Architecture)
set(TKET_AAS_DIR ${TKET_SRC_DIR}/ArchAwareSynth)
set(TKET_CIRCUIT_DIR ${TKET_SRC_DIR}/Circuit)
set(TKET_CLIFFORD_DIR ${TKET_SRC_DIR}/Clifford)
set(TKET_DIAGONALISATION_DIR ${TKET_SRC_DIR}/Diagonalisation)
set(TKET_GRAPH_DIR ${TKET_SRC_DIR}/Graphs)
set(TKET_UTILS_DIR ${TKET_SRC_DIR}/Utils)
set(TKET_OPTYPE_DIR ${TKET_SRC_DIR}/OpType)
set(TKET_OPS_DIR ${TKET_SRC_DIR}/Ops)
set(TKET_GATE_DIR ${TKET_SRC_DIR}/Gate)
set(TKET_SIMULATION_DIR ${TKET_SRC_DIR}/Simulation)
set(TKET_ROUTING_DIR ${TKET_SRC_DIR}/Routing)
set(TKET_MAPPING_DIR ${TKET_SRC_DIR}/Mapping)
set(TKET_TOKEN_SWAPPING_DIR ${TKET_SRC_DIR}/TokenSwapping)
set(TKET_TRANSFORM_DIR ${TKET_SRC_DIR}/Transformations)
set(TKET_CHARACTERISATION_DIR ${TKET_SRC_DIR}/Characterisation)
set(TKET_PREDS_DIR ${TKET_SRC_DIR}/Predicates)
set(TKET_PAULIGRAPH_DIR ${TKET_SRC_DIR}/PauliGraph)
set(TKET_CONVERTERS_DIR ${TKET_SRC_DIR}/Converters)
set(TKET_PROGRAM_DIR ${TKET_SRC_DIR}/Program)
set(TKET_MEASUREMENT_DIR ${TKET_SRC_DIR}/MeasurementSetup)
set(TKET_ZX_DIR ${TKET_SRC_DIR}/ZX)

set(TKET_SOURCES

    # OpType
    ${TKET_OPTYPE_DIR}/OpDesc.cpp
    ${TKET_OPTYPE_DIR}/OpTypeInfo.cpp
    ${TKET_OPTYPE_DIR}/OpTypeFunctions.cpp
    ${TKET_OPTYPE_DIR}/OpTypeJson.cpp

    # Ops
    ${TKET_OPS_DIR}/Conditional.cpp
    ${TKET_OPS_DIR}/FlowOp.cpp
    ${TKET_OPS_DIR}/MetaOp.cpp
    ${TKET_OPS_DIR}/Op.cpp
    ${TKET_OPS_DIR}/ClassicalOps.cpp
    ${TKET_OPS_DIR}/OpJsonFactory.cpp

    # Gate
    ${TKET_GATE_DIR}/Gate.cpp
    ${TKET_GATE_DIR}/GatePtr.cpp
    ${TKET_GATE_DIR}/OpPtrFunctions.cpp
    ${TKET_GATE_DIR}/Rotation.cpp
    ${TKET_GATE_DIR}/SymTable.cpp
    ${TKET_GATE_DIR}/GateUnitaryMatrix.cpp
    ${TKET_GATE_DIR}/GateUnitaryMatrixComposites.cpp
    ${TKET_GATE_DIR}/GateUnitaryMatrixError.cpp
    ${TKET_GATE_DIR}/GateUnitaryMatrixFixedMatrices.cpp
    ${TKET_GATE_DIR}/GateUnitaryMatrixPrimitives.cpp
    ${TKET_GATE_DIR}/GateUnitaryMatrixUtils.cpp
    ${TKET_GATE_DIR}/GateUnitaryMatrixVariableQubits.cpp
    ${TKET_GATE_DIR}/GateUnitarySparseMatrix.cpp

    # Circuit
    ${TKET_CIRCUIT_DIR}/Boxes.cpp
    ${TKET_CIRCUIT_DIR}/Circuit.cpp
    ${TKET_CIRCUIT_DIR}/CircuitJson.cpp
    ${TKET_CIRCUIT_DIR}/CommandJson.cpp
    ${TKET_CIRCUIT_DIR}/macro_manipulation.cpp
    ${TKET_CIRCUIT_DIR}/basic_circ_manip.cpp
    ${TKET_CIRCUIT_DIR}/latex_drawing.cpp
    ${TKET_CIRCUIT_DIR}/macro_circ_info.cpp
    ${TKET_CIRCUIT_DIR}/setters_and_getters.cpp
    ${TKET_CIRCUIT_DIR}/CircUtils.cpp
    ${TKET_CIRCUIT_DIR}/ThreeQubitConversion.cpp
    ${TKET_CIRCUIT_DIR}/AssertionSynthesis.cpp
    ${TKET_CIRCUIT_DIR}/CircPool.cpp
    ${TKET_CIRCUIT_DIR}/DAGProperties.cpp
    ${TKET_CIRCUIT_DIR}/OpJson.cpp

    # Simulation
    ${TKET_SIMULATION_DIR}/BitOperations.cpp
    ${TKET_SIMULATION_DIR}/CircuitSimulator.cpp
    ${TKET_SIMULATION_DIR}/DecomposeCircuit.cpp
    ${TKET_SIMULATION_DIR}/GateNode.cpp
    ${TKET_SIMULATION_DIR}/GateNodesBuffer.cpp
    ${TKET_SIMULATION_DIR}/PauliExpBoxUnitaryCalculator.cpp

    # Clifford
    ${TKET_CLIFFORD_DIR}/CliffTableau.cpp

    # Diagonalisation
    ${TKET_DIAGONALISATION_DIR}/DiagUtils.cpp
    ${TKET_DIAGONALISATION_DIR}/Diagonalisation.cpp
    ${TKET_DIAGONALISATION_DIR}/PauliPartition.cpp

    # Graph algorithms
    ${TKET_GRAPH_DIR}/AdjacencyData.cpp
    ${TKET_GRAPH_DIR}/BruteForceColouring.cpp
    ${TKET_GRAPH_DIR}/ColouringPriority.cpp
    ${TKET_GRAPH_DIR}/GraphColouring.cpp
    ${TKET_GRAPH_DIR}/GraphRoutines.cpp
    ${TKET_GRAPH_DIR}/LargeCliquesResult.cpp
    ${TKET_GRAPH_DIR}/ArticulationPoints.cpp

    # Token swapping
    ${TKET_TOKEN_SWAPPING_DIR}/TSAUtils/DebugFunctions.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/TSAUtils/DistanceFunctions.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/TSAUtils/GeneralFunctions.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/TSAUtils/SwapFunctions.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/TSAUtils/VertexMappingFunctions.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/TSAUtils/VertexSwapResult.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/ArchitectureMapping.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/BestFullTsa.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/CyclesCandidateManager.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/CyclesGrowthManager.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/CyclesPartialTsa.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/CyclicShiftCostEstimate.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/DistancesFromArchitecture.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/DistancesInterface.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/DynamicTokenTracker.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/HybridTsa00.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/main_entry_functions.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/NeighboursFromArchitecture.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/NeighboursInterface.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/PartialTsaInterface.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/PathFinderInterface.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/RiverFlowPathFinder.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/RNG.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/SwapListOptimiser.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/TrivialTSA.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/VectorListHybridSkeleton.cpp

    # Token swapping table lookup
    ${TKET_TOKEN_SWAPPING_DIR}/TableLookup/CanonicalRelabelling.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/TableLookup/ExactMappingLookup.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/TableLookup/FilteredSwapSequences.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/TableLookup/PartialMappingLookup.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/TableLookup/SwapConversion.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/TableLookup/SwapListSegmentOptimiser.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/TableLookup/SwapListTableOptimiser.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/TableLookup/SwapSequenceTable.cpp
    ${TKET_TOKEN_SWAPPING_DIR}/TableLookup/VertexMapResizing.cpp

    # Transformations
    ${TKET_TRANSFORM_DIR}/Combinator.cpp
    ${TKET_TRANSFORM_DIR}/Rebase.cpp
    ${TKET_TRANSFORM_DIR}/BasicOptimisation.cpp
    ${TKET_TRANSFORM_DIR}/PauliOptimisation.cpp
    ${TKET_TRANSFORM_DIR}/CliffordOptimisation.cpp
    ${TKET_TRANSFORM_DIR}/CliffordReductionPass.cpp
    ${TKET_TRANSFORM_DIR}/OptimisationPass.cpp
    ${TKET_TRANSFORM_DIR}/PhaseOptimisation.cpp
    ${TKET_TRANSFORM_DIR}/ControlledGates.cpp
    ${TKET_TRANSFORM_DIR}/Decomposition.cpp
    ${TKET_TRANSFORM_DIR}/Replacement.cpp
    ${TKET_TRANSFORM_DIR}/MeasurePass.cpp
    ${TKET_TRANSFORM_DIR}/ContextualReduction.cpp
    ${TKET_TRANSFORM_DIR}/ThreeQubitSquash.cpp
    ${TKET_TRANSFORM_DIR}/CZReordering.cpp

    # Routing
    ${TKET_ROUTING_DIR}/PlacementGraphClasses.cpp
    ${TKET_ROUTING_DIR}/Qubit_Placement.cpp
    ${TKET_ROUTING_DIR}/Swap_Analysis.cpp
    ${TKET_ROUTING_DIR}/Board_Analysis.cpp
    ${TKET_ROUTING_DIR}/Routing.cpp
    ${TKET_ROUTING_DIR}/Slice_Manipulation.cpp
    ${TKET_ROUTING_DIR}/subgraph_mapping.cpp
    ${TKET_ROUTING_DIR}/Placement.cpp
    ${TKET_ROUTING_DIR}/Verification.cpp

    # Architecture
    ${TKET_ARCHITECTURE_DIR}/ArchitectureGraphClasses.cpp
    ${TKET_ARCHITECTURE_DIR}/Architecture.cpp


    # Mapping
    ${TKET_MAPPING_DIR}/MappingFrontier.cpp
    ${TKET_MAPPING_DIR}/RoutingMethodCircuit.cpp
    ${TKET_MAPPING_DIR}/MappingManager.cpp
    ${TKET_MAPPING_DIR}/LexicographicalComparison.cpp
    ${TKET_MAPPING_DIR}/LexiRoute.cpp
    ${TKET_MAPPING_DIR}/RoutingMethodJson.cpp


    # Architecture Aware Synthesis
    ${TKET_AAS_DIR}/Path.cpp
    ${TKET_AAS_DIR}/SteinerTree.cpp
    ${TKET_AAS_DIR}/SteinerForest.cpp

    # Utils
    ${TKET_UTILS_DIR}/TketLog.cpp
    ${TKET_UTILS_DIR}/UnitID.cpp
    ${TKET_UTILS_DIR}/HelperFunctions.cpp
    ${TKET_UTILS_DIR}/MatrixAnalysis.cpp
    ${TKET_UTILS_DIR}/PauliStrings.cpp
    ${TKET_UTILS_DIR}/CosSinDecomposition.cpp
    ${TKET_UTILS_DIR}/Expression.cpp

    # Predicates
    ${TKET_PREDS_DIR}/Predicates.cpp
    ${TKET_PREDS_DIR}/CompilationUnit.cpp
    ${TKET_PREDS_DIR}/CompilerPass.cpp
    ${TKET_PREDS_DIR}/PassGenerators.cpp
    ${TKET_PREDS_DIR}/PassLibrary.cpp

    # PauliGraph
    ${TKET_PAULIGRAPH_DIR}/ConjugatePauliFunctions.cpp
    ${TKET_PAULIGRAPH_DIR}/PauliGraph.cpp

    # Converters
    ${TKET_CONVERTERS_DIR}/CliffTableauConverters.cpp
    ${TKET_CONVERTERS_DIR}/PauliGadget.cpp
    ${TKET_CONVERTERS_DIR}/PauliGraphConverters.cpp
    ${TKET_CONVERTERS_DIR}/Gauss.cpp
    ${TKET_CONVERTERS_DIR}/PhasePoly.cpp

    # Program
    ${TKET_PROGRAM_DIR}/Program_accessors.cpp
    ${TKET_PROGRAM_DIR}/Program_analysis.cpp
    ${TKET_PROGRAM_DIR}/Program_iteration.cpp
    ${TKET_PROGRAM_DIR}/Program_manipulation.cpp
    ${TKET_PROGRAM_DIR}/Program_units.cpp

    # MeasurementSetup
    ${TKET_MEASUREMENT_DIR}/MeasurementSetup.cpp
    ${TKET_MEASUREMENT_DIR}/MeasurementReduction.cpp

    # Characterisation
    ${TKET_CHARACTERISATION_DIR}/Cycles.cpp
    ${TKET_CHARACTERISATION_DIR}/FrameRandomisation.cpp
    ${TKET_CHARACTERISATION_DIR}/DeviceCharacterisation.cpp

    # ZX
    ${TKET_ZX_DIR}/ZXDConstructors.cpp
    ${TKET_ZX_DIR}/ZXDExpansions.cpp
    ${TKET_ZX_DIR}/ZXDGettersSetters.cpp
    ${TKET_ZX_DIR}/ZXDManipulation.cpp
    ${TKET_ZX_DIR}/ZXGenerator.cpp
    ${TKET_ZX_DIR}/ZXRWCombinators.cpp
    ${TKET_ZX_DIR}/ZXRWAxioms.cpp
    ${TKET_ZX_DIR}/ZXRWDecompositions.cpp
    ${TKET_ZX_DIR}/ZXRWGraphLikeForm.cpp
    ${TKET_ZX_DIR}/ZXRWGraphLikeSimplification.cpp
)

add_library(${TKET} SHARED ${TKET_SOURCES})
# ----- Location of header files ----------------------------------------------
target_include_directories(${TKET} PRIVATE ${TKET_SRC_DIR})
target_link_libraries(${TKET} ${CONAN_LIBS})
if(WIN32)
    target_link_libraries(${TKET} bcrypt)
ENDIF()
