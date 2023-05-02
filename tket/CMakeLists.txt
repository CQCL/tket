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
project(tket CXX)

find_package(Boost CONFIG REQUIRED)
find_package(symengine CONFIG REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)
find_package(tklog CONFIG REQUIRED)
find_package(tkassert CONFIG REQUIRED)
find_package(tkrng CONFIG REQUIRED)
find_package(tktokenswap CONFIG REQUIRED)
find_package(tkwsm CONFIG REQUIRED)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
endif()

if(WIN32)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS yes)
elseif(APPLE)
    # set correct install_name
    set(CMAKE_INSTALL_NAME_DIR "@loader_path")
    set(CMAKE_BUILD_WITH_INSTALL_NAME_DIR ON)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

IF (WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /WX /EHsc")
ELSE()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -Wunreachable-code -Wunused")
ENDIF()
if(CMAKE_CXX_COMPILER_ID MATCHES "(Apple)?Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Winconsistent-missing-override -Wloop-analysis -Wno-deprecated-declarations")
    # remove -Wno-deprecated-declarations once https://github.com/boostorg/boost/issues/688 is resolved
endif()

add_library(tket)

set(PROFILE_COVERAGE no CACHE BOOL "Build library with profiling for test coverage")
IF (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    IF (PROFILE_COVERAGE)
        target_compile_options(tket PRIVATE -g --coverage)
        target_link_options(tket PUBLIC --coverage)
        # Bug in gcc 10: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=95353
        IF (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11)
            target_compile_options(tket PRIVATE -Wno-stringop-overflow)
        ENDIF()
    ENDIF()
ENDIF()

target_include_directories(tket PUBLIC include/tket)
target_link_libraries(tket PRIVATE boost::boost)
target_link_libraries(tket PRIVATE symengine::symengine)
target_link_libraries(tket PRIVATE Eigen3::Eigen)
target_link_libraries(tket PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(tket PRIVATE tklog::tklog)
target_link_libraries(tket PRIVATE tkassert::tkassert)
target_link_libraries(tket PRIVATE tkrng::tkrng)
target_link_libraries(tket PRIVATE tktokenswap::tktokenswap)
target_link_libraries(tket PRIVATE tkwsm::tkwsm)

target_sources(tket
    PRIVATE
        src/Utils/UnitID.cpp
        src/Utils/HelperFunctions.cpp
        src/Utils/MatrixAnalysis.cpp
        src/Utils/PauliStrings.cpp
        src/Utils/CosSinDecomposition.cpp
        src/Utils/Expression.cpp
        src/OpType/OpDesc.cpp
        src/OpType/OpTypeInfo.cpp
        src/OpType/OpTypeFunctions.cpp
        src/OpType/OpTypeJson.cpp
        src/Clifford/CliffTableau.cpp
        src/Clifford/ChoiMixTableau.cpp
        src/Clifford/SymplecticTableau.cpp
        src/Clifford/UnitaryTableau.cpp
        src/Ops/FlowOp.cpp
        src/Ops/MetaOp.cpp
        src/Ops/Op.cpp
        src/Ops/ClassicalOps.cpp
        src/Ops/OpJsonFactory.cpp
        src/Graphs/AdjacencyData.cpp
        src/Graphs/BruteForceColouring.cpp
        src/Graphs/ColouringPriority.cpp
        src/Graphs/GraphColouring.cpp
        src/Graphs/GraphRoutines.cpp
        src/Graphs/LargeCliquesResult.cpp
        src/Graphs/ArticulationPoints.cpp
        src/Gate/Gate.cpp
        src/Gate/GatePtr.cpp
        src/Gate/OpPtrFunctions.cpp
        src/Gate/Rotation.cpp
        src/Gate/SymTable.cpp
        src/Gate/GateUnitaryMatrix.cpp
        src/Gate/GateUnitaryMatrixComposites.cpp
        src/Gate/GateUnitaryMatrixError.cpp
        src/Gate/GateUnitaryMatrixFixedMatrices.cpp
        src/Gate/GateUnitaryMatrixPrimitives.cpp
        src/Gate/GateUnitaryMatrixUtils.cpp
        src/Gate/GateUnitaryMatrixVariableQubits.cpp
        src/Gate/GateUnitarySparseMatrix.cpp
        src/PauliGraph/ConjugatePauliFunctions.cpp
        src/PauliGraph/PauliGraph.cpp
        src/Circuit/Boxes.cpp
        src/Circuit/Circuit.cpp
        src/Circuit/CircuitJson.cpp
        src/Circuit/CommandJson.cpp
        src/Circuit/macro_manipulation.cpp
        src/Circuit/basic_circ_manip.cpp
        src/Circuit/latex_drawing.cpp
        src/Circuit/macro_circ_info.cpp
        src/Circuit/setters_and_getters.cpp
        src/Circuit/CircUtils.cpp
        src/Circuit/ThreeQubitConversion.cpp
        src/Circuit/AssertionSynthesis.cpp
        src/Circuit/CircPool.cpp
        src/Circuit/DAGProperties.cpp
        src/Circuit/OpJson.cpp
        src/Circuit/Conditional.cpp
        src/Circuit/ControlledGates.cpp
        src/Circuit/Multiplexor.cpp
        src/Circuit/StatePreparation.cpp
        src/Circuit/DiagonalBox.cpp
        src/Circuit/ToffoliBox.cpp
        src/Circuit/PauliExpBoxes.cpp
        src/Architecture/Architecture.cpp
        src/Architecture/ArchitectureGraphClasses.cpp
        src/Architecture/ArchitectureMapping.cpp
        src/Architecture/BestTsaWithArch.cpp
        src/Architecture/DistancesFromArchitecture.cpp
        src/Architecture/NeighboursFromArchitecture.cpp
        src/Architecture/SubgraphMonomorphisms.cpp
        src/Simulation/BitOperations.cpp
        src/Simulation/CircuitSimulator.cpp
        src/Simulation/DecomposeCircuit.cpp
        src/Simulation/GateNode.cpp
        src/Simulation/GateNodesBuffer.cpp
        src/Simulation/PauliExpBoxUnitaryCalculator.cpp
        src/Diagonalisation/DiagUtils.cpp
        src/Diagonalisation/Diagonalisation.cpp
        src/Diagonalisation/PauliPartition.cpp
        src/Characterisation/Cycles.cpp
        src/Characterisation/FrameRandomisation.cpp
        src/Characterisation/DeviceCharacterisation.cpp
        src/ZX/ZXDConstructors.cpp
        src/ZX/ZXDExpansions.cpp
        src/ZX/ZXDFormats.cpp
        src/ZX/ZXDGettersSetters.cpp
        src/ZX/ZXDManipulation.cpp
        src/ZX/ZXDSubdiagram.cpp
        src/ZX/ZXGenerator.cpp
        src/ZX/ZXRWCombinators.cpp
        src/ZX/ZXRWAxioms.cpp
        src/ZX/ZXRWDecompositions.cpp
        src/ZX/ZXRWGraphLikeForm.cpp
        src/ZX/ZXRWGraphLikeSimplification.cpp
        src/ZX/Flow.cpp
        src/ZX/MBQCRewrites.cpp
        src/ZX/ZXRWSequences.cpp
        src/Converters/CliffTableauConverters.cpp
        src/Converters/ChoiMixTableauConverters.cpp
        src/Converters/PauliGadget.cpp
        src/Converters/PauliGraphConverters.cpp
        src/Converters/Gauss.cpp
        src/Converters/PhasePoly.cpp
        src/Converters/UnitaryTableauBox.cpp
        src/Converters/UnitaryTableauConverters.cpp
        src/Converters/ZXConverters.cpp
        src/Placement/Frontier.cpp
        src/Placement/Placement.cpp
        src/Placement/GraphPlacement.cpp
        src/Placement/LinePlacement.cpp
        src/Placement/NoiseAwarePlacement.cpp
        src/Placement/PlacementGraphClasses.cpp
        src/Placement/NeighbourPlacements.cpp
        src/Placement/MonomorphismCalculation.cpp
        src/ArchAwareSynth/Path.cpp
        src/ArchAwareSynth/SteinerTree.cpp
        src/ArchAwareSynth/SteinerForest.cpp
        src/Mapping/AASRoute.cpp
        src/Mapping/AASLabelling.cpp
        src/Mapping/LexicographicalComparison.cpp
        src/Mapping/LexiRoute.cpp
        src/Mapping/LexiRouteRoutingMethod.cpp
        src/Mapping/LexiLabelling.cpp
        src/Mapping/MappingFrontier.cpp
        src/Mapping/MappingManager.cpp
        src/Mapping/MultiGateReorder.cpp
        src/Mapping/BoxDecomposition.cpp
        src/Mapping/RoutingMethodCircuit.cpp
        src/Mapping/RoutingMethodJson.cpp
        src/Mapping/Verification.cpp
        src/MeasurementSetup/MeasurementSetup.cpp
        src/MeasurementSetup/MeasurementReduction.cpp
        src/Transformations/Combinator.cpp
        src/Transformations/Rebase.cpp
        src/Transformations/BasicOptimisation.cpp
        src/Transformations/RedundancyRemoval.cpp
        src/Transformations/PauliOptimisation.cpp
        src/Transformations/CliffordOptimisation.cpp
        src/Transformations/CliffordReductionPass.cpp
        src/Transformations/OptimisationPass.cpp
        src/Transformations/PhaseOptimisation.cpp
        src/Transformations/Decomposition.cpp
        src/Transformations/Replacement.cpp
        src/Transformations/MeasurePass.cpp
        src/Transformations/ContextualReduction.cpp
        src/Transformations/ThreeQubitSquash.cpp
        src/Transformations/SingleQubitSquash.cpp
        src/Transformations/PhasedXFrontier.cpp
        src/Transformations/PQPSquash.cpp
        src/Transformations/RzPhasedXSquash.cpp
        src/Transformations/StandardSquash.cpp
        src/Predicates/Predicates.cpp
        src/Predicates/CompilationUnit.cpp
        src/Predicates/CompilerPass.cpp
        src/Predicates/PassGenerators.cpp
        src/Predicates/PassLibrary.cpp
    PUBLIC FILE_SET HEADERS
    BASE_DIRS ${PROJECT_SOURCE_DIR}/include
    FILES
        include/Utils/BiMapHeaders.hpp
        include/Utils/Constants.hpp
        include/Utils/CosSinDecomposition.hpp
        include/Utils/EigenConfig.hpp
        include/Utils/Expression.hpp
        include/Utils/GraphHeaders.hpp
        include/Utils/HelperFunctions.hpp
        include/Utils/Json.hpp
        include/Utils/MatrixAnalysis.hpp
        include/Utils/PauliStrings.hpp
        include/Utils/SequencedContainers.hpp
        include/Utils/Symbols.hpp
        include/Utils/UnitID.hpp
        include/OpType/EdgeType.hpp
        include/OpType/OpDesc.hpp
        include/OpType/OpTypeFunctions.hpp
        include/OpType/OpType.hpp
        include/OpType/OpTypeInfo.hpp
        include/Clifford/ChoiMixTableau.hpp
        include/Clifford/CliffTableau.hpp
        include/Clifford/SymplecticTableau.hpp
        include/Clifford/UnitaryTableau.hpp
        include/Ops/ClassicalOps.hpp
        include/Ops/FlowOp.hpp
        include/Ops/MetaOp.hpp
        include/Ops/Op.hpp
        include/Ops/OpJsonFactory.hpp
        include/Ops/OpPtr.hpp
        include/Graphs/AbstractGraph.hpp
        include/Graphs/AdjacencyData.hpp
        include/Graphs/ArticulationPoints.hpp
        include/Graphs/ArticulationPoints_impl.hpp
        include/Graphs/CompleteGraph.hpp
        include/Graphs/DirectedGraph.hpp
        include/Graphs/GraphColouring.hpp
        include/Graphs/GraphRoutines.hpp
        include/Graphs/LargeCliquesResult.hpp
        include/Graphs/TreeSearch.hpp
        include/Graphs/TreeSearch_impl.hpp
        include/Graphs/Utils.hpp
        include/Graphs/Utils_impl.hpp
        include/Gate/Gate.hpp
        include/Gate/GatePtr.hpp
        include/Gate/GateUnitaryMatrixError.hpp
        include/Gate/GateUnitaryMatrix.hpp
        include/Gate/GateUnitaryMatrixImplementations.hpp
        include/Gate/GateUnitaryMatrixUtils.hpp
        include/Gate/OpPtrFunctions.hpp
        include/Gate/Rotation.hpp
        include/Gate/SymTable.hpp
        include/PauliGraph/ConjugatePauliFunctions.hpp
        include/PauliGraph/PauliGraph.hpp
        include/Circuit/AssertionSynthesis.hpp
        include/Circuit/Boxes.hpp
        include/Circuit/CircPool.hpp
        include/Circuit/Circuit.hpp
        include/Circuit/CircUtils.hpp
        include/Circuit/ClassicalExpBox.hpp
        include/Circuit/Command.hpp
        include/Circuit/Conditional.hpp
        include/Circuit/DAGDefs.hpp
        include/Circuit/DiagonalBox.hpp
        include/Circuit/Multiplexor.hpp
        include/Circuit/StatePreparation.hpp
        include/Circuit/ThreeQubitConversion.hpp
        include/Circuit/ToffoliBox.hpp
        include/Circuit/PauliExpBoxes.hpp
        include/Architecture/Architecture.hpp
        include/Architecture/ArchitectureMapping.hpp
        include/Architecture/BestTsaWithArch.hpp
        include/Architecture/DistancesFromArchitecture.hpp
        include/Architecture/NeighboursFromArchitecture.hpp
        include/Architecture/SubgraphMonomorphisms.hpp
        include/Simulation/CircuitSimulator.hpp
        include/Simulation/PauliExpBoxUnitaryCalculator.hpp
        include/Diagonalisation/Diagonalisation.hpp
        include/Diagonalisation/DiagUtils.hpp
        include/Diagonalisation/PauliPartition.hpp
        include/Characterisation/Cycles.hpp
        include/Characterisation/DeviceCharacterisation.hpp
        include/Characterisation/ErrorTypes.hpp
        include/Characterisation/FrameRandomisation.hpp
        include/ZX/Flow.hpp
        include/ZX/Rewrite.hpp
        include/ZX/Types.hpp
        include/ZX/ZXDiagram.hpp
        include/ZX/ZXDiagramImpl.hpp
        include/ZX/ZXGenerator.hpp
        include/Converters/Converters.hpp
        include/Converters/Gauss.hpp
        include/Converters/PauliGadget.hpp
        include/Converters/PhasePoly.hpp
        include/Converters/UnitaryTableauBox.hpp
        include/Placement/NeighbourPlacements.hpp
        include/Placement/Placement.hpp
        include/Placement/QubitGraph.hpp
        include/ArchAwareSynth/Path.hpp
        include/ArchAwareSynth/SteinerForest.hpp
        include/ArchAwareSynth/SteinerTree.hpp
        include/Mapping/AASLabelling.hpp
        include/Mapping/AASRoute.hpp
        include/Mapping/BoxDecomposition.hpp
        include/Mapping/LexicographicalComparison.hpp
        include/Mapping/LexiLabelling.hpp
        include/Mapping/LexiRoute.hpp
        include/Mapping/LexiRouteRoutingMethod.hpp
        include/Mapping/MappingFrontier.hpp
        include/Mapping/MappingManager.hpp
        include/Mapping/MultiGateReorder.hpp
        include/Mapping/RoutingMethodCircuit.hpp
        include/Mapping/RoutingMethod.hpp
        include/Mapping/RoutingMethodJson.hpp
        include/Mapping/Verification.hpp
        include/MeasurementSetup/MeasurementReduction.hpp
        include/MeasurementSetup/MeasurementSetup.hpp
        include/Transformations/BasicOptimisation.hpp
        include/Transformations/CliffordOptimisation.hpp
        include/Transformations/CliffordReductionPass.hpp
        include/Transformations/Combinator.hpp
        include/Transformations/ContextualReduction.hpp
        include/Transformations/Decomposition.hpp
        include/Transformations/MeasurePass.hpp
        include/Transformations/OptimisationPass.hpp
        include/Transformations/PauliOptimisation.hpp
        include/Transformations/PhasedXFrontier.hpp
        include/Transformations/PhaseOptimisation.hpp
        include/Transformations/PQPSquash.hpp
        include/Transformations/Rebase.hpp
        include/Transformations/Replacement.hpp
        include/Transformations/RzPhasedXSquash.hpp
        include/Transformations/SingleQubitSquash.hpp
        include/Transformations/StandardSquash.hpp
        include/Transformations/ThreeQubitSquash.hpp
        include/Transformations/Transform.hpp
        include/Predicates/CompilationUnit.hpp
        include/Predicates/CompilerPass.hpp
        include/Predicates/PassGenerators.hpp
        include/Predicates/PassLibrary.hpp
        include/Predicates/Predicates.hpp
    )

install(TARGETS tket FILE_SET HEADERS)

if(MSVC)
    install(TARGETS tket
        RUNTIME DESTINATION bin
        RUNTIME DESTINATION lib
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)
endif()
