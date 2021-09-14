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

# file to store all the files for tket
# new files should be added here

set(BUBBLE_ARCHITECTURE_DIR ${BUBBLE_SRC_DIR}/Architecture)
set(BUBBLE_AAS_DIR ${BUBBLE_SRC_DIR}/ArchAwareSynth)
set(BUBBLE_CIRCUIT_DIR ${BUBBLE_SRC_DIR}/Circuit)
set(BUBBLE_CLIFFORD_DIR ${BUBBLE_SRC_DIR}/Clifford)
set(BUBBLE_DIAGONALISATION_DIR ${BUBBLE_SRC_DIR}/Diagonalisation)
set(BUBBLE_GRAPH_DIR ${BUBBLE_SRC_DIR}/Graphs)
set(BUBBLE_UTILS_DIR ${BUBBLE_SRC_DIR}/Utils)
set(BUBBLE_OPTYPE_DIR ${BUBBLE_SRC_DIR}/OpType)
set(BUBBLE_OPS_DIR ${BUBBLE_SRC_DIR}/Ops)
set(BUBBLE_GATE_DIR ${BUBBLE_SRC_DIR}/Gate)
set(BUBBLE_SIMULATION_DIR ${BUBBLE_SRC_DIR}/Simulation)
set(BUBBLE_ROUTING_DIR ${BUBBLE_SRC_DIR}/Routing)
set(BUBBLE_TRANSFORM_DIR ${BUBBLE_SRC_DIR}/Transformations)
set(BUBBLE_CHARACTERISATION_DIR ${BUBBLE_SRC_DIR}/Characterisation)
set(BUBBLE_PREDS_DIR ${BUBBLE_SRC_DIR}/Predicates)
set(BUBBLE_PAULIGRAPH_DIR ${BUBBLE_SRC_DIR}/PauliGraph)
set(BUBBLE_CONVERTERS_DIR ${BUBBLE_SRC_DIR}/Converters)
set(BUBBLE_PROGRAM_DIR ${BUBBLE_SRC_DIR}/Program)
set(BUBBLE_MEASUREMENT_DIR ${BUBBLE_SRC_DIR}/MeasurementSetup)

set(TKET_SOURCES

    # OpType
    ${BUBBLE_OPTYPE_DIR}/OpDesc.cpp
    ${BUBBLE_OPTYPE_DIR}/OpTypeInfo.cpp
    ${BUBBLE_OPTYPE_DIR}/OpTypeFunctions.cpp
    ${BUBBLE_OPTYPE_DIR}/OpTypeJson.cpp

    # Ops
    ${BUBBLE_OPS_DIR}/Conditional.cpp
    ${BUBBLE_OPS_DIR}/ConjugatePauliFunctions.cpp
    ${BUBBLE_OPS_DIR}/FlowOp.cpp
    ${BUBBLE_OPS_DIR}/MetaOp.cpp
    ${BUBBLE_OPS_DIR}/Op.cpp
    ${BUBBLE_OPS_DIR}/ClassicalOps.cpp
    ${BUBBLE_OPS_DIR}/OpJsonFactory.cpp

    # Gate
    ${BUBBLE_GATE_DIR}/Gate.cpp
    ${BUBBLE_GATE_DIR}/GatePtr.cpp
    ${BUBBLE_GATE_DIR}/OpPtrFunctions.cpp
    ${BUBBLE_GATE_DIR}/Rotation.cpp
    ${BUBBLE_GATE_DIR}/SymTable.cpp
    ${BUBBLE_GATE_DIR}/GateUnitaryMatrix.cpp
    ${BUBBLE_GATE_DIR}/GateUnitaryMatrixComposites.cpp
    ${BUBBLE_GATE_DIR}/GateUnitaryMatrixError.cpp
    ${BUBBLE_GATE_DIR}/GateUnitaryMatrixFixedMatrices.cpp
    ${BUBBLE_GATE_DIR}/GateUnitaryMatrixPrimitives.cpp
    ${BUBBLE_GATE_DIR}/GateUnitaryMatrixUtils.cpp
    ${BUBBLE_GATE_DIR}/GateUnitaryMatrixVariableQubits.cpp
    ${BUBBLE_GATE_DIR}/GateUnitarySparseMatrix.cpp

    # Circuit
    ${BUBBLE_CIRCUIT_DIR}/Boxes.cpp
    ${BUBBLE_CIRCUIT_DIR}/Circuit.cpp
    ${BUBBLE_CIRCUIT_DIR}/CircuitJson.cpp
    ${BUBBLE_CIRCUIT_DIR}/CommandJson.cpp
    ${BUBBLE_CIRCUIT_DIR}/macro_manipulation.cpp
    ${BUBBLE_CIRCUIT_DIR}/basic_circ_manip.cpp
    ${BUBBLE_CIRCUIT_DIR}/latex_drawing.cpp
    ${BUBBLE_CIRCUIT_DIR}/macro_circ_info.cpp
    ${BUBBLE_CIRCUIT_DIR}/setters_and_getters.cpp
    ${BUBBLE_CIRCUIT_DIR}/CircUtils.cpp
    ${BUBBLE_CIRCUIT_DIR}/ThreeQubitConversion.cpp
    ${BUBBLE_CIRCUIT_DIR}/AssertionSynthesis.cpp
    ${BUBBLE_CIRCUIT_DIR}/CircPool.cpp
    ${BUBBLE_CIRCUIT_DIR}/DAGProperties.cpp
    ${BUBBLE_CIRCUIT_DIR}/OpJson.cpp

    # Simulation
    ${BUBBLE_SIMULATION_DIR}/BitOperations.cpp
    ${BUBBLE_SIMULATION_DIR}/CircuitSimulator.cpp
    ${BUBBLE_SIMULATION_DIR}/DecomposeCircuit.cpp
    ${BUBBLE_SIMULATION_DIR}/GateNode.cpp
    ${BUBBLE_SIMULATION_DIR}/GateNodesBuffer.cpp
    ${BUBBLE_SIMULATION_DIR}/PauliExpBoxUnitaryCalculator.cpp

    # Clifford
    ${BUBBLE_CLIFFORD_DIR}/CliffTableau.cpp

    # Diagonalisation
    ${BUBBLE_DIAGONALISATION_DIR}/DiagUtils.cpp
    ${BUBBLE_DIAGONALISATION_DIR}/Diagonalisation.cpp
    ${BUBBLE_DIAGONALISATION_DIR}/PauliPartition.cpp

    # Graph algorithms
    ${BUBBLE_GRAPH_DIR}/AdjacencyData.cpp
    ${BUBBLE_GRAPH_DIR}/BruteForceColouring.cpp
    ${BUBBLE_GRAPH_DIR}/ColouringPriority.cpp
    ${BUBBLE_GRAPH_DIR}/GraphColouring.cpp
    ${BUBBLE_GRAPH_DIR}/GraphRoutines.cpp
    ${BUBBLE_GRAPH_DIR}/LargeCliquesResult.cpp
    ${BUBBLE_GRAPH_DIR}/ArticulationPoints.cpp
    ${BUBBLE_GRAPH_DIR}/UIDConnectivity.cpp

    # Transformations
    ${BUBBLE_TRANSFORM_DIR}/Combinator.cpp
    ${BUBBLE_TRANSFORM_DIR}/Rebase.cpp
    ${BUBBLE_TRANSFORM_DIR}/BasicOptimisation.cpp
    ${BUBBLE_TRANSFORM_DIR}/PauliOptimisation.cpp
    ${BUBBLE_TRANSFORM_DIR}/CliffordOptimisation.cpp
    ${BUBBLE_TRANSFORM_DIR}/CliffordReductionPass.cpp
    ${BUBBLE_TRANSFORM_DIR}/OptimisationPass.cpp
    ${BUBBLE_TRANSFORM_DIR}/PhaseOptimisation.cpp
    ${BUBBLE_TRANSFORM_DIR}/ControlledGates.cpp
    ${BUBBLE_TRANSFORM_DIR}/Decomposition.cpp
    ${BUBBLE_TRANSFORM_DIR}/Replacement.cpp
    ${BUBBLE_TRANSFORM_DIR}/MeasurePass.cpp
    ${BUBBLE_TRANSFORM_DIR}/ContextualReduction.cpp
    ${BUBBLE_TRANSFORM_DIR}/ThreeQubitSquash.cpp

    # Routing
    ${BUBBLE_ROUTING_DIR}/Qubit_Placement.cpp
    ${BUBBLE_ROUTING_DIR}/Swap_Analysis.cpp
    ${BUBBLE_ROUTING_DIR}/Board_Analysis.cpp
    ${BUBBLE_ROUTING_DIR}/Routing.cpp
    ${BUBBLE_ROUTING_DIR}/Slice_Manipulation.cpp
    ${BUBBLE_ROUTING_DIR}/subgraph_mapping.cpp
    ${BUBBLE_ROUTING_DIR}/Placement.cpp
    ${BUBBLE_ROUTING_DIR}/Verification.cpp

    # Architecture
    ${BUBBLE_ARCHITECTURE_DIR}/Architectures.cpp

    # Architecture Aware Synthesis
    ${BUBBLE_AAS_DIR}/Path.cpp
    ${BUBBLE_AAS_DIR}/SteinerTree.cpp
    ${BUBBLE_AAS_DIR}/SteinerForest.cpp    

    # Utils
    ${BUBBLE_UTILS_DIR}/TketLog.cpp
    ${BUBBLE_UTILS_DIR}/UnitID.cpp
    ${BUBBLE_UTILS_DIR}/HelperFunctions.cpp
    ${BUBBLE_UTILS_DIR}/MatrixAnalysis.cpp
    ${BUBBLE_UTILS_DIR}/PauliStrings.cpp
    ${BUBBLE_UTILS_DIR}/CosSinDecomposition.cpp
    ${BUBBLE_UTILS_DIR}/Expression.cpp

    # Predicates
    ${BUBBLE_PREDS_DIR}/Predicates.cpp
    ${BUBBLE_PREDS_DIR}/CompilationUnit.cpp
    ${BUBBLE_PREDS_DIR}/CompilerPass.cpp
    ${BUBBLE_PREDS_DIR}/PassGenerators.cpp
    ${BUBBLE_PREDS_DIR}/PassLibrary.cpp

    # PauliGraph
    ${BUBBLE_PAULIGRAPH_DIR}/PauliGraph.cpp

    # Converters
    ${BUBBLE_CONVERTERS_DIR}/CliffTableauConverters.cpp
    ${BUBBLE_CONVERTERS_DIR}/PauliGadget.cpp
    ${BUBBLE_CONVERTERS_DIR}/PauliGraphConverters.cpp
    ${BUBBLE_CONVERTERS_DIR}/Gauss.cpp
    ${BUBBLE_CONVERTERS_DIR}/PhasePoly.cpp

    # Program
    ${BUBBLE_PROGRAM_DIR}/Program_accessors.cpp
    ${BUBBLE_PROGRAM_DIR}/Program_analysis.cpp
    ${BUBBLE_PROGRAM_DIR}/Program_iteration.cpp
    ${BUBBLE_PROGRAM_DIR}/Program_manipulation.cpp
    ${BUBBLE_PROGRAM_DIR}/Program_units.cpp

    # MeasurementSetup
    ${BUBBLE_MEASUREMENT_DIR}/MeasurementSetup.cpp
    ${BUBBLE_MEASUREMENT_DIR}/MeasurementReduction.cpp

    # Characterisation
    ${BUBBLE_CHARACTERISATION_DIR}/Cycles.cpp
    ${BUBBLE_CHARACTERISATION_DIR}/FrameRandomisation.cpp
    ${BUBBLE_CHARACTERISATION_DIR}/DeviceCharacterisation.cpp
)
