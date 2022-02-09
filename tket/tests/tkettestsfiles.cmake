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

# file to store all the files for the tket unit tests
# new files should be added here

set(TEST_SOURCES
    # The ordering here is the order in which the tests are run
    # We should test simpler modules (e.g. Op, Circuit) before
    # the more complicated things that rely on them (e.g. Routing,
    # Transform) to help identify exactly where stuff breaks
    ${TKET_TESTS_DIR}/tests_main.cpp
    ${TKET_TESTS_DIR}/testutil.cpp
    ${TKET_TESTS_DIR}/CircuitsForTesting.cpp
    ${TKET_TESTS_DIR}/Utils/test_MatrixAnalysis.cpp
    ${TKET_TESTS_DIR}/Utils/test_CosSinDecomposition.cpp
    ${TKET_TESTS_DIR}/Graphs/EdgeSequence.cpp
    ${TKET_TESTS_DIR}/Graphs/EdgeSequenceColouringParameters.cpp
    ${TKET_TESTS_DIR}/Graphs/GraphTestingRoutines.cpp
    ${TKET_TESTS_DIR}/Graphs/RandomGraphGeneration.cpp
    ${TKET_TESTS_DIR}/Graphs/RandomPlanarGraphs.cpp
    ${TKET_TESTS_DIR}/Graphs/RNG.cpp
    ${TKET_TESTS_DIR}/Graphs/test_GraphColouring.cpp
    ${TKET_TESTS_DIR}/Graphs/test_GraphFindComponents.cpp
    ${TKET_TESTS_DIR}/Graphs/test_GraphFindMaxClique.cpp
    ${TKET_TESTS_DIR}/Graphs/test_RNG.cpp
    ${TKET_TESTS_DIR}/Graphs/test_GraphUtils.cpp
    ${TKET_TESTS_DIR}/Graphs/test_DirectedGraph.cpp
    ${TKET_TESTS_DIR}/Graphs/test_ArticulationPoints.cpp
    ${TKET_TESTS_DIR}/Graphs/test_TreeSearch.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/Data/FixedCompleteSolutions.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/Data/FixedSwapSequences.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/TableLookup/NeighboursFromEdges.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/TableLookup/PermutationTestUtils.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/TableLookup/SwapSequenceReductionTester.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/TableLookup/test_CanonicalRelabelling.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/TableLookup/test_ExactMappingLookup.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/TableLookup/test_FilteredSwapSequences.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/TableLookup/test_SwapSequenceReductions.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/TableLookup/test_SwapSequenceTable.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/TestUtils/ArchitectureEdgesReimplementation.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/TestUtils/BestTsaTester.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/TestUtils/DecodedProblemData.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/TestUtils/FullTsaTesting.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/TestUtils/PartialTsaTesting.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/TestUtils/ProblemGeneration.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/TestUtils/TestStatsStructs.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/TSAUtils/test_SwapFunctions.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/test_ArchitectureMappingEndToEnd.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/test_BestTsaFixedSwapSequences.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/test_DistancesFromArchitecture.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/test_FullTsa.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/test_main_entry_functions.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/test_RiverFlowPathFinder.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/test_SwapList.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/test_SwapListOptimiser.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/test_VariousPartialTsa.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/test_VectorListHybrid.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/test_VectorListHybridSkeleton.cpp
    ${TKET_TESTS_DIR}/test_PauliString.cpp
    ${TKET_TESTS_DIR}/Ops/test_ClassicalOps.cpp
    ${TKET_TESTS_DIR}/Ops/test_Expression.cpp
    ${TKET_TESTS_DIR}/Ops/test_Ops.cpp
    ${TKET_TESTS_DIR}/Gate/GatesData.cpp
    ${TKET_TESTS_DIR}/Gate/test_GateUnitaryMatrix.cpp
    ${TKET_TESTS_DIR}/Simulation/ComparisonFunctions.cpp
    ${TKET_TESTS_DIR}/Simulation/test_CircuitSimulator.cpp
    ${TKET_TESTS_DIR}/Simulation/test_PauliExpBoxUnitaryCalculator.cpp
    ${TKET_TESTS_DIR}/test_Utils.cpp
    ${TKET_TESTS_DIR}/Circuit/test_Boxes.cpp
    ${TKET_TESTS_DIR}/Circuit/test_Circ.cpp
    ${TKET_TESTS_DIR}/Circuit/test_Symbolic.cpp
    ${TKET_TESTS_DIR}/Circuit/test_ThreeQubitConversion.cpp
    ${TKET_TESTS_DIR}/test_Program.cpp
    ${TKET_TESTS_DIR}/test_CliffTableau.cpp
    ${TKET_TESTS_DIR}/test_UnitaryTableau.cpp
    ${TKET_TESTS_DIR}/test_PhasePolynomials.cpp
    ${TKET_TESTS_DIR}/test_PauliGraph.cpp
    ${TKET_TESTS_DIR}/test_Architectures.cpp
    ${TKET_TESTS_DIR}/test_Placement.cpp
    ${TKET_TESTS_DIR}/test_MappingFrontier.cpp
    ${TKET_TESTS_DIR}/test_RoutingMethod.cpp
    ${TKET_TESTS_DIR}/test_MappingManager.cpp
    ${TKET_TESTS_DIR}/test_LexicographicalComparison.cpp
    ${TKET_TESTS_DIR}/test_LexiRoute.cpp
    ${TKET_TESTS_DIR}/test_AASRoute.cpp
    ${TKET_TESTS_DIR}/test_MultiGateReorder.cpp
    ${TKET_TESTS_DIR}/test_DeviceCharacterisation.cpp
    ${TKET_TESTS_DIR}/test_Clifford.cpp
    ${TKET_TESTS_DIR}/test_MeasurementSetup.cpp
    ${TKET_TESTS_DIR}/test_Partition.cpp
    ${TKET_TESTS_DIR}/test_MeasurementReduction.cpp
    ${TKET_TESTS_DIR}/test_PhaseGadget.cpp
    ${TKET_TESTS_DIR}/test_Rebase.cpp
    ${TKET_TESTS_DIR}/test_Synthesis.cpp
    ${TKET_TESTS_DIR}/test_TwoQubitCanonical.cpp
    ${TKET_TESTS_DIR}/test_ControlDecomp.cpp
    ${TKET_TESTS_DIR}/test_Combinators.cpp
    ${TKET_TESTS_DIR}/test_Predicates.cpp
    ${TKET_TESTS_DIR}/test_CompilerPass.cpp
    ${TKET_TESTS_DIR}/test_ContextOpt.cpp
    ${TKET_TESTS_DIR}/test_FrameRandomisation.cpp
    ${TKET_TESTS_DIR}/test_Assertion.cpp
    ${TKET_TESTS_DIR}/test_json.cpp
    ${TKET_TESTS_DIR}/test_Path.cpp
    ${TKET_TESTS_DIR}/test_SteinerTree.cpp
    ${TKET_TESTS_DIR}/test_SteinerForest.cpp
    ${TKET_TESTS_DIR}/ZX/test_ZXDiagram.cpp
    ${TKET_TESTS_DIR}/ZX/test_ZXAxioms.cpp
    ${TKET_TESTS_DIR}/ZX/test_ZXSimp.cpp
)
