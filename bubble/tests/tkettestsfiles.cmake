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

# file to store all the files for the tket unit tests
# new files should be added here

set(TEST_SOURCES
    # The ordering here is the order in which the tests are run
    # We should test simpler modules (e.g. Op, Circuit) before
    # the more complicated things that rely on them (e.g. Routing,
    # Transform) to help identify exactly where stuff breaks
    ${BUBBLE_TESTS_DIR}/tests_main.cpp
    ${BUBBLE_TESTS_DIR}/testutil.cpp
    ${BUBBLE_TESTS_DIR}/CircuitsForTesting.cpp
    ${BUBBLE_TESTS_DIR}/Utils/test_MatrixAnalysis.cpp
    ${BUBBLE_TESTS_DIR}/Utils/test_CosSinDecomposition.cpp
    ${BUBBLE_TESTS_DIR}/Graphs/EdgeSequence.cpp
    ${BUBBLE_TESTS_DIR}/Graphs/EdgeSequenceColouringParameters.cpp
    ${BUBBLE_TESTS_DIR}/Graphs/GraphTestingRoutines.cpp
    ${BUBBLE_TESTS_DIR}/Graphs/RandomGraphGeneration.cpp
    ${BUBBLE_TESTS_DIR}/Graphs/RandomPlanarGraphs.cpp
    ${BUBBLE_TESTS_DIR}/Graphs/RNG.cpp
    ${BUBBLE_TESTS_DIR}/Graphs/test_GraphColouring.cpp
    ${BUBBLE_TESTS_DIR}/Graphs/test_GraphFindComponents.cpp
    ${BUBBLE_TESTS_DIR}/Graphs/test_GraphFindMaxClique.cpp
    ${BUBBLE_TESTS_DIR}/Graphs/test_RNG.cpp
    ${BUBBLE_TESTS_DIR}/Graphs/test_GraphUtils.cpp
    ${BUBBLE_TESTS_DIR}/Graphs/test_UIDConnectivity.cpp
    ${BUBBLE_TESTS_DIR}/Graphs/test_ArticulationPoints.cpp
    ${BUBBLE_TESTS_DIR}/Graphs/test_TreeSearch.cpp
    ${BUBBLE_TESTS_DIR}/test_PauliString.cpp
    ${BUBBLE_TESTS_DIR}/Ops/test_ClassicalOps.cpp
    ${BUBBLE_TESTS_DIR}/Ops/test_Expression.cpp
    ${BUBBLE_TESTS_DIR}/Ops/test_Ops.cpp
    ${BUBBLE_TESTS_DIR}/Gate/GatesData.cpp
    ${BUBBLE_TESTS_DIR}/Gate/test_GateUnitaryMatrix.cpp
    ${BUBBLE_TESTS_DIR}/Simulation/ComparisonFunctions.cpp
    ${BUBBLE_TESTS_DIR}/Simulation/test_CircuitSimulator.cpp
    ${BUBBLE_TESTS_DIR}/Simulation/test_PauliExpBoxUnitaryCalculator.cpp
    ${BUBBLE_TESTS_DIR}/test_Utils.cpp
    ${BUBBLE_TESTS_DIR}/Circuit/test_Boxes.cpp
    ${BUBBLE_TESTS_DIR}/Circuit/test_Circ.cpp
    ${BUBBLE_TESTS_DIR}/Circuit/test_Symbolic.cpp
    ${BUBBLE_TESTS_DIR}/Circuit/test_ThreeQubitConversion.cpp
    ${BUBBLE_TESTS_DIR}/test_Program.cpp
    ${BUBBLE_TESTS_DIR}/test_CliffTableau.cpp
    ${BUBBLE_TESTS_DIR}/test_PhasePolynomials.cpp
    ${BUBBLE_TESTS_DIR}/test_PauliGraph.cpp
    ${BUBBLE_TESTS_DIR}/test_Architectures.cpp
    ${BUBBLE_TESTS_DIR}/test_Placement.cpp
    ${BUBBLE_TESTS_DIR}/test_Routing.cpp
    ${BUBBLE_TESTS_DIR}/test_DeviceCharacterisation.cpp
    ${BUBBLE_TESTS_DIR}/test_Clifford.cpp
    ${BUBBLE_TESTS_DIR}/test_MeasurementSetup.cpp
    ${BUBBLE_TESTS_DIR}/test_Partition.cpp
    ${BUBBLE_TESTS_DIR}/test_MeasurementReduction.cpp
    ${BUBBLE_TESTS_DIR}/test_PhaseGadget.cpp
    ${BUBBLE_TESTS_DIR}/test_Rebase.cpp
    ${BUBBLE_TESTS_DIR}/test_Synthesis.cpp
    ${BUBBLE_TESTS_DIR}/test_TwoQubitCanonical.cpp
    ${BUBBLE_TESTS_DIR}/test_ControlDecomp.cpp
    ${BUBBLE_TESTS_DIR}/test_Combinators.cpp
    ${BUBBLE_TESTS_DIR}/test_Predicates.cpp
    ${BUBBLE_TESTS_DIR}/test_CompilerPass.cpp
    ${BUBBLE_TESTS_DIR}/test_ContextOpt.cpp
    ${BUBBLE_TESTS_DIR}/test_FrameRandomisation.cpp
    ${BUBBLE_TESTS_DIR}/test_Assertion.cpp
    ${BUBBLE_TESTS_DIR}/test_json.cpp
    ${BUBBLE_TESTS_DIR}/test_Path.cpp
    ${BUBBLE_TESTS_DIR}/test_SteinerTree.cpp
    ${BUBBLE_TESTS_DIR}/test_SteinerForest.cpp
)
