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

# file to store all the files that serve as utils for the tket unit tests
# new files should be added here

set(TESTUTILS_SOURCES
    ${TKET_TESTS_DIR}/tests_main.cpp
    ${TKET_TESTS_DIR}/testutil.cpp
    ${TKET_TESTS_DIR}/CircuitsForTesting.cpp
    ${TKET_TESTS_DIR}/Graphs/EdgeSequence.cpp
    ${TKET_TESTS_DIR}/Graphs/EdgeSequenceColouringParameters.cpp
    ${TKET_TESTS_DIR}/Graphs/GraphTestingRoutines.cpp
    ${TKET_TESTS_DIR}/Graphs/RandomGraphGeneration.cpp
    ${TKET_TESTS_DIR}/Graphs/RandomPlanarGraphs.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/Data/FixedCompleteSolutions.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/Data/FixedSwapSequences.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/TableLookup/NeighboursFromEdges.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/TableLookup/PermutationTestUtils.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/TableLookup/SwapSequenceReductionTester.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/TestUtils/ArchitectureEdgesReimplementation.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/TestUtils/BestTsaTester.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/TestUtils/DebugFunctions.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/TestUtils/DecodedProblemData.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/TestUtils/FullTsaTesting.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/TestUtils/GetRandomSet.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/TestUtils/PartialTsaTesting.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/TestUtils/ProblemGeneration.cpp
    ${TKET_TESTS_DIR}/TokenSwapping/TestUtils/TestStatsStructs.cpp
    ${TKET_TESTS_DIR}/Gate/GatesData.cpp
    ${TKET_TESTS_DIR}/Simulation/ComparisonFunctions.cpp
)
