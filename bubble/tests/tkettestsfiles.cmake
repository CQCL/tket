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

set(TEST_SOURCES
    ${BUBBLE_TESTS_DIR}/tests_main.cpp
    # Single file containing all tests to guarantee test order
    ${BUBBLE_TESTS_DIR}/tkettestseq.cpp
    # Additional utilities to support testing
    ${BUBBLE_TESTS_DIR}/testutil.cpp
    ${BUBBLE_TESTS_DIR}/CircuitsForTesting.cpp
    ${BUBBLE_TESTS_DIR}/Graphs/EdgeSequence.cpp
    ${BUBBLE_TESTS_DIR}/Graphs/EdgeSequenceColouringParameters.cpp
    ${BUBBLE_TESTS_DIR}/Graphs/GraphTestingRoutines.cpp
    ${BUBBLE_TESTS_DIR}/Graphs/RandomGraphGeneration.cpp
    ${BUBBLE_TESTS_DIR}/Graphs/RandomPlanarGraphs.cpp
    ${BUBBLE_TESTS_DIR}/Graphs/RNG.cpp
    ${BUBBLE_TESTS_DIR}/Gate/GatesData.cpp
    ${BUBBLE_TESTS_DIR}/Simulation/ComparisonFunctions.cpp
)
