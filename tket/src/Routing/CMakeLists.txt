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

project(tket-Routing)

add_library(tket-Routing
    Qubit_Placement.cpp
    Swap_Analysis.cpp
    Board_Analysis.cpp
    Routing.cpp
    Slice_Manipulation.cpp
    subgraph_mapping.cpp
    Placement.cpp
    PlacementGraphClasses.cpp
    Verification.cpp)

target_include_directories(tket-Routing
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${TKET_Routing_INCLUDE_DIR}
    ${TKET_Routing_INCLUDE_DIR}/Routing
    ${TKET_Utils_INCLUDE_DIR}
    ${TKET_OpType_INCLUDE_DIR}
    ${TKET_Ops_INCLUDE_DIR}
    ${TKET_Gate_INCLUDE_DIR}
    ${TKET_Graphs_INCLUDE_DIR}
    ${TKET_Circuit_INCLUDE_DIR}
    ${TKET_Architecture_INCLUDE_DIR}
    ${TKET_Characterisation_INCLUDE_DIR})

target_link_libraries(tket-Routing
    PRIVATE
    tket-Utils
    tket-Gate
    tket-Graphs
    tket-Circuit
    tket-Architecture
    tket-Characterisation
    ${CONAN_LIBS})
