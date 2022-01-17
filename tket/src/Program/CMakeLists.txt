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

project(tket-program)

add_library(tket-program
    Program_accessors.cpp
    Program_analysis.cpp
    Program_iteration.cpp
    Program_manipulation.cpp
    Program_units.cpp)

target_include_directories(tket-program
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${TKET_Program_INCLUDE_DIR}
    ${TKET_Program_INCLUDE_DIR}/Program
    ${TKET_Utils_INCLUDE_DIR}
    ${TKET_OpType_INCLUDE_DIR}
    ${TKET_Ops_INCLUDE_DIR}
    ${TKET_Gate_INCLUDE_DIR}
    ${TKET_Circuit_INCLUDE_DIR})

target_link_libraries(tket-program
    PRIVATE
    tket-utils
    tket-ops
    tket-circuit
    ${CONAN_LIBS})
