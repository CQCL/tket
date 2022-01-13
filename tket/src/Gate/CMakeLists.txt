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

project(tket-gate)

add_library(tket-gate
    Gate.cpp
    GatePtr.cpp
    OpPtrFunctions.cpp
    Rotation.cpp
    SymTable.cpp
    GateUnitaryMatrix.cpp
    GateUnitaryMatrixComposites.cpp
    GateUnitaryMatrixError.cpp
    GateUnitaryMatrixFixedMatrices.cpp
    GateUnitaryMatrixPrimitives.cpp
    GateUnitaryMatrixUtils.cpp
    GateUnitaryMatrixVariableQubits.cpp
    GateUnitarySparseMatrix.cpp)

target_include_directories(tket-gate PRIVATE ${TKET_SRC_DIR})

target_link_libraries(tket-gate
    PRIVATE
    tket-utils
    tket-optype
    tket-ops
    ${CONAN_LIBS})
