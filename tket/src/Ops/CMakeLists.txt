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

project(tket-Ops)

add_library(tket-Ops
    FlowOp.cpp
    MetaOp.cpp
    Op.cpp
    ClassicalOps.cpp
    OpJsonFactory.cpp)

target_include_directories(tket-Ops
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${TKET_Ops_INCLUDE_DIR}
    ${TKET_Ops_INCLUDE_DIR}/Ops
    ${TKET_Utils_INCLUDE_DIR}
    ${TKET_OpType_INCLUDE_DIR})

target_link_libraries(tket-Ops
    PRIVATE
    tket-Utils
    tket-OpType
    ${CONAN_LIBS})
