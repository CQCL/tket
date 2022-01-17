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

project(tket-utils)

add_library(tket-utils
    TketLog.cpp
    UnitID.cpp
    HelperFunctions.cpp
    MatrixAnalysis.cpp
    PauliStrings.cpp
    CosSinDecomposition.cpp
    Expression.cpp)

target_include_directories(tket-utils
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${TKET_Utils_INCLUDE_DIR}
    ${TKET_Utils_INCLUDE_DIR}/Utils)

target_link_libraries(tket-utils
    ${CONAN_LIBS})
