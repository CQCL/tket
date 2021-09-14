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

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

cmake_minimum_required(VERSION 3.16.0)

set(BUBBLE_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

project(tket-compiler)

include(tketCMake.cmake)

include(tketsrcfiles.cmake)

add_library(${TKET} SHARED ${TKET_SOURCES})
# ----- Location of header files ----------------------------------------------
target_include_directories(${TKET} PRIVATE ${BUBBLE_SRC_DIR})
target_link_libraries(${TKET} ${CONAN_LIBS})
if(WIN32)
    target_link_libraries(${TKET} bcrypt)
ENDIF()
