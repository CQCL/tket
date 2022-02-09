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

cmake_minimum_required(VERSION 3.16.0)

set (CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

cmake_minimum_required(VERSION 3.16.0)

project(PackageTest CXX)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

IF (WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc")
ENDIF()

include_directories(${CONAN_INCLUDE_DIRS_TKET})

add_executable(test test.cpp)

target_link_libraries(test PRIVATE
    tket-Circuit
    tket-Gate
    tket-OpType
    tket-Transformations
    tket-Utils)

target_link_libraries(test PRIVATE
    ${CONAN_LIBS_FMT} ${CONAN_LIBS_SPDLOG} ${CONAN_LIBS_SYMENGINE})
