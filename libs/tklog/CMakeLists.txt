# Copyright 2019-2023 Cambridge Quantum Computing
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

cmake_minimum_required(VERSION 3.23)
project(tklog CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
endif()

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

if(WIN32)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS yes)
endif()

if(APPLE)
    # set correct install_name
    set(CMAKE_INSTALL_NAME_DIR "@loader_path")
    set(CMAKE_BUILD_WITH_INSTALL_NAME_DIR ON)
endif()

add_library(tklog TketLog.cpp)
target_include_directories(tklog PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include/tklog)
target_include_directories(tklog INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)


if(MSVC)
  target_compile_options(tklog PRIVATE /W4 /WX)
else()
  target_compile_options(tklog PRIVATE -Wall -Wextra -Wpedantic -Werror)
endif()

set(PROFILE_COVERAGE no CACHE BOOL "Build library with profiling for test coverage")
IF (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    IF (PROFILE_COVERAGE)
        target_compile_options(tklog PRIVATE -g --coverage)
        target_link_options(tklog PUBLIC --coverage)
        # Bug in gcc 10: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=95353
        IF (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11)
            target_compile_options(tklog PRIVATE -Wno-stringop-overflow)
        ENDIF()
    ENDIF()
ENDIF()
