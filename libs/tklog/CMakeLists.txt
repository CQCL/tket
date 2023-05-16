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

list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_SOURCE_DIR}/cmake)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
endif()

if(WIN32)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS yes)
endif()

if(APPLE)
    set(CMAKE_INSTALL_NAME_DIR "@loader_path")
    set(CMAKE_BUILD_WITH_INSTALL_NAME_DIR ON)
endif()

add_library(tklog)

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

target_include_directories(tklog PUBLIC
    $<INSTALL_INTERFACE:include/tklog>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/tklog>)

target_sources(tklog
    PRIVATE src/TketLog.cpp
    PUBLIC FILE_SET HEADERS
    BASE_DIRS ${PROJECT_SOURCE_DIR}/include
    FILES include/tklog/TketLog.hpp)

include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/tklog)

install(TARGETS tklog
    EXPORT tklog-targets
    FILE_SET HEADERS)

if(MSVC)
    install(TARGETS tklog
        RUNTIME DESTINATION bin
        RUNTIME DESTINATION lib
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)
endif()

install(EXPORT tklog-targets
    FILE tklogTargets.cmake
    NAMESPACE tklog::
    DESTINATION ${INSTALL_CONFIGDIR}
)

include(CMakePackageConfigHelpers)

configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/tklogConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/tklogConfig.cmake
    INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/tklogConfig.cmake
    DESTINATION ${INSTALL_CONFIGDIR}
)

export(EXPORT tklog-targets
    FILE ${CMAKE_CURRENT_BINARY_DIR}/tklogTargets.cmake
    NAMESPACE tklog::)

export(PACKAGE tklog)
