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
project(tket CXX)

list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_SOURCE_DIR}/cmake)

#find_package(Boost CONFIG REQUIRED)
#find_package(SymEngine CONFIG REQUIRED)
#find_package(Eigen3 CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)
#find_package(tklog CONFIG REQUIRED)
#find_package(tkassert CONFIG REQUIRED)
#find_package(tkrng CONFIG REQUIRED)
#find_package(tktokenswap CONFIG REQUIRED)
#find_package(tkwsm CONFIG REQUIRED)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_program(ccache_exe ccache)
if(ccache_exe)
    message("ccache_exe ${ccache_exe}")
    file(REAL_PATH ${ccache_exe} ccache_real_exe)
    message("real ccache_exe ${ccache_real_exe}")
    get_filename_component(ccache_dir ${ccache_real_exe})
    message("ccache_dir ${ccache_dir}")
    file(COPY_FILE
            ${ccache_real_exe} ${ccache_dir}/cl.exe
            ONLY_IF_DIFFERENT)

    set(CMAKE_VS_GLOBALS
            "CLToolExe=cl.exe"
            "CLToolPath=${ccache_dir}"
            "TrackFileAccess=false"
            "UseMultiToolTask=true"
            "DebugInformationFormat=OldStyle"
            )
endif()
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    message("CCACHE is found")
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
else(CCACHE_PROGRAM)
    message("CCACHE is NOT found")
endif()

if(WIN32)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS yes)
elseif(APPLE)
    # set correct install_name
    set(CMAKE_INSTALL_NAME_DIR "@loader_path")
    set(CMAKE_BUILD_WITH_INSTALL_NAME_DIR ON)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

IF (WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /WX /EHsc")
ELSE()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -Wunreachable-code -Wunused")
ENDIF()
if(CMAKE_CXX_COMPILER_ID MATCHES "(Apple)?Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Winconsistent-missing-override -Wloop-analysis -Wno-deprecated-declarations")
    # remove -Wno-deprecated-declarations once https://github.com/boostorg/boost/issues/688 is resolved
endif()

add_library(tket)

set(PROFILE_COVERAGE no CACHE BOOL "Build library with profiling for test coverage")
IF (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    IF (PROFILE_COVERAGE)
        target_compile_options(tket PRIVATE -g --coverage)
        target_link_options(tket PUBLIC --coverage)
        # Bug in gcc 10: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=95353
        IF (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11)
            target_compile_options(tket PRIVATE -Wno-stringop-overflow)
        ENDIF()
    ENDIF()
ENDIF()

#if (NOT TARGET symengine::symengine)
#    add_library(symengine::symengine ALIAS symengine)
#endif()

target_include_directories(tket PUBLIC
    $<INSTALL_INTERFACE:include/tket>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/tket>)
#target_link_libraries(tket PRIVATE Boost::headers)
#target_link_libraries(tket PRIVATE symengine::symengine)
#target_link_libraries(tket PRIVATE Eigen3::Eigen)
target_link_libraries(tket PRIVATE nlohmann_json::nlohmann_json)
#target_link_libraries(tket PRIVATE tklog::tklog)
#target_link_libraries(tket PRIVATE tkassert::tkassert)
#target_link_libraries(tket PRIVATE tkrng::tkrng)
#target_link_libraries(tket PRIVATE tktokenswap::tktokenswap)
#target_link_libraries(tket PRIVATE tkwsm::tkwsm)

target_sources(tket
    PRIVATE
        src/Utils/UnitID.cpp
      #  src/Utils/CosSinDecomposition.cpp
      #  src/Utils/HelperFunctions.cpp
      #  src/Utils/MatrixAnalysis.cpp
      #  src/Utils/PauliStrings.cpp
    PUBLIC FILE_SET HEADERS
    BASE_DIRS ${PROJECT_SOURCE_DIR}/include
    FILES
        include/tket/Utils/UnitID.hpp
      #  include/tket/Utils/CosSinDecomposition.hpp
      #  include/tket/Utils/HelperFunctions.hpp
      #  include/tket/Utils/MatrixAnalysis.hpp
      #  include/tket/Utils/PauliStrings.hpp
      #  include/tket/Utils/BiMapHeaders.hpp
      #  include/tket/Utils/Constants.hpp
      #  include/tket/Utils/EigenConfig.hpp
      #  include/tket/Utils/GraphHeaders.hpp
      #  include/tket/Utils/MatrixAnalysis.hpp
      #  include/tket/Utils/SequencedContainers.hpp
      #  include/tket/Utils/Symbols.hpp
        include/tket/Utils/Json.hpp
    )

#if (BUILD_TKET_TEST)
#    add_subdirectory(test)
#endif()
#if (BUILD_TKET_PROPTEST)
#    add_subdirectory(proptest)
#endif()

include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/tket)

install(TARGETS tket
    EXPORT tket-targets
    FILE_SET HEADERS)

if(MSVC)
    install(TARGETS tket
        RUNTIME DESTINATION bin
        RUNTIME DESTINATION lib
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)
endif()

install(EXPORT tket-targets
    FILE tketTargets.cmake
    NAMESPACE tket::
    DESTINATION ${INSTALL_CONFIGDIR}
)

include(CMakePackageConfigHelpers)

configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/tketConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/tketConfig.cmake
    INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/tketConfig.cmake
    DESTINATION ${INSTALL_CONFIGDIR}
)

export(EXPORT tket-targets
    FILE ${CMAKE_CURRENT_BINARY_DIR}/tketTargets.cmake
    NAMESPACE tket::)

export(PACKAGE tket)
