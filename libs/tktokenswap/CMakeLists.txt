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
project(tktokenswap CXX)

list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_SOURCE_DIR}/cmake)

find_package(tklog CONFIG REQUIRED)
find_package(tkassert CONFIG REQUIRED)
find_package(tkrng CONFIG REQUIRED)
find_package(Boost CONFIG REQUIRED)

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

add_library(tktokenswap)

if(MSVC)
  target_compile_options(tktokenswap PRIVATE /W4 /WX /wd4267)
else()
  target_compile_options(tktokenswap PRIVATE -Wall -Wextra -Wpedantic -Werror)
endif()

set(PROFILE_COVERAGE no CACHE BOOL "Build library with profiling for test coverage")
IF (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    IF (PROFILE_COVERAGE)
        target_compile_options(tktokenswap PRIVATE -g --coverage)
        target_link_options(tktokenswap PUBLIC --coverage)
        # Bug in gcc 10: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=95353
        IF (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11)
            target_compile_options(tktokenswap PRIVATE -Wno-stringop-overflow)
        ENDIF()
    ENDIF()
ENDIF()

target_include_directories(tktokenswap PUBLIC
    $<INSTALL_INTERFACE:include/tktokenswap>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/tktokenswap>)
target_link_libraries(tktokenswap PRIVATE tklog::tklog)
target_link_libraries(tktokenswap PUBLIC tkassert::tkassert)
target_link_libraries(tktokenswap PRIVATE tkrng::tkrng)
target_link_libraries(tktokenswap PRIVATE Boost::headers)

target_sources(tktokenswap
    PRIVATE
        src/BestFullTsa.cpp
        src/CyclesCandidateManager.cpp
        src/CyclesGrowthManager.cpp
        src/CyclesPartialTsa.cpp
        src/CyclicShiftCostEstimate.cpp
        src/DistancesInterface.cpp
        src/DynamicTokenTracker.cpp
        src/HybridTsa.cpp
        src/NeighboursInterface.cpp
        src/PartialTsaInterface.cpp
        src/RiverFlowPathFinder.cpp
        src/SwapListOptimiser.cpp
        src/TrivialTSA.cpp
        src/VectorListHybridSkeleton.cpp
        src/TSAUtils/DistanceFunctions.cpp
        src/TSAUtils/SwapFunctions.cpp
        src/TSAUtils/VertexMappingFunctions.cpp
        src/TSAUtils/VertexSwapResult.cpp
        src/TableLookup/CanonicalRelabelling.cpp
        src/TableLookup/ExactMappingLookup.cpp
        src/TableLookup/FilteredSwapSequences.cpp
        src/TableLookup/PartialMappingLookup.cpp
        src/TableLookup/SwapConversion.cpp
        src/TableLookup/SwapListSegmentOptimiser.cpp
        src/TableLookup/SwapListTableOptimiser.cpp
        src/TableLookup/SwapSequenceTable.cpp
        src/TableLookup/VertexMapResizing.cpp
    PUBLIC FILE_SET HEADERS
    BASE_DIRS ${PROJECT_SOURCE_DIR}/include
    FILES
        include/tktokenswap/BestFullTsa.hpp
        include/tktokenswap/CanonicalRelabelling.hpp
        include/tktokenswap/CyclesCandidateManager.hpp
        include/tktokenswap/CyclesGrowthManager.hpp
        include/tktokenswap/CyclesPartialTsa.hpp
        include/tktokenswap/DistanceFunctions.hpp
        include/tktokenswap/DistancesInterface.hpp
        include/tktokenswap/DynamicTokenTracker.hpp
        include/tktokenswap/ExactMappingLookup.hpp
        include/tktokenswap/FilteredSwapSequences.hpp
        include/tktokenswap/GeneralFunctions.hpp
        include/tktokenswap/HybridTsa.hpp
        include/tktokenswap/NeighboursInterface.hpp
        include/tktokenswap/PartialMappingLookup.hpp
        include/tktokenswap/PartialTsaInterface.hpp
        include/tktokenswap/RiverFlowPathFinder.hpp
        include/tktokenswap/SwapConversion.hpp
        include/tktokenswap/SwapFunctions.hpp
        include/tktokenswap/SwapListOptimiser.hpp
        include/tktokenswap/SwapListSegmentOptimiser.hpp
        include/tktokenswap/SwapListTableOptimiser.hpp
        include/tktokenswap/SwapSequenceTable.hpp
        include/tktokenswap/TrivialTSA.hpp
        include/tktokenswap/VectorListHybrid.hpp
        include/tktokenswap/VectorListHybridSkeleton.hpp
        include/tktokenswap/VertexMappingFunctions.hpp
        include/tktokenswap/VertexMapResizing.hpp
        include/tktokenswap/VertexSwapResult.hpp
        )

include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/tktokenswap)

install(TARGETS tktokenswap
    EXPORT tktokenswap-targets
    FILE_SET HEADERS)

if(MSVC)
    install(TARGETS tktokenswap
        RUNTIME DESTINATION bin
        RUNTIME DESTINATION lib
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)
endif()

install(EXPORT tktokenswap-targets
    FILE tktokenswapTargets.cmake
    NAMESPACE tktokenswap::
    DESTINATION ${INSTALL_CONFIGDIR}
)

include(CMakePackageConfigHelpers)

configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/tktokenswapConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/tktokenswapConfig.cmake
    INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/tktokenswapConfig.cmake
    DESTINATION ${INSTALL_CONFIGDIR}
)

export(EXPORT tktokenswap-targets
    FILE ${CMAKE_CURRENT_BINARY_DIR}/tktokenswapTargets.cmake
    NAMESPACE tktokenswap::)

export(PACKAGE tktokenswap)
