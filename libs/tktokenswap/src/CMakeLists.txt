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

cmake_minimum_required(VERSION 3.16.0)
project(tktokenswap CXX)

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

add_library(tktokenswap
    BestFullTsa.cpp
    CyclesCandidateManager.cpp
    CyclesGrowthManager.cpp
    CyclesPartialTsa.cpp
    CyclicShiftCostEstimate.cpp
    DistancesInterface.cpp
    DynamicTokenTracker.cpp
    HybridTsa.cpp
    NeighboursInterface.cpp
    PartialTsaInterface.cpp
    RiverFlowPathFinder.cpp
    SwapListOptimiser.cpp
    TrivialTSA.cpp
    VectorListHybridSkeleton.cpp
    TSAUtils/DistanceFunctions.cpp
    TSAUtils/SwapFunctions.cpp
    TSAUtils/VertexMappingFunctions.cpp
    TSAUtils/VertexSwapResult.cpp
    TableLookup/CanonicalRelabelling.cpp
    TableLookup/ExactMappingLookup.cpp
    TableLookup/FilteredSwapSequences.cpp
    TableLookup/PartialMappingLookup.cpp
    TableLookup/SwapConversion.cpp
    TableLookup/SwapListSegmentOptimiser.cpp
    TableLookup/SwapListTableOptimiser.cpp
    TableLookup/SwapSequenceTable.cpp
    TableLookup/VertexMapResizing.cpp
    )
target_include_directories(tktokenswap PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(tktokenswap PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include/tktokenswap)
target_include_directories(tktokenswap INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(tktokenswap PRIVATE
    ${CONAN_LIBS_TKLOG} ${CONAN_LIBS_TKASSERT} ${CONAN_LIBS_TKRNG})

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
