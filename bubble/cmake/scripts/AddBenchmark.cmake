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

# Needed to build and run additional benchmarks without re-building the whole library
if(NOT TARGET benchmarks)
  add_custom_target(benchmarks-run)
endif()

function(add_benchmark targetname)

  cmake_parse_arguments(benchmark
    ""
    "WORKING_DIRECTORY"
    "INPUT_FILES;LIBRARIES;BENCHMARK;COMPILE_FLAGS;INCLUDES"
    ${ARGN}
    )

  # Set SOURCE from potential targetname extensions
  unset(SOURCE)
  if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${targetname}.cc")
    set(SOURCE ${targetname}.cc)
  elseif (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${targetname}.cpp")
    set(SOURCE ${targetname}.cpp)
  elseif("${benchmark_UNPARSED_ARGUMENTS}" STREQUAL "")
    message(FATAL_ERROR "No source given or found for ${targetname}.")
  endif()

  # Create benchmark executable, links and rename
  add_executable(benchmark_${targetname} ${SOURCE} ${benchmark_UNPARSED_ARGUMENTS})
  if (benchmark_INPUT_FILES)
    add_custom_command(
        TARGET benchmark_${targetname} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/input_files/${benchmark_INPUT_FILES}
                ${CMAKE_CURRENT_BINARY_DIR}/input_files/${benchmark_INPUT_FILES})
  endif()
  if (benchmark_LIBRARIES)
    target_link_libraries(benchmark_${targetname} ${benchmark_LIBRARIES})
  endif()
  if(benchmark_COMPILE_FLAGS)
    set_target_properties(benchmark_${targetname} PROPERTIES COMPILE_FLAGS ${benchmark_COMPILE_FLAGS})
  endif()
  if (benchmark_BENCHMARK AND ${benchmark_BENCHMARK} STREQUAL "googlebenchmark")
    target_link_libraries(benchmark_${targetname} benchmark)
  endif()
  if (benchmark_INCLUDES)
    target_include_directories(benchmark_${targetname} PRIVATE ${benchmark_INCLUDES})
  endif()
  set_target_properties(benchmark_${targetname} PROPERTIES OUTPUT_NAME ${targetname})

  # Build and run additional benchmarks either individually
  add_custom_target(benchmark_${targetname}-run COMMAND benchmark_${targetname}
    COMMENT "Running benchmark ${targetname}" )
  add_dependencies(benchmark_${targetname}-run benchmark_${targetname})
  # Or as a batch
  add_dependencies(benchmarks-run benchmark_${targetname}-run)
  
endfunction()
