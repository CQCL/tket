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

# INCLUDES are PRIVATE
add_benchmark(circuit
  LIBRARIES
    tket
  BENCHMARK			# Already adds benchmark specific includes
    googlebenchmark
  INCLUDES
    ${TKET_SRC_DIR} ${TKET_INCLUDE_DIR})
# INCLUDES are PRIVATE
add_benchmark(circuit_random
  INPUT_FILES
    circuit_random_nb_qubits=20_nb_layers=200_example.tkc
  LIBRARIES
    tket
  BENCHMARK			# Already adds benchmark specific includes
    googlebenchmark
  INCLUDES
    ${TKET_SRC_DIR} ${TKET_INCLUDE_DIR})
