// Copyright 2019-2022 Cambridge Quantum Computing
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once

#include "Circuit/Circuit.hpp"

namespace tket {

/** These circuits are not necessarily fundamentally significant,
 *  but they do appear in multiple places in the tests.
 */
struct CircuitsForTesting {
  Circuit uccsd;

  /** This is often accompanied by the comment
   *  "add some arbitrary rotations to get away from |00> state".
   */
  Circuit prepend_2qb_circuit;

  /** The same initial ops as for "prepend_2qb_circuit",
   *  which is the case N=2, but for N >= 2 qubits.
   *  @param qubits The number of qubits in the circuit; must be >= 2.
   *  @return A circuit with rotations to get away from |00> state.
   */
  static Circuit get_prepend_circuit(unsigned qubits);

  /** Adds the same ops as get_prepend_circuit and prepend_2qb_circuit
   *  to an already constructed circuit.
   */
  static void add_initial_prepend_ops(Circuit& circ);

  /** Constructed once and shared between tests.
   *  The caller need not make a copy if not required.
   */
  static const CircuitsForTesting& get();

  CircuitsForTesting();
};

}  // namespace tket
