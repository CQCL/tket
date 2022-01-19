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

#include "CircuitsForTesting.hpp"

#include "Utils/Assert.hpp"

namespace tket {

const CircuitsForTesting& CircuitsForTesting::get() {
  static const CircuitsForTesting circuits;
  return circuits;
}

static void add_to_uccsd(Circuit& uccsd) {
  uccsd.add_op<unsigned>(OpType::Rx, 0.5, {0});
  uccsd.add_op<unsigned>(OpType::Rx, 0.5, {1});
  uccsd.add_op<unsigned>(OpType::H, {2});
  uccsd.add_op<unsigned>(OpType::Rx, 0.5, {3});
  uccsd.add_op<unsigned>(OpType::CX, {3, 2});
  uccsd.add_op<unsigned>(OpType::CX, {2, 1});
  uccsd.add_op<unsigned>(OpType::CX, {1, 0});
  uccsd.add_op<unsigned>(OpType::Rz, 0.356, {0});
  uccsd.add_op<unsigned>(OpType::CX, {1, 0});
  uccsd.add_op<unsigned>(OpType::CX, {2, 1});
  uccsd.add_op<unsigned>(OpType::CX, {3, 2});
  uccsd.add_op<unsigned>(OpType::Rx, 1.5, {0});
  uccsd.add_op<unsigned>(OpType::Rx, 1.5, {1});
  uccsd.add_op<unsigned>(OpType::H, {2});
  uccsd.add_op<unsigned>(OpType::Rx, 1.5, {3});
  uccsd.add_op<unsigned>(OpType::H, {0});
  uccsd.add_op<unsigned>(OpType::Rx, 0.5, {1});
  uccsd.add_op<unsigned>(OpType::Rx, 0.5, {2});
  uccsd.add_op<unsigned>(OpType::Rx, 0.5, {3});
  uccsd.add_op<unsigned>(OpType::CX, {3, 2});
  uccsd.add_op<unsigned>(OpType::CX, {2, 1});
  uccsd.add_op<unsigned>(OpType::CX, {1, 0});
  uccsd.add_op<unsigned>(OpType::Rz, 1.183, {0});
  uccsd.add_op<unsigned>(OpType::CX, {1, 0});
  uccsd.add_op<unsigned>(OpType::CX, {2, 1});
  uccsd.add_op<unsigned>(OpType::CX, {3, 2});
  uccsd.add_op<unsigned>(OpType::H, {0});
  uccsd.add_op<unsigned>(OpType::Rx, 1.5, {1});
  uccsd.add_op<unsigned>(OpType::Rx, 1.5, {2});
  uccsd.add_op<unsigned>(OpType::Rx, 1.5, {3});
}

void CircuitsForTesting::add_initial_prepend_ops(Circuit& circ) {
  TKET_ASSERT(circ.n_qubits() >= 2);
  circ.add_op<unsigned>(OpType::Rx, 0.333, {0});
  circ.add_op<unsigned>(OpType::Rz, 1.2, {0});
  circ.add_op<unsigned>(OpType::Rx, -0.1111, {1});
  circ.add_op<unsigned>(OpType::Rz, 0.973, {1});
}

Circuit CircuitsForTesting::get_prepend_circuit(unsigned qubits) {
  Circuit prepend(qubits);
  add_initial_prepend_ops(prepend);
  return prepend;
}

CircuitsForTesting::CircuitsForTesting() : uccsd(4) {
  add_to_uccsd(uccsd);
  prepend_2qb_circuit = get_prepend_circuit(2);
}

}  // namespace tket
