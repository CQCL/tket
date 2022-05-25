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

#include <Circuit/Circuit.hpp>
#include <Transformations/OptimisationPass.hpp>
#include <Utils/Assert.hpp>
#include <iostream>

using namespace tket;

int main() {
  Circuit circ(4);

  circ.add_op<unsigned>(OpType::CZ, {0, 2});
  circ.add_op<unsigned>(OpType::CZ, {3, 1});
  circ.add_op<unsigned>(OpType::V, {2});
  circ.add_op<unsigned>(OpType::V, {3});
  circ.add_op<unsigned>(OpType::CZ, {0, 3});
  circ.add_op<unsigned>(OpType::V, {3});
  circ.add_op<unsigned>(OpType::CZ, {3, 1});
  circ.add_op<unsigned>(OpType::CZ, {2, 1});
  circ.add_op<unsigned>(OpType::V, {2});
  circ.add_op<unsigned>(OpType::CZ, {0, 2});
  circ.add_op<unsigned>(OpType::X, {2});
  circ.add_op<unsigned>(OpType::V, {1});
  circ.add_op<unsigned>(OpType::CZ, {3, 1});
  circ.add_op<unsigned>(OpType::CZ, {2, 1});
  circ.add_op<unsigned>(OpType::CZ, {3, 1});
  circ.add_op<unsigned>(OpType::V, {2});
  circ.add_op<unsigned>(OpType::V, {1});
  circ.add_op<unsigned>(OpType::CZ, {2, 1});
  circ.add_op<unsigned>(OpType::X, {2});
  circ.add_op<unsigned>(OpType::CZ, {2, 1});
  circ.add_op<unsigned>(OpType::V, {2});
  circ.add_op<unsigned>(OpType::CZ, {2, 1});
  circ.add_op<unsigned>(OpType::CZ, {0, 2});
  circ.add_op<unsigned>(OpType::CZ, {2, 1});

  Transforms::clifford_simp().apply(circ);

  unsigned n = circ.n_qubits();
  TKET_ASSERT(n == 4);
  Circuit newcirc;  // TKET-800
  unsigned n_zero = newcirc.n_qubits();
  TKET_ASSERT(n_zero == 0);

  std::cout << "success" << std::endl;
}
