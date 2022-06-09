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
#include <Transformations/BasicOptimisation.hpp>
#include <Transformations/OptimisationPass.hpp>
#include <Transformations/PauliOptimisation.hpp>
#include <Utils/Assert.hpp>
#include <cmath>
#include <iostream>

#include "Utils/Expression.hpp"

using namespace tket;

int main() {
  Circuit circ(3);

  for (unsigned i = 0; i < circ.n_qubits(); ++i) {
    circ.add_op<unsigned>(OpType::Rz, 0.3, {i});
  }
  circ.add_op<unsigned>(OpType::Rz, 0.4, {2});
  circ.add_op<unsigned>(OpType::NPhasedX, {0.5, 0.2}, {0, 1});
  for (unsigned i = 0; i < circ.n_qubits(); ++i) {
    circ.add_op<unsigned>(OpType::Rz, i * 0.2, {i});
  }

  std::cout << circ << std::endl;

  Transforms::absorb_Rz_NPhasedX().apply(circ);

  std::cout << circ << std::endl;
  std::cout << "n gates: " << circ.count_gates(OpType::NPhasedX) << std::endl;
  circ.to_graphviz_file("test.dot");
}
