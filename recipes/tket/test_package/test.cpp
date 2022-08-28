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
#include <iostream>
#include <tkassert/Assert.hpp>

using namespace tket;

int main() {
  Circuit circ(2);

  circ.add_op<unsigned>(OpType::H, {0});
  circ.add_op<unsigned>(OpType::CX, {0, 1});
  circ.add_op<unsigned>(OpType::Rz, 0.2, {1});
  circ.add_op<unsigned>(OpType::CX, {0, 1});

  Transforms::synthesise_tk().apply(circ);
  std::cout << circ << std::endl;
  std::cout << "============== (SythesiseTK)" << std::endl;

  Transforms::two_qubit_squash(OpType::TK2).apply(circ);
  std::cout << circ << std::endl;
  std::cout << "============== (KAK Decomposition)" << std::endl;
}
