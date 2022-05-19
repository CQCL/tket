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
#include <Transformations/Decomposition.hpp>
#include <Transformations/OptimisationPass.hpp>
#include <Utils/Assert.hpp>
#include <iostream>
#include <optional>

using namespace tket;

int main() {
  Circuit circ(2);

  circ.add_op<unsigned>(OpType::TK2, {0.4, 0.1, 0.}, {0, 1});
  // circ.add_op<unsigned>(OpType::XXPhase, 0.2, {0, 1});
  // circ.add_op<unsigned>(OpType::YYPhase, 0.2, {0, 1});

  // Transforms::decompose_ZZPhase().apply(circ);
  Transforms::TwoQbFidelities fid{0.99, std::nullopt, std::nullopt};
  Transforms::decompose_TK2(fid).apply(circ);

  std::cout << circ << std::endl;
}
