// Copyright 2019-2023 Cambridge Quantum Computing
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

#include "Boxes.hpp"

namespace tket {

/**
 * Exception indicating that dummy boxes cannot be decomposed.
 */
class DummyBoxNotDecomposable : public std::logic_error {
 public:
  DummyBoxNotDecomposable()
      : std::logic_error("Cannot generate circuit from DummyBox") {}
};

/**
 * @brief A placeholder operation that holds resource data.
 *
 * This box type cannot be decomposed into a circuit. It only serves to record
 * resource data for a region of a circuit: for example, upper and lower bounds
 * on gate counts and depth. A circuit containing such a box cannot be executed.
 */
class DummyBox : public Box {
 protected:
  /**
   * @brief Throw an exception.
   *
   * This box does not correspond to any actual circuit.
   *
   * @throws DummyBoxNotDecomposable
   */
  void generate_circuit() const override;
};

}  // namespace tket
