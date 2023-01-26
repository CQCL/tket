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
#include "Circuit.hpp"
#include "Utils/Json.hpp"

namespace tket {
/**
 * Box to synthesise a statevector
 */
class StatePreparationBox : public Box {
 public:
  /**
   * Construct a circuit that prepares the given statevector
   * from the computational-basis zero state
   *
   * @param statevector a normalised statevector
   * @param is_inverse whether to output the inverse of the preparation circuit
   */
  explicit StatePreparationBox(
      const Eigen::VectorXcd &statevector, bool is_inverse = false);

  /**
   * Copy constructor
   */
  StatePreparationBox(const StatePreparationBox &other);
  ~StatePreparationBox() override {}

  Op_ptr symbol_substitution(
      const SymEngine::map_basic_basic &) const override {
    return Op_ptr();
  }

  SymSet free_symbols() const override { return {}; }

  /**
   * Equality check between two StatePreparationBox instances
   */
  bool is_equal(const Op &op_other) const override {
    const StatePreparationBox &other =
        dynamic_cast<const StatePreparationBox &>(op_other);
    return id_ == other.get_id();
  }

  Op_ptr dagger() const override;

  op_signature_t get_signature() const override;

  static Op_ptr from_json(const nlohmann::json &j);

  static nlohmann::json to_json(const Op_ptr &op);

  Eigen::VectorXcd get_statevector() const;

  bool is_inverse() const;

 protected:
  /**
   * @brief Generate the state preparation circuit using
   * multiplexed-Ry gates and multiplexed-Rz gates
   *
   */
  void generate_circuit() const override;

  StatePreparationBox()
      : Box(OpType::StatePreparationBox), statevector_(), is_inverse_(false) {}

 private:
  const Eigen::VectorXcd statevector_;
  const bool is_inverse_;
};
}  // namespace tket