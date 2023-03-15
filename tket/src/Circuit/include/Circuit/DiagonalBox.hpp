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
 * Box to synthesise a diagonal operator
 */
class DiagonalBox : public Box {
 public:
  /**
   * Construct a circuit that synthesise the given unitary diagonal operator
   *
   * @param diagonal the digonal entries of the operator
   */
  explicit DiagonalBox(const Eigen::VectorXcd &diagonal);

  /**
   * Copy constructor
   */
  DiagonalBox(const DiagonalBox &other);
  ~DiagonalBox() override {}

  Op_ptr symbol_substitution(
      const SymEngine::map_basic_basic &) const override {
    return Op_ptr();
  }

  SymSet free_symbols() const override { return {}; }

  /**
   * Equality check between two DiagonalBox instances
   */
  bool is_equal(const Op &op_other) const override {
    const DiagonalBox &other = dynamic_cast<const DiagonalBox &>(op_other);
    return id_ == other.get_id();
  }

  Op_ptr dagger() const override;
  Op_ptr transpose() const override;

  op_signature_t get_signature() const override;

  static Op_ptr from_json(const nlohmann::json &j);

  static nlohmann::json to_json(const Op_ptr &op);

  Eigen::VectorXcd get_diagonal() const;

 protected:
  /**
   * @brief Generate the decomposed circuit using
   * multiplexed-Rz gates
   *
   */
  void generate_circuit() const override;

  DiagonalBox() : Box(OpType::DiagonalBox), diagonal_() {}

 private:
  const Eigen::VectorXcd diagonal_;
};
}  // namespace tket