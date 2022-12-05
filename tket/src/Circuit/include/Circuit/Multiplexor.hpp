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

#include "Boxes.hpp"
#include "Circuit.hpp"

namespace tket {

/**
 * @brief Map bitstrings to Ops
 *
 *
 *
 */
typedef std::map<std::vector<bool>, Op_ptr> ctrl_op_map_t;

/**
 * Uniformly controlled quantum ops
 */
class UniformQControlBox : public Box {
 public:
  /**
   * Construct from a given ctrl_op_map_t.
   *
   * The ctrl_op_map_t must meet the following requirements:
   * 1. all bitstrings have the same length
   * 		and the length should be no bigger than 32
   * 2. all ops have the same number of wires
   * 3. all ops only have quantum wires
   *
   * @param op_map
   */
  explicit UniformQControlBox(const ctrl_op_map_t &op_map);

  /**
   * Copy constructor
   */
  UniformQControlBox(const UniformQControlBox &other);

  ~UniformQControlBox() override {}

  Op_ptr symbol_substitution(
      const SymEngine::map_basic_basic &sub_map) const override;

  SymSet free_symbols() const override;

  ctrl_op_map_t get_ops() const { return op_map_; }

  /**
   * Equality check between two UniformQControlBox instances
   */
  bool is_equal(const Op &op_other) const override {
    const UniformQControlBox &other =
        dynamic_cast<const UniformQControlBox &>(op_other);
    return id_ == other.get_id();
  }

  Op_ptr dagger() const override;

  Op_ptr transpose() const override;

 protected:
  /**
   * @brief Implement the multiplexor naively using X gates and QControlBoxes
   *
   */
  void generate_circuit() const override;
  UniformQControlBox()
      : Box(OpType::UniformQControlBox),
        n_controls_(0),
        n_targets_(0),
        op_map_() {}

 private:
  unsigned n_controls_;
  unsigned n_targets_;
  ctrl_op_map_t op_map_;
};

/**
 * Uniformly controlled single-axis rotations
 */
class UniformQControlRotationBox : public Box {
 public:
  /**
   * @brief Construct from a op_map. All ops should be of the same type.
   * They can be either Rx, Ry or Rz.
   *
   * @param op_map
   */
  explicit UniformQControlRotationBox(const ctrl_op_map_t &op_map);
  /**
   * Copy constructor
   */
  UniformQControlRotationBox(const UniformQControlRotationBox &other);

  ~UniformQControlRotationBox() override {}

  Op_ptr symbol_substitution(
      const SymEngine::map_basic_basic &sub_map) const override;

  SymSet free_symbols() const override;

  ctrl_op_map_t get_ops() const { return op_map_; }

  /**
   * Equality check between two UniformQControlBox instances
   */
  bool is_equal(const Op &op_other) const override {
    const UniformQControlRotationBox &other =
        dynamic_cast<const UniformQControlRotationBox &>(op_other);
    return id_ == other.get_id();
  }

  Op_ptr dagger() const override;

  Op_ptr transpose() const override;

 protected:
  /**
   * @brief Implement uniformly controlled same-axis rotations (UCR)
   * with 2^ctrl_qubits SQ rotations, 2^ctrl_qubits CXs, and 2 H gates for
   * X-axis rotations.
   *
   * https://arxiv.org/abs/quant-ph/0410066
   */
  void generate_circuit() const override;

 private:
  unsigned n_controls_;
  unsigned n_targets_;
  ctrl_op_map_t op_map_;
  OpType axis_;
};

}  // namespace tket