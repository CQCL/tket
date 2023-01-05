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
#include "Utils/Json.hpp"

namespace tket {

/**
 * @brief Map bitstrings to Ops
 *
 */
typedef std::map<std::vector<bool>, Op_ptr> ctrl_op_map_t;

/**
 * Multiplexed ops
 */
class MultiplexorBox : public Box {
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
  explicit MultiplexorBox(const ctrl_op_map_t &op_map);

  /**
   * Copy constructor
   */
  MultiplexorBox(const MultiplexorBox &other);

  ~MultiplexorBox() override {}

  Op_ptr symbol_substitution(
      const SymEngine::map_basic_basic &sub_map) const override;

  SymSet free_symbols() const override;

  ctrl_op_map_t get_ops() const { return op_map_; }

  /**
   * Equality check between two MultiplexorBox instances
   */
  bool is_equal(const Op &op_other) const override {
    const MultiplexorBox &other =
        dynamic_cast<const MultiplexorBox &>(op_other);
    return id_ == other.get_id();
  }

  Op_ptr dagger() const override;

  Op_ptr transpose() const override;

  op_signature_t get_signature() const override;

  ctrl_op_map_t get_op_map() const { return op_map_; }

  static Op_ptr from_json(const nlohmann::json &j);

  static nlohmann::json to_json(const Op_ptr &op);

 protected:
  /**
   * @brief Implement the multiplexor naively using X gates and QControlBoxes
   *
   */
  void generate_circuit() const override;
  MultiplexorBox()
      : Box(OpType::MultiplexorBox), n_controls_(0), n_targets_(0), op_map_() {}

 private:
  unsigned n_controls_;
  unsigned n_targets_;
  ctrl_op_map_t op_map_;
};

/**
 * Multiplexed single-axis rotations
 */
class MultiplexedRotationBox : public Box {
 public:
  /**
   * @brief Construct from a op_map. All ops should be of the same type.
   * They can be either Rx, Ry or Rz.
   *
   * @param op_map
   */
  explicit MultiplexedRotationBox(const ctrl_op_map_t &op_map);
  /**
   * Copy constructor
   */
  MultiplexedRotationBox(const MultiplexedRotationBox &other);

  ~MultiplexedRotationBox() override {}

  Op_ptr symbol_substitution(
      const SymEngine::map_basic_basic &sub_map) const override;

  SymSet free_symbols() const override;

  ctrl_op_map_t get_ops() const { return op_map_; }

  /**
   * Equality check between two MultiplexedRotationBox instances
   */
  bool is_equal(const Op &op_other) const override {
    const MultiplexedRotationBox &other =
        dynamic_cast<const MultiplexedRotationBox &>(op_other);
    return id_ == other.get_id();
  }

  Op_ptr dagger() const override;

  Op_ptr transpose() const override;

  op_signature_t get_signature() const override;

  ctrl_op_map_t get_op_map() const { return op_map_; }

  static Op_ptr from_json(const nlohmann::json &j);

  static nlohmann::json to_json(const Op_ptr &op);

 protected:
  /**
   * @brief Implement multiplexed rotation
   * (i.e. uniformly controlled same-axis rotations (UCR)),
   * with 2^ctrl_qubits SQ rotations, 2^ctrl_qubits CXs, and 2 H gates for
   * X-axis rotations.
   *
   * https://arxiv.org/abs/quant-ph/0410066
   */
  void generate_circuit() const override;

 private:
  unsigned n_controls_;
  ctrl_op_map_t op_map_;
  OpType axis_;
};

/**
 * Multiplexed U2 gate
 */
class MultiplexedU2Box : public Box {
 public:
  /**
   * @brief Construct from a op_map. Ops must be single-qubit unitary gate types
   * or Unitary1QBox.
   *
   * @param op_map
   * @param impl_diag whether to implement the final MultiplexedRotationBox,
   * default to true
   */
  explicit MultiplexedU2Box(const ctrl_op_map_t &op_map, bool impl_diag = true);
  /**
   * Copy constructor
   */
  MultiplexedU2Box(const MultiplexedU2Box &other);

  ~MultiplexedU2Box() override {}

  Op_ptr symbol_substitution(
      const SymEngine::map_basic_basic &sub_map) const override;

  SymSet free_symbols() const override;

  ctrl_op_map_t get_ops() const { return op_map_; }

  /**
   * Equality check between two MultiplexedU2Box instances
   */
  bool is_equal(const Op &op_other) const override {
    const MultiplexedU2Box &other =
        dynamic_cast<const MultiplexedU2Box &>(op_other);
    return id_ == other.get_id();
  }

  Op_ptr dagger() const override;

  Op_ptr transpose() const override;

  op_signature_t get_signature() const override;

  ctrl_op_map_t get_op_map() const { return op_map_; }

  bool get_impl_diag() const { return impl_diag_; }

  static Op_ptr from_json(const nlohmann::json &j);

  static nlohmann::json to_json(const Op_ptr &op);

 protected:
  /**
   * @brief Implement multiplexed U2 gate
   * (i.e. uniformly controlled U2 gate (UCU2))
   * with 2^ctrl_qubits SQ gates, 2^ctrl_qubits CXs, and a
   * MultiplexedRotationBox at the end
   *
   * https://arxiv.org/abs/quant-ph/0410066
   */
  void generate_circuit() const override;

 private:
  unsigned n_controls_;
  ctrl_op_map_t op_map_;
  bool impl_diag_;
};
}  // namespace tket