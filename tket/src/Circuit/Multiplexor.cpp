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

#include "Circuit/Multiplexor.hpp"

namespace tket {

static const unsigned MAX_N_CONTROLS = 32;

/**
 * @brief Implement a multiplexor by sequentially applying QControlBoxes.
 * Assume all ops have width n_targets, and all bitstrings have size n_controls
 * @param op_map
 * @param n_controls
 * @param n_targets
 * @return Circuit
 */
static Circuit multiplexor_sequential_decomp(
    const ctrl_op_map_t &op_map, unsigned n_controls, unsigned n_targets) {
  Circuit c(n_controls + n_targets);
  std::vector<unsigned> qubits(n_controls + n_targets);
  std::iota(std::begin(qubits), std::end(qubits), 0);
  for (auto it = op_map.begin(); it != op_map.end(); it++) {
    std::vector<unsigned> zero_ctrls;
    for (unsigned i = 0; i < it->first.size(); i++) {
      if (!it->first[i]) {
        zero_ctrls.push_back(i);
        c.add_op<unsigned>(OpType::X, {i});
      }
    }
    QControlBox qcbox(it->second, n_controls);
    c.add_box(qcbox, qubits);
    for (unsigned i : zero_ctrls) {
      c.add_op<unsigned>(OpType::X, {i});
    }
  }
  return c;
}

/**
 * @brief Convert an unsigned to its binary representation
 * e.g. 4 -> [1, 0, 0]
 *
 * @param dec
 * @param width used to pad zeros
 * @return std::vector<bool>
 */
static std::vector<bool> dec_to_bin(unsigned dec, unsigned width) {
  auto bs = std::bitset<MAX_N_CONTROLS>(dec);
  std::vector<bool> bits(width);
  for (unsigned i = 0; i < width; i++) {
    bits[width - i - 1] = bs[i];
  }
  return bits;
}

/**
 * @brief Implement uniformly controlled same-axis rotations (UCR)
 * with 2^ctrl_qubits SQ rotations, 2^ctrl_qubits CXs, and 2 H gates for X-axis
 * rotations.
 *
 * https://arxiv.org/abs/quant-ph/0410066
 * This is a special case derived from equation (3)
 * A UCR gate controlled by n qubits have the decomposition UCR = CX P CX Q
 * (multiplication order), where P and Q are themselves UCR gates controlled by
 * n-1 qubits.
 *
 * Also notice that CX P CX Q = Q CX P CX, therefore we can control the
 * direction of each decomposition to avoid adding adjacent CX gates.
 * e.g. UCR = CX P CX Q = (CX Q' CX P' CX*) CX (CX* P'' CX Q'')
 * The two CX* can be cancelled,
 * hence UCR = CX P CX Q = (CX Q' CX P') CX (P''CX Q'')
 *
 *
 * @param angles list of 2^ctrl_qubits angles, angles[i] is the angle activated
 * by bitstring binary(i)
 * @param axis can be either Ry or Rz
 * @param total_qubits the total number of qubits in the final output circuit
 * @param circ circuit to update
 * @param direction controls the decomposition direction of each demultiplex
 * step. Will be implemented as P CX Q if true, and Q CX P if false. If nullopt,
 * CX P CX Q will be implemented as the root step.
 *
 */
static void recursive_demultiplex_rotation(
    const std::vector<Expr> &angles, const OpType &axis, unsigned total_qubits,
    Circuit &circ, std::optional<bool> direction) {
  unsigned n_rotations = angles.size();
  unsigned n_qubits = (unsigned)log2(n_rotations) + 1;
  unsigned mid = (unsigned)(n_rotations / 2);
  std::vector<Expr> p_angles;
  std::vector<Expr> q_angles;
  for (unsigned i = 0; i < mid; i++) {
    p_angles.push_back((angles[i] - angles[mid + i]) / 2);
    q_angles.push_back((angles[i] + angles[mid + i]) / 2);
  }
  if (direction != std::nullopt && !direction.value()) {
    std::swap(p_angles, q_angles);
  }
  if (q_angles.size() == 1) {
    circ.add_op<unsigned>(axis, q_angles[0], {total_qubits - 1});
  } else {
    recursive_demultiplex_rotation(q_angles, axis, total_qubits, circ, true);
  }
  circ.add_op<unsigned>(
      OpType::CX, {total_qubits - n_qubits, total_qubits - 1});
  if (p_angles.size() == 1) {
    circ.add_op<unsigned>(axis, p_angles[0], {total_qubits - 1});
  } else {
    recursive_demultiplex_rotation(p_angles, axis, total_qubits, circ, false);
  }
  if (direction == std::nullopt) {
    circ.add_op<unsigned>(
        OpType::CX, {total_qubits - n_qubits, total_qubits - 1});
  }
}

static void op_map_validate(const ctrl_op_map_t &op_map) {
  unsigned n_controls = 0;
  unsigned n_targets = 0;
  for (auto it = op_map.begin(); it != op_map.end(); it++) {
    op_signature_t op_sig = it->second->get_signature();
    if ((unsigned long)std::count(
            op_sig.begin(), op_sig.end(), EdgeType::Quantum) != op_sig.size()) {
      throw BadOpType(
          "Quantum control of classical wires not supported",
          it->second->get_type());
    }
    if (it == op_map.begin()) {
      n_controls = (unsigned)it->first.size();
      if (n_controls > MAX_N_CONTROLS) {
        throw std::invalid_argument(
            "Bitstrings longer than " + std::to_string(MAX_N_CONTROLS) +
            " are not supported.");
      }
      n_targets = (unsigned)op_sig.size();
    } else {
      if (it->first.size() != n_controls) {
        throw std::invalid_argument("Bitstrings must have the same width.");
      }
      if (op_sig.size() != n_targets) {
        throw std::invalid_argument("Ops must have the same width.");
      }
    }
  }
}

static ctrl_op_map_t op_map_symbol_sub(
    const SymEngine::map_basic_basic &sub_map, const ctrl_op_map_t &op_map) {
  ctrl_op_map_t new_op_map;
  for (auto it = op_map.begin(); it != op_map.end(); it++) {
    new_op_map.insert({it->first, it->second->symbol_substitution(sub_map)});
  }
  return new_op_map;
}

static SymSet op_map_free_symbols(const ctrl_op_map_t &op_map) {
  SymSet all_symbols;
  for (auto it = op_map.begin(); it != op_map.end(); it++) {
    SymSet op_symbols = it->second->free_symbols();
    all_symbols.insert(op_symbols.begin(), op_symbols.end());
  }
  return all_symbols;
}

static ctrl_op_map_t op_map_dagger(const ctrl_op_map_t &op_map) {
  ctrl_op_map_t new_op_map;
  for (auto it = op_map.begin(); it != op_map.end(); it++) {
    new_op_map.insert({it->first, it->second->dagger()});
  }
  return new_op_map;
}

static ctrl_op_map_t op_map_transpose(const ctrl_op_map_t &op_map) {
  ctrl_op_map_t new_op_map;
  for (auto it = op_map.begin(); it != op_map.end(); it++) {
    new_op_map.insert({it->first, it->second->transpose()});
  }
  return new_op_map;
}

UniformQControlBox::UniformQControlBox(const ctrl_op_map_t &op_map)
    : Box(OpType::UniformQControlBox), op_map_(op_map) {
  auto it = op_map.begin();
  if (it == op_map.end()) {
    throw std::invalid_argument("No Ops provided.");
  }
  n_controls_ = (unsigned)it->first.size();
  n_targets_ = it->second->n_qubits();
  op_map_validate(op_map);
}

UniformQControlBox::UniformQControlBox(const UniformQControlBox &other)
    : Box(other),
      n_controls_(other.n_controls_),
      n_targets_(other.n_targets_),
      op_map_(other.op_map_) {}

Op_ptr UniformQControlBox::symbol_substitution(
    const SymEngine::map_basic_basic &sub_map) const {
  ctrl_op_map_t new_op_map = op_map_symbol_sub(sub_map, op_map_);
  return std::make_shared<UniformQControlBox>(new_op_map);
}

SymSet UniformQControlBox::free_symbols() const {
  return op_map_free_symbols(op_map_);
}

Op_ptr UniformQControlBox::dagger() const {
  return std::make_shared<UniformQControlBox>(op_map_dagger(op_map_));
}

Op_ptr UniformQControlBox::transpose() const {
  return std::make_shared<UniformQControlBox>(op_map_transpose(op_map_));
}

void UniformQControlBox::generate_circuit() const {
  circ_ = std::make_shared<Circuit>(
      multiplexor_sequential_decomp(op_map_, n_controls_, n_targets_));
}

UniformQControlRotationBox::UniformQControlRotationBox(
    const ctrl_op_map_t &op_map)
    : Box(OpType::UniformQControlRotationBox), op_map_(op_map) {
  auto it = op_map.begin();
  if (it == op_map.end()) {
    throw std::invalid_argument("No Ops provided.");
  }
  n_targets_ = 1;
  for (; it != op_map.end(); it++) {
    if (it == op_map.begin()) {
      n_controls_ = (unsigned)it->first.size();
      axis_ = it->second->get_type();
      if (axis_ != OpType::Rx && axis_ != OpType::Ry && axis_ != OpType::Rz) {
        throw std::invalid_argument("Ops must be either Rx, Ry, or Rz.");
      }
    } else {
      if (it->second->get_type() != axis_) {
        throw std::invalid_argument("Ops must have the same rotation type.");
      }
    }
  }
  op_map_validate(op_map);
}

UniformQControlRotationBox::UniformQControlRotationBox(
    const UniformQControlRotationBox &other)
    : Box(other),
      n_controls_(other.n_controls_),
      n_targets_(other.n_targets_),
      op_map_(other.op_map_),
      axis_(other.axis_) {}

Op_ptr UniformQControlRotationBox::symbol_substitution(
    const SymEngine::map_basic_basic &sub_map) const {
  ctrl_op_map_t new_op_map = op_map_symbol_sub(sub_map, op_map_);
  return std::make_shared<UniformQControlRotationBox>(new_op_map);
}

SymSet UniformQControlRotationBox::free_symbols() const {
  return op_map_free_symbols(op_map_);
}

Op_ptr UniformQControlRotationBox::dagger() const {
  return std::make_shared<UniformQControlRotationBox>(op_map_dagger(op_map_));
}

Op_ptr UniformQControlRotationBox::transpose() const {
  return std::make_shared<UniformQControlRotationBox>(
      op_map_transpose(op_map_));
}

void UniformQControlRotationBox::generate_circuit() const {
  Circuit circ(n_controls_ + 1);
  if (n_controls_ == 0) {
    auto it = op_map_.begin();
    circ.add_op<unsigned>(axis_, it->second->get_params()[0], {0});
    circ_ = std::make_shared<Circuit>(circ);
    return;
  }

  std::vector<Expr> rotations(1 << n_controls_);
  // convert op_map to a vector of 2^n_controls_ angles
  for (unsigned i = 0; i < (1 << n_controls_); i++) {
    auto it = op_map_.find(dec_to_bin(i, n_controls_));
    if (it == op_map_.end()) {
      rotations[i] = 0;
    } else {
      rotations[i] = it->second->get_params()[0];
    }
  }
  OpType axis = axis_;
  if (axis_ == OpType::Rx) {
    circ.add_op<unsigned>(OpType::H, {n_controls_});
    axis = OpType::Rz;
  }
  recursive_demultiplex_rotation(
      rotations, axis, n_controls_ + 1, circ, std::nullopt);
  if (axis_ == OpType::Rx) {
    circ.add_op<unsigned>(OpType::H, {n_controls_});
  }
  circ_ = std::make_shared<Circuit>(circ);
}

}  // namespace tket