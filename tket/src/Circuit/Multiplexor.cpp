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
    const UniformQControlBox::op_map_t &op_map, unsigned n_controls,
    unsigned n_targets) {
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

UniformQControlBox::UniformQControlBox(const op_map_t &op_map)
    : Box(OpType::UniformQControlBox), op_map_(op_map) {
  auto it = op_map.begin();
  if (it == op_map.end()) {
    throw std::invalid_argument("No Ops provided.");
  }
  n_controls_ = 0;
  n_targets_ = 0;
  for (; it != op_map.end(); it++) {
    op_signature_t op_sig = it->second->get_signature();
    if ((unsigned long)std::count(
            op_sig.begin(), op_sig.end(), EdgeType::Quantum) != op_sig.size()) {
      throw BadOpType(
          "Quantum control of classical wires not supported",
          it->second->get_type());
    }
    if (it == op_map.begin()) {
      n_controls_ = (unsigned)it->first.size();
      if (n_controls_ > MAX_N_CONTROLS) {
        throw std::invalid_argument(
            "Bitstrings longer than " + std::to_string(MAX_N_CONTROLS) +
            " are not supported.");
      }
      n_targets_ = (unsigned)op_sig.size();
    } else {
      if (it->first.size() != n_controls_) {
        throw std::invalid_argument("Bitstrings must have the same width.");
      }
      if (op_sig.size() != n_targets_) {
        throw std::invalid_argument("Ops must have the same width.");
      }
    }
  }
}

UniformQControlBox::UniformQControlBox(const UniformQControlBox &other)
    : Box(other),
      n_controls_(other.n_controls_),
      n_targets_(other.n_targets_),
      op_map_(other.op_map_) {}

Op_ptr UniformQControlBox::symbol_substitution(
    const SymEngine::map_basic_basic &sub_map) const {
  op_map_t new_op_map;
  for (auto it = op_map_.begin(); it != op_map_.end(); it++) {
    new_op_map.insert({it->first, it->second->symbol_substitution(sub_map)});
  }
  return std::make_shared<UniformQControlBox>(new_op_map);
}

SymSet UniformQControlBox::free_symbols() const {
  SymSet all_symbols;
  for (auto it = op_map_.begin(); it != op_map_.end(); it++) {
    SymSet op_symbols = it->second->free_symbols();
    all_symbols.insert(op_symbols.begin(), op_symbols.end());
  }
  return all_symbols;
}

void UniformQControlBox::generate_circuit() const {
  circ_ = std::make_shared<Circuit>(
      multiplexor_sequential_decomp(op_map_, n_controls_, n_targets_));
}

Op_ptr UniformQControlBox::dagger() const {
  op_map_t new_op_map;
  for (auto it = op_map_.begin(); it != op_map_.end(); it++) {
    new_op_map.insert({it->first, it->second->dagger()});
  }
  return std::make_shared<UniformQControlBox>(new_op_map);
}

Op_ptr UniformQControlBox::transpose() const {
  op_map_t new_op_map;
  for (auto it = op_map_.begin(); it != op_map_.end(); it++) {
    new_op_map.insert({it->first, it->second->transpose()});
  }
  return std::make_shared<UniformQControlBox>(new_op_map);
}
}  // namespace tket