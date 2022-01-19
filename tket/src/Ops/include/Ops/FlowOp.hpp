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

#include "Op.hpp"

namespace tket {

class FlowOp : public Op {
 public:
  explicit FlowOp(OpType type, std::optional<std::string> label = std::nullopt);

  Op_ptr symbol_substitution(
      const SymEngine::map_basic_basic &sub_map) const override;

  SymSet free_symbols() const override;

  std::string get_name(bool latex = false) const override;

  op_signature_t get_signature() const override;

  std::optional<std::string> get_label() const;

  ~FlowOp() override;

  /**
   * Equality check between two FlowOp instances
   */
  bool is_equal(const Op &other) const override;

 protected:
  FlowOp();

 private:
  std::optional<std::string> label_;
};

}  // namespace tket
