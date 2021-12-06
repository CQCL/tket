// Copyright 2019-2021 Cambridge Quantum Computing
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

#include "Circuit/Boxes.hpp"
#include "Converters.hpp"

namespace tket {

class UnitaryTableauBox : public Box {
 public:
  explicit UnitaryTableauBox(const UnitaryTableau& tab);
  explicit UnitaryTableauBox(
      const MatrixXb& xx, const MatrixXb& xz, const VectorXb& xph,
      const MatrixXb& zx, const MatrixXb& zz, const VectorXb& zph);

  Op_ptr dagger() const override;
  Op_ptr transpose() const override;

  Op_ptr symbol_substitution(
      const SymEngine::map_basic_basic& sub_map) const override;

  SymSet free_symbols() const override;

  bool is_equal(const Op& op_other) const override;

  const UnitaryTableau& get_tableau() const;

  static Op_ptr from_json(const nlohmann::json& j);
  static nlohmann::json to_json(const Op_ptr& op);

  op_signature_t get_signature() const override;

 protected:
  void generate_circuit() const override;

 private:
  UnitaryTableau tab_;
};

}  // namespace tket
