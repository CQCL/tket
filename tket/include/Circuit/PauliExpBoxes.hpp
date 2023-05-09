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
 * Operation defined as the exponential of a tensor of Pauli operators
 */
class PauliExpBox : public Box {
 public:
  /**
   * The operation implements the unitary operator
   * \f$ e^{-\frac12 i \pi t \sigma_0 \otimes \sigma_1 \otimes \cdots} \f$
   * where \f$ \sigma_i \in \{I,X,Y,Z\} \f$ are the Pauli operators.
   */
  PauliExpBox(const std::vector<Pauli> &paulis, const Expr &t, CXConfigType cx_config_type = CXConfigType::Tree);

  /**
   * Construct from the empty vector
   */
  PauliExpBox();

  /**
   * Copy constructor
   */
  PauliExpBox(const PauliExpBox &other);

  ~PauliExpBox() override = default;

  bool is_clifford() const override;

  SymSet free_symbols() const override;

  /**
   * Equality check between two PauliExpBox instances
   */
  bool is_equal(const Op &op_other) const override {
    const auto &other = dynamic_cast<const PauliExpBox &>(op_other);
    return id_ == other.get_id();
  }

  /** Get the Pauli string */
  std::vector<Pauli> get_paulis() const { return paulis_; }

  /** Get the phase parameter */
  Expr get_phase() const { return t_; }

  /** Get the cx_config parameter (affects box decomposition) */
  CXConfigType get_cx_config() const { return cx_config_; }

  Op_ptr dagger() const override;

  Op_ptr transpose() const override;

  Op_ptr symbol_substitution(
      const SymEngine::map_basic_basic &sub_map) const override;

  static Op_ptr from_json(const nlohmann::json &j);

  static nlohmann::json to_json(const Op_ptr &op);

 protected:
  void generate_circuit() const override;

 private:
  std::vector<Pauli> paulis_;
  Expr t_;
  CXConfigType cx_config_;
};

// class PairwisePauliExpBox : public Box {
//  public:
//   /**
//    * The operation implements the unitary operator
//    * \f$ e^{-\frac12 i \pi t \sigma_0 \otimes \sigma_1 \otimes \cdots} \f$
//    * where \f$ \sigma_i \in \{I,X,Y,Z\} \f$ are the Pauli operators.
//    */
//   PairwisePauliExpBox(const std::vector<Pauli> &paulis, const Expr &t);
//
//   /**
//    * Construct from the empty vector
//    */
//   PairwisePauliExpBox();
//
//   /**
//    * Copy constructor
//    */
//   PairwisePauliExpBox(const PairwisePauliExpBox &other);
//
//   ~PairwisePauliExpBox() override {}
//
//   /** Get the Pauli string */
//   std::vector<Pauli> get_paulis() const { return paulis_; }
//
//   /** Get the phase parameter */
//   Expr get_phase() const { return t_; }
//
//   Op_ptr dagger() const override;
//
//   Op_ptr transpose() const override;
//
//   Op_ptr symbol_substitution(
//       const SymEngine::map_basic_basic &sub_map) const override;
//
//   static Op_ptr from_json(const nlohmann::json &j);
//
//   static nlohmann::json to_json(const Op_ptr &op);
//
//  protected:
//   void generate_circuit() const override;
//
//  private:
//   std::vector<Pauli> paulis_;
//   Expr t_;
// };
//
// class CommutingPauliExpBox : public Box {
//  public:
//   /**
//    * The operation implements the unitary operator
//    * \f$ e^{-\frac12 i \pi t \sigma_0 \otimes \sigma_1 \otimes \cdots} \f$
//    * where \f$ \sigma_i \in \{I,X,Y,Z\} \f$ are the Pauli operators.
//    */
//   CommutingPauliExpBox(const std::vector<Pauli> &paulis, const Expr &t);
//
//   /**
//    * Construct from the empty vector
//    */
//   CommutingPauliExpBox();
//
//   /**
//    * Copy constructor
//    */
//   CommutingPauliExpBox(const CommutingPauliExpBox &other);
//
//   ~CommutingPauliExpBox() override {}
// };

}  // namespace tket
