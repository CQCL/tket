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

class PauliExpBoxInvalidity : public std::logic_error {
 public:
  explicit PauliExpBoxInvalidity(const std::string &message)
      : std::logic_error(message) {}
};

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
  PauliExpBox(
      const std::vector<Pauli> &paulis, const Expr &t,
      CXConfigType cx_config_type = CXConfigType::Tree);

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

class PauliExpPairBox : public Box {
 public:
  PauliExpPairBox(
      const std::vector<Pauli> &paulis0, const Expr &t0,
      const std::vector<Pauli> &paulis1, const Expr &t1,
      CXConfigType cx_config_type = CXConfigType::Tree);

  /**
   * Construct from the empty vector
   */
  PauliExpPairBox();

  /**
   * Copy constructor
   */
  PauliExpPairBox(const PauliExpPairBox &other);

  ~PauliExpPairBox() override {}

  bool is_clifford() const override;

  SymSet free_symbols() const override;

  /**
   * Equality check between two instances
   */
  bool is_equal(const Op &op_other) const override {
    const auto &other = dynamic_cast<const PauliExpPairBox &>(op_other);
    return id_ == other.get_id();
  }

  /** Get Pauli strings for the pair */
  std::pair<std::vector<Pauli>, std::vector<Pauli>> get_paulis_pair() const {
    return std::make_pair(paulis0_, paulis1_);
  }

  /** Get phase parameters for the pair */
  std::pair<Expr, Expr> get_phase_pair() const {
    return std::make_pair(t0_, t1_);
  }

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
  std::vector<Pauli> paulis0_;
  Expr t0_;
  std::vector<Pauli> paulis1_;
  Expr t1_;
  CXConfigType cx_config_;
};

class PauliExpCommutingSetBox : public Box {
 public:
  PauliExpCommutingSetBox(
      const std::vector<std::pair<std::vector<Pauli>, Expr>> &pauli_gadgets,
      CXConfigType cx_config_type = CXConfigType::Tree);

  /**
   * Construct from the empty vector
   */
  PauliExpCommutingSetBox();

  /**
   * Copy constructor
   */
  PauliExpCommutingSetBox(const PauliExpCommutingSetBox &other);

  ~PauliExpCommutingSetBox() override {}

  bool paulis_commute() const;

  bool is_clifford() const override;

  SymSet free_symbols() const override;

  /**
   * Equality check between two instances
   */
  bool is_equal(const Op &op_other) const override {
    const auto &other = dynamic_cast<const PauliExpCommutingSetBox &>(op_other);
    return id_ == other.get_id();
  }

  /** Get the pauli gadgets */
  auto get_pauli_gadgets() const { return pauli_gadgets_; }

  /** Get the cx_config parameter (affects box decomposition) */
  auto get_cx_config() const { return cx_config_; }

  Op_ptr dagger() const override;

  Op_ptr transpose() const override;

  Op_ptr symbol_substitution(
      const SymEngine::map_basic_basic &sub_map) const override;

  static Op_ptr from_json(const nlohmann::json &j);

  static nlohmann::json to_json(const Op_ptr &op);

 protected:
  void generate_circuit() const override;

 private:
  std::vector<std::pair<std::vector<Pauli>, Expr>> pauli_gadgets_;
  CXConfigType cx_config_;
};

/**
 * Constructs a PauliExpBox for a single pauli gadget and appends it to a
 * circuit.
 *
 * @param circ The circuit to append the box to
 * @param pauli The pauli operator of the gadget
 * @param angle The rotation angle of the Pauli gadget in half-turns
 * @param cx_config The CX configuration to be used during synthesis
 */
void append_single_pauli_gadget_as_pauli_exp_box(
    Circuit &circ, const QubitPauliTensor &pauli, Expr angle,
    CXConfigType cx_config);

/**
 * Constructs a PauliExpPairBox for a pair of pauli gadgets and appends it to a
 * circuit. The pauli gadgets may or may not commute, so the ordering matters.
 *
 * @param circ The circuit to append the box to
 * @param pauli0 The pauli operator of the first gadget
 * @param angle0 The rotation angle of the first Pauli gadget in half-turns
 * @param pauli1 The pauli operator of the second gadget
 * @param angle1 The rotation angle of the second Pauli gadget in half-turns
 * @param cx_config The CX configuration to be used during synthesis
 */
void append_pauli_gadget_pair_as_box(
    Circuit &circ, const QubitPauliTensor &pauli0, Expr angle0,
    const QubitPauliTensor &pauli1, Expr angle1, CXConfigType cx_config);

/**
 * Constructs a PauliExpCommutingSetBox for a set of mutually commuting pauli
 * gadgets and appends it to a circuit. As the pauli gadgets all commute, the
 * ordering does not matter semantically, but may yield different synthesised
 * circuits.
 *
 * @param circ The circuit to append the box to
 * @param gadgets Description of the pauli gadgets as pairs of the angle of
 * rotation (in half-turns) and the Pauli operator the rotation acts about
 * @param cx_config The CX configuration to be used during synthesis
 */
void append_commuting_pauli_gadget_set_as_box(
    Circuit &circ, const std::list<std::pair<QubitPauliTensor, Expr>> &gadgets,
    CXConfigType cx_config);

}  // namespace tket
