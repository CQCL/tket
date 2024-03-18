// Copyright 2019-2024 Cambridge Quantum Computing
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
#include "tket/Diagonalisation/PauliPartition.hpp"
#include "tket/Transformations/PauliOptimisation.hpp"

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
   * where \f$ \sigma_i \in \{I,X,Y,Z\} \f$ are the Pauli operators and \f$ t
   * \f$ is the coefficient.
   */
  PauliExpBox(
      const SymPauliTensor &paulis,
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
  bool is_equal(const Op &op_other) const override;

  /** Get the Pauli string */
  std::vector<Pauli> get_paulis() const { return paulis_.string; }

  /** Get the phase parameter */
  Expr get_phase() const { return paulis_.coeff; }

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
  SymPauliTensor paulis_;
  CXConfigType cx_config_;
};

class PauliExpPairBox : public Box {
 public:
  PauliExpPairBox(
      const SymPauliTensor &paulis0, const SymPauliTensor &paulis1,
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
  bool is_equal(const Op &op_other) const override;

  /** Get Pauli strings for the pair */
  std::pair<std::vector<Pauli>, std::vector<Pauli>> get_paulis_pair() const {
    return std::make_pair(paulis0_.string, paulis1_.string);
  }

  /** Get phase parameters for the pair */
  std::pair<Expr, Expr> get_phase_pair() const {
    return std::make_pair(paulis0_.coeff, paulis1_.coeff);
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
  SymPauliTensor paulis0_;
  SymPauliTensor paulis1_;
  CXConfigType cx_config_;
};

class PauliExpCommutingSetBox : public Box {
 public:
  PauliExpCommutingSetBox(
      const std::vector<SymPauliTensor> &pauli_gadgets,
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
  bool is_equal(const Op &op_other) const override;

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
  std::vector<SymPauliTensor> pauli_gadgets_;
  CXConfigType cx_config_;
};

class TermSequenceBox : public Box {
 public:
  TermSequenceBox(
      const std::vector<SymPauliTensor> &pauli_gadgets,
      Transforms::PauliSynthStrat synth_strategy =
          Transforms::PauliSynthStrat::Sets,
      PauliPartitionStrat partition_strategy =
          PauliPartitionStrat::CommutingSets,
      GraphColourMethod graph_colouring = GraphColourMethod::Lazy,
      CXConfigType cx_configuration = CXConfigType::Tree);

  /**
   * Construct from the empty vector
   */
  TermSequenceBox();

  /**
   * Copy constructor
   */
  TermSequenceBox(const TermSequenceBox &other);

  ~TermSequenceBox() override {}

  bool is_clifford() const override;

  SymSet free_symbols() const override;

  /**
   * Equality check between two instances
   */
  bool is_equal(const Op &op_other) const override;

  /** Get the circuit synthesis strategy parameter (affects box decomposition)
   */
  Transforms::PauliSynthStrat get_synth_strategy() const;

  /** Get the pauli partitioning strategy parameter (affects box decomposition)
   */
  PauliPartitionStrat get_partition_strategy() const;

  /** Get the graph colouring parameter (affects box decomposition) */
  GraphColourMethod get_graph_colouring() const;

  /** Get the pauli gadgets */
  std::vector<SymPauliTensor> get_pauli_gadgets() const;

  /** Get the cx config parameter (affects box decomposition) */
  CXConfigType get_cx_config() const;

  Op_ptr dagger() const override;

  Op_ptr transpose() const override;

  Op_ptr symbol_substitution(
      const SymEngine::map_basic_basic &sub_map) const override;

  static Op_ptr from_json(const nlohmann::json &j);

  static nlohmann::json to_json(const Op_ptr &op);

 protected:
  void generate_circuit() const override;

 private:
  std::vector<SymPauliTensor> pauli_gadgets_;
  Transforms::PauliSynthStrat synth_strategy_;
  PauliPartitionStrat partition_strategy_;
  GraphColourMethod graph_colouring_;
  CXConfigType cx_configuration_;
};

}  // namespace tket
