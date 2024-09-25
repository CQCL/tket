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

#include "Transform.hpp"
#include "tket/Circuit/Circuit.hpp"

namespace tket {

namespace Transforms {

namespace GreedyPauliSimp {

class GreedyPauliSimpError : public std::logic_error {
 public:
  explicit GreedyPauliSimpError(const std::string& message)
      : std::logic_error(message) {}
};

/**
 * @brief Types of 2-qubit entangled Clifford gates
 *
 */
enum class TQEType : unsigned {
  XX,
  XY,
  XZ,
  YX,
  YY,
  YZ,
  ZX,
  ZY,
  ZZ,
};

enum class PauliNodeType {
  // Pauli rotation
  Rotation,
  // Defines how a Pauli X and a Pauli Z on the same qubit
  // get propagated from right to left through a Clifford operator.
  Propagation,
};

/**
 * @brief The type of a pair of Pauli letters defined by
    their commutation relation
 *
 */
enum class CommuteType : unsigned {
  // Both are (I)dentity
  I,
  // (A)nti-commute
  A,
  // (C)ommute and not both identity
  C,
};

/**
 * @brief Type for 2-qubit entangled Clifford gates
 *
 */
using TQE = std::tuple<TQEType, unsigned, unsigned>;

class PauliNode {
 public:
  virtual PauliNodeType get_type() const = 0;
  virtual unsigned tqe_cost() const = 0;
  virtual int tqe_cost_increase(const TQE& tqe) const = 0;
  virtual void update(const TQE& tqe) = 0;
  virtual void update(const OpType& sq_cliff, const unsigned& a);
  virtual void swap(const unsigned& a, const unsigned& b);
  virtual std::vector<TQE> reduction_tqes() const = 0;
  virtual ~PauliNode();
};

typedef std::shared_ptr<PauliNode> PauliNode_ptr;

/**
 * @brief A node defined by a single Pauli string
 */
class SingleNode : public PauliNode {
 public:
  /**
   * @brief Construct a new SinglePauliNode object.
   *
   * @param string the Pauli string
   */
  SingleNode(std::vector<Pauli> string, bool sign);

  /**
   * @brief Number of TQEs required to reduce the weight to 1
   *
   * @return unsigned
   */
  unsigned tqe_cost() const override;

  /**
   * @brief Number of TQEs would required to reduce the weight to 1
   * after the given TQE is applied
   *
   * @return unsigned
   */
  int tqe_cost_increase(const TQE& tqe) const override;

  /**
   * @brief Update the Pauli string with a TQE gate
   *
   * @param tqe
   */
  void update(const TQE& tqe) override;

  /**
   * @brief Return all possible TQE gates that will reduce the tqe cost by 1
   *
   * @return std::vector<std::tuple<TQEType, unsigned, unsigned>>
   */
  std::vector<TQE> reduction_tqes() const override;

  /**
   * @brief Return the index and value of the first non-identity
   *
   * @return std::pair<unsigned, Pauli>
   */
  std::pair<unsigned, Pauli> first_support() const;

  bool sign() const { return sign_; };

 protected:
  std::vector<Pauli> string_;
  bool sign_;
  // extra cached data used by greedy synthesis
  unsigned weight_;
};

/**
 * @brief Node consists of a pair of anti-commuting Pauli strings
 */
class ACPairNode : public PauliNode {
 public:
  /**
   * @brief Construct a new ACPairNode object
   *
   * @param z_propagation
   * @param x_propagation
   * @param z_sign
   * @param x_sign
   */
  ACPairNode(
      std::vector<Pauli> z_propagation, std::vector<Pauli> x_propagation,
      bool z_sign, bool x_sign);

  /**
   * @brief Number of TQEs required to reduce the weight to 1
   *
   * @return unsigned
   */
  unsigned tqe_cost() const override;

  /**
   * @brief Number of TQEs would required to reduce the weight to 1
   * after the given TQE is applied
   *
   * @return unsigned
   */
  int tqe_cost_increase(const TQE& tqe) const override;

  /**
   * @brief Update the support vector with a TQE gate
   *
   * @param tqe
   */
  void update(const TQE& tqe) override;

  /**
   * @brief Update the support vector with a single-qubit Clifford gate
   *
   * @param tqe
   */
  void update(const OpType& sq_cliff, const unsigned& a) override;

  /**
   * @brief Update the support vector with a SWAP gate
   *
   * @param tqe
   */
  void swap(const unsigned& a, const unsigned& b) override;

  /**
   * @brief Return all possible TQE gates that will reduce the tqe cost
   *
   * @return std::vector<TQE>
   */
  std::vector<TQE> reduction_tqes() const override;

  /**
   * @brief Return the index and value of the first anti-commute entry
   */
  std::tuple<unsigned, Pauli, Pauli> first_support() const;

  bool z_sign() const { return z_sign_; };

  bool x_sign() const { return x_sign_; };

  std::vector<Pauli> z_propagation() const { return z_propagation_; };

  std::vector<Pauli> x_propagation() const { return x_propagation_; };

 protected:
  std::vector<Pauli> z_propagation_;
  std::vector<Pauli> x_propagation_;
  bool z_sign_;
  bool x_sign_;
  // extra cached data used by greedy synthesis
  std::vector<CommuteType> commute_type_vec_;
  unsigned n_commute_entries_;
  unsigned n_anti_commute_entries_;
};

/**
 * @brief A Pauli exponential defined by a padded Pauli string
 * and a rotation angle
 */
class PauliRotation : public SingleNode {
 public:
  /**
   * @brief Construct a new PauliRotation object.
   *
   * @param string the Pauli string
   * @param theta the rotation angle in half-turns
   */
  PauliRotation(std::vector<Pauli> string, Expr theta);

  PauliNodeType get_type() const { return PauliNodeType::Rotation; };

  Expr theta() const { return theta_; };

 private:
  Expr theta_;
};

/**
 * @brief Defines how a Pauli X and a Pauli Z on the same qubit
 * get propagated from right to left through a Clifford operator.
 * A n-qubit Clifford operator is completely defined by n such propagations
 * with one on each qubit. A PauliPropagation also corresponds to a row in
 * a Clifford tableau
 */
class PauliPropagation : public ACPairNode {
 public:
  /**
   * @brief Construct a new PauliPropagation object
   *
   * @param z_propagation
   * @param x_propagation
   * @param z_sign
   * @param x_sign
   * @param qubit_index
   */
  PauliPropagation(
      std::vector<Pauli> z_propagation, std::vector<Pauli> x_propagation,
      bool z_sign, bool x_sign, unsigned qubit_index);

  PauliNodeType get_type() const { return PauliNodeType::Propagation; };

  unsigned qubit_index() const { return qubit_index_; };

 private:
  unsigned qubit_index_;
};

/**
 * @brief Convert a unordered set of SymPauliTensor into a set of PauliRotations
 * followed by a set of PauliPropagations
 *
 * @param unordered_set
 * @return std::tuple<std::vector<PauliNode_ptr>, std::vector<PauliNode_ptr>>
 */
std::tuple<std::vector<PauliNode_ptr>, std::vector<PauliNode_ptr>>
gpg_from_unordered_set(const std::vector<SymPauliTensor>& unordered_set);

/**
 * @brief Transforms a circuit of PauliExpBoxes followed by Clifford gates and
 * end-of-circuit measurements into a sequence of PauliRotations,
 * PauliPropagations, and an end-of-circuit measurement circuit.
 *
 * The first returned object is an empty circuit onto which the synthesized
 * circuit will be built. The final element of the returned tuple is a map from
 * integers to the original UnitIDs, allowing the final circuit to use the
 * orignal UnitIDs.
 *
 * @param circ
 * @return std::tuple<
 * Circuit,
 * std::vector<std::vector<PauliNode_ptr>>, std::vector<PauliNode_ptr>,
 * Circuit, unit_map_t>
 */
std::tuple<
    Circuit, std::vector<std::vector<PauliNode_ptr>>,
    std::vector<PauliNode_ptr>, Circuit, unit_map_t>
gpg_from_circuit(const Circuit& circ);

/**
 * @brief Given a circuit consists of PauliExpBoxes followed by clifford gates,
 * and end-of-circuit measurements, implement the PauliExpBoxes and the final
 * clifford subcircuit by applying Clifford gates and single qubit rotations in
 * a greedy fashion.
 *
 * @param circ
 * @param discount_rate
 * @param depth_weight
 * @return Circuit
 */
Circuit greedy_pauli_graph_synthesis(
    const Circuit& circ, double discount_rate = 0.7, double depth_weight = 0.3);

/**
 * @brief Synthesise a set of unordered Pauli exponentials by applying Clifford
 * gates and single qubit rotations in a greedy fashion. Assume all Pauli
 * strings have the same length.
 *
 * @param unordered_set
 * @param depth_weight
 * @return Circuit
 */
Circuit greedy_pauli_set_synthesis(
    const std::vector<SymPauliTensor>& unordered_set,
    double depth_weight = 0.3);

}  // namespace GreedyPauliSimp

Transform greedy_pauli_optimisation(
    double discount_rate = 0.7, double depth_weight = 0.3);

}  // namespace Transforms

}  // namespace tket
