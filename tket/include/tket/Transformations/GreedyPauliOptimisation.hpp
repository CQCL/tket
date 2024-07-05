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
#include "tket/Architecture/Architecture.hpp"
#include "tket/Circuit/Circuit.hpp"

namespace tket {

namespace Transforms {

namespace GreedyPauliSimp {

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

/**
 * @brief Local Clifford
 *
 */
enum class LocalCliffordType {
  H,
  S,
  V,
};

/**
 * @brief Type for 2-qubit entangled Clifford gates
 *
 */
using TQE = std::tuple<TQEType, unsigned, unsigned>;

typedef std::vector<unsigned> pauli_letter_distances_t;

/**
 * @brief A Pauli exponential described by its commutation relations
 * with the rows in a reference Clifford tableau.
 * We store the commutation relations using an n-dimensional
 * vector with entries in {0,1,2,3}, where
 * 0: commute with ith Z row and ith X row
 * 1: commute with ith Z row and anti-commute with ith X row
 * 2: anti-commute with ith Z row and commute with ith X row
 * 3: anti-commute with ith Z row and anti-commute with ith X row
 * We call such vector a support vector
 */
class PauliExpNode {
 public:
  /**
   * @brief Construct a new PauliExpNode object.
   *
   * @param support_vec the support vector
   * @param theta the rotation angle in half-turns
   */
  PauliExpNode(std::vector<unsigned> support_vec, Expr theta);

  /**
   * @brief Number of TQEs required to reduce the weight to 1
   *
   * @return unsigned
   */
  unsigned tqe_cost() const { return tqe_cost_; }

  /**
   * @brief Number of TQEs would required to reduce the weight to 1
   * after the given TQE is applied
   *
   * @return unsigned
   */
  int tqe_cost_increase(const TQE& tqe) const;

  /**
   * @brief Weighted sum over number of Pauli letters
   * with no adjacent non-identities and distances
   * between letters.
   *
   * @return unsigned
   */
  double aa_tqe_cost_increase(
      const TQE& tqe, std::shared_ptr<Architecture> architecture,
      const std::map<unsigned, Node>& node_mapping,
      unsigned n_comparisons) const;

  std::vector<unsigned> get_updated_support(const TQE& tqe) const;

  bool updates_support(const TQE& tqe) const;

  pauli_letter_distances_t get_updated_distance(
      const TQE& tqe, std::shared_ptr<Architecture> architecture,
      const std::map<unsigned, Node>& node_mapping) const;

  /**
   * @brief For each pair of indices in the support_vec_, returns
   * a vector where value n at index d gives the number of pairs at
   * distance index on the architecure graph
   */
  pauli_letter_distances_t all_distances(
      const std::vector<unsigned>& support,
      std::shared_ptr<Architecture> architecture,
      const std::map<unsigned, Node>& node_mapping) const;

  /**
   * @brief Update the support vector with a TQE gate
   *
   * @param tqe
   */
  void update(const TQE& tqe);

  Expr theta() const { return theta_; };

  /**
   * @brief Return all possible TQE gates that will reduce the tqe cost by 1
   *
   * @return std::vector<std::tuple<TQEType, unsigned, unsigned>>
   */
  std::vector<TQE> reduction_tqes() const;

  /**
   * @brief return all TQE gates
   *
   * @return std::vector<std::pair<unsigned, unsigned>>
   */
  std::vector<TQE> reduction_tqes_all_letters(
      std::shared_ptr<Architecture> architecture,
      const std::map<unsigned, Node>& node_mapping) const;

  /**
   * @brief Return the index and value of the first support
   *
   * @return std::pair<unsigned, unsigned>
   */
  std::pair<unsigned, unsigned> first_support() const;

  void pad_support_vector(unsigned width);

  void print() const {
    for (auto p : support_vec_) {
      std::cout << p;
    }
    std::cout << std::endl;
  }

  unsigned support() const;

 private:
  std::vector<unsigned> support_vec_;
  Expr theta_;
  unsigned tqe_cost_;
};

/**
 * @brief Each row of a Clifford tableau consists a pair of anti-commuting
 * Pauli strings (p0,p1). Similar to the PauliExpNode, such pairs can be
 * alternatively described by their commutation relations with the rows in a
 * reference Clifford tableau. Let Xi and Zi be the ith X row and the ith Z row
 * in a reference Tableau T, then the commutation relation between (p0, p1) and
 * the ith row of T is defined by how p0, p1 commute with Xi and Zi. That's 4
 * bits of information. We store such information using an n-dimensional vector
 * with entries in {0,1,2,...,15}. The 4 bits from the most significant to the
 * least are: f(p0, Xi), f(p0,Zi), f(q,Xi), f(q,Zi) where f(p,q)==1 if p,q
 * anti-commute and 0 otherwise
 */
class TableauRowNode {
 public:
  /**
   * @brief Construct a new TableauRowNode object.
   *
   * @param support_vec the support vector
   */
  TableauRowNode(std::vector<unsigned> support_vec);

  /**
   * @brief Number of TQEs required to reduce the weight to 1
   *
   * @return unsigned
   */
  unsigned tqe_cost() const { return tqe_cost_; };

  /**
   * @brief Number of TQEs would required to reduce the weight to 1
   * after the given TQE is applied
   *
   * @return unsigned
   */
  int tqe_cost_increase(const TQE& tqe) const;

  /**
   * @brief Update the support vector with a TQE gate
   *
   * @param tqe
   */
  void update(const TQE& tqe);

  /**
   * @brief Return all possible TQE gates that will reduce the tqe cost
   *
   * @return std::vector<std::pair<unsigned, unsigned>>
   */
  std::vector<TQE> reduction_tqes() const;

  /**
   * @brief Return the index and value of the first support
   */
  std::pair<unsigned, unsigned> first_support() const;

 private:
  std::vector<unsigned> support_vec_;
  unsigned n_weaks_;
  unsigned n_strongs_;
  unsigned tqe_cost_;
};

/**
 * @brief The commutation relation between a TableauRowNode (p0,p1) and the ith
 * row of the reference Tableau can be further classified as Strong, Weak or
 * No-support.
 */
enum class SupportType : unsigned {
  Strong,
  Weak,
  No,
};

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
    const Circuit& circ, double discount_rate = 0.7, double depth_weight = 0.3,
    std::optional<std::shared_ptr<Architecture>> architecture = std::nullopt);
}  // namespace GreedyPauliSimp

Transform greedy_pauli_optimisation(
    double discount_rate = 0.7, double depth_weight = 0.3);

Transform aa_greedy_pauli_optimisation(
    std::shared_ptr<Architecture> architecture, double discount_rate = 0.7,
    double depth_weight = 0.3);

}  // namespace Transforms

}  // namespace tket
