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

/**
 * @file
 * @brief Utilities for the test programs
 */

#include <cstdlib>

#include "Circuit/Circuit.hpp"
#include "Utils/EigenConfig.hpp"

namespace tket {

/** Default tolerance for testing */
constexpr double ERR_EPS = 1e-10;

/** Random number in a given range */
static inline double frand(double fMin, double fMax) {
  double f = static_cast<double>(rand()) / RAND_MAX;
  return fMin + f * (fMax - fMin);
}

/** Compare the statevectors of two Circuits, assuming
 * they both start with |00...0> input
 */
bool test_statevector_comparison(
    const Circuit& circ1, const Circuit& circ2, bool projective = false);
bool test_unitary_comparison(
    const Circuit& circ1, const Circuit& circ2, bool projective = false);

/** Check that n_qubits() gives the right answer for all operations */
bool verify_n_qubits_for_ops(const Circuit& circ);

static inline bool test_equiv_expr(
    const Expr& e0, const Expr& e1, unsigned n = 2) {
  return equiv_expr(e0, e1, n, ERR_EPS);
}

static inline bool test_equiv_expr_c(const Expr& e0, const Expr& e1) {
  std::optional<Complex> c0 = eval_expr_c(e0);
  std::optional<Complex> c1 = eval_expr_c(e1);
  if (!c0 || !c1) return e0 == e1;
  Complex diff = c0.value() - c1.value();
  return std::abs(diff) < ERR_EPS;
}

static inline bool test_equiv_val(const Expr& e, double x, unsigned n = 2) {
  return equiv_val(e, x, n, ERR_EPS);
}

static inline bool test_equiv_0(const Expr& e, unsigned n = 2) {
  return equiv_0(SymEngine::expand(e), n, ERR_EPS);
}

typedef std::vector<unsigned> uvec;

/** Adds the same two-qubit op to the circuit multiple times,
 *  acting between the given sequence of qubit pairs.
 */
void add_2qb_gates(
    Circuit& circ, OpType op_type,
    const std::vector<std::pair<unsigned, unsigned>>& qubit_pairs);

/** Adds the same one-qubit op to the circuit multiple times,
 *  onto the given sequence of qubits.
 */
void add_1qb_gates(
    Circuit& circ, OpType op_type, const std::vector<unsigned>& qubits);

/**
 * remaps UnitID of `circ` to use provided `nodes`
 * `nodes` defaults to [Node(0), ..., Node(n)]
 */
void reassign_boundary(
    Circuit& circ_, std::optional<node_vector_t> nodes = std::nullopt);

/** REQUIREs that the list of commands from the circuit has types
 *  exactly matching the expected list, in that order.
 */
void check_command_types(
    const Circuit& circ, const std::vector<OpType>& expected_types);

/** Normally it is very naughty to test for exact equality with doubles,
 *  but there are occasionally cases where it is justified.
 *  However, the Eigen matrix equality operator == is giving trouble
 *  with some compilers, so replace it with this simple function.
 *  Intended only for dense Eigen matrices.
 */
template <class Matr1, class Matr2>
bool matrices_are_equal(const Matr1& mat1, const Matr2& mat2) {
  if (!(mat1.rows() == mat2.rows() && mat1.cols() == mat2.cols())) {
    return false;
  }
  for (unsigned jj = 0; jj < mat1.cols(); ++jj) {
    for (unsigned ii = 0; ii < mat1.rows(); ++ii) {
      if (mat1(ii, jj) != mat2(ii, jj)) {
        return false;
      }
    }
  }
  return true;
}

/**
 * Generate a random unitary matrix of a given size.
 *
 * Results are _not_ distributed uniformly in the Haar measure, so statistical
 * properties should not be relied upon. (However, all regions of positive
 * Haar measure have positive probability.)
 *
 * @param n matrix dimension
 * @param seed RNG seed
 *
 * @return unitary matrix generated from the seed
 */
Eigen::MatrixXcd random_unitary(unsigned n, int seed);

}  // namespace tket
