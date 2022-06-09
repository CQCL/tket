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

#include <array>
#include <vector>

#include "EigenConfig.hpp"
#include "UnitID.hpp"
#include "Utils/Constants.hpp"

namespace tket {

typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> MatrixXb;
typedef Eigen::Matrix<bool, Eigen::Dynamic, 1> VectorXb;
typedef Eigen::Matrix<bool, 2, 1> Vector2b;
// Eigen::Matrix8cd doesn't exist.
typedef Eigen::Matrix<std::complex<double>, 8, 8> Matrix8cd;

/**
 * Test matrix for unitarity.
 *
 * The tolerance refers to the l2 (Frobenius) norm of I - U*U.
 *
 * @param U matrix
 * @param tol tolerance
 *
 * @retval whether the matrix is unitary to within some tolerance
 */
bool is_unitary(const Eigen::MatrixXcd &U, double tol = EPS);
bool is_projector(const Eigen::MatrixXcd &P, double tol = EPS);

/**
 * Lift a permutation of \f$ [0,n) \f$ to a permutation of \f$ [0,2^n) \f$ as a
 * matrix.
 *
 * The result is a \f$ 2^n \times 2^n \f$ permutation matrix representing the
 * derived permutation on the powerset of \f$ [0,n) \f$ using big-endian
 * encoding, i.e. a subset \f$ S \subseteq [0,n) \f$ is represented by
 * \f$ f(S) = \sum_{i \in S} 2^{n-1-i} \f$ and the \f$ (f(S),f(T)) \f$ entry of
 * the matrix is 1 iff \f$ T = p(S) \f$.
 *
 * @param p permutation of \f$ [0,n) \f$
 *
 * @return \f$ 2^n \times 2^n \f$ permutation matrix
 */
Eigen::PermutationMatrix<Eigen::Dynamic> lift_perm(
    const std::map<unsigned, unsigned> &p);

// Convert between ILO-BE and DLO-BE conventions
Eigen::Matrix4cd reverse_indexing(const Eigen::Matrix4cd &m);
Matrix8cd reverse_indexing(const Matrix8cd &m);
Eigen::MatrixXcd reverse_indexing(const Eigen::MatrixXcd &m);
Eigen::VectorXcd reverse_indexing(const Eigen::VectorXcd &v);

Eigen::MatrixXcd apply_qubit_permutation(
    const Eigen::MatrixXcd &m, const qubit_map_t &perm);
Eigen::VectorXcd apply_qubit_permutation(
    const Eigen::VectorXcd &v, const qubit_map_t &perm);

std::pair<MatrixXb, MatrixXb> binary_LLT_decomposition(const MatrixXb &a);
std::vector<std::pair<unsigned, unsigned>> gaussian_elimination_col_ops(
    const MatrixXb &a, unsigned blocksize = 6);
std::vector<std::pair<unsigned, unsigned>> gaussian_elimination_row_ops(
    const MatrixXb &a, unsigned blocksize = 6);

/**
 * Performs KAK decomposition.
 *
 * Given a unitary \f$ X \f$, returns
 * \f$ (K_1, (k_\mathrm{XX},k_\mathrm{YY},k_\mathrm{ZZ}, K_2) \f$
 * such that
 * \f$ X = K_1 \cdot \exp(-\frac12 i \pi \sum_{w \in \{\mathrm{XX}, \mathrm{YY},
 * \mathrm{ZZ}\}} k_w \sigma_w) \cdot K_2 \f$.
 * The \f$ k_w \f$ are called information content and partition
 * \f$ \mathrm{SU}(4) \f$ into equivalence classes modulo local transformations.
 *
 * See arXiv quant-ph/0507171 for details
 *
 * @param X unitary 4x4 matrix
 *
 * @return KAK decomposition
 */
std::tuple<Eigen::Matrix4cd, std::array<double, 3>, Eigen::Matrix4cd>
get_information_content(const Eigen::Matrix4cd &X);

// given a 4x4 unitary matrix (ILO-BE), returns two 2x2 unitaries that
// approximately make the input by kronecker product
std::pair<Eigen::Matrix2cd, Eigen::Matrix2cd> kronecker_decomposition(
    Eigen::Matrix4cd &U);

/** Returns 2^n, but throws if n is too large (causes overflow).
 *  (Warning: simply taking 1<<n and testing against zero
 *  can give incorrect results with large n).
 */
unsigned get_matrix_size(unsigned number_of_qubits);

/** We have a matrix size, which should be 2^n. Return n, or throw if invalid.
 */
unsigned get_number_of_qubits(unsigned matrix_size);

/** It is sometimes more convenient to deal with Triplets directly,
 *  rather than sparse matrices.
 */
typedef Eigen::Triplet<std::complex<double>> TripletCd;

typedef Eigen::SparseMatrix<std::complex<double>> SparseMatrixXcd;

SparseMatrixXcd get_sparse_matrix(
    const std::vector<TripletCd> &triplets, unsigned rows, unsigned cols);

SparseMatrixXcd get_sparse_square_matrix(
    const std::vector<TripletCd> &triplets, unsigned rows);

/** abs_epsilon is used to decide if a near-zero entry should be
 *  set to zero exactly. Thus, if std::abs(z) <= abs_epsilon, then
 *  z is treated as zero.
 */
std::vector<TripletCd> get_triplets(
    const SparseMatrixXcd &matr, double abs_epsilon = EPS);

/** Convert a matrix M into a list of tuples (i,j,z), meaning that M(i,j)=z,
 *  used in sparse representations of M.
 *  @param matr The dense matrix.
 *  @param abs_epsilon If std::abs(z) <= abs_epsilon, then z is so small
 *              that it is treated as zero exactly and hence not added to
 *              the list of triplets.
 *  @return A vector of (i,j,z) tuples, meaning that M(i,j)=z.
 */
std::vector<TripletCd> get_triplets(
    const Eigen::MatrixXcd &matr, double abs_epsilon = EPS);

// TODO: move get_CX_fidelity to Decomposition.cpp when no longer needed here
/** Given TK2 angles, computes the fidelity that can be achieved using
 *  nb_cx CX gates.
 *
 *  @param k The TK2 angles.
 *  @param nb_cx The number of CX gates to be used for decomposition.
 *  @return The fidelity.
 */
double get_CX_fidelity(const std::array<double, 3> &k, unsigned nb_cx);

/** Similarity measure of TK2(a, b, c) to SU(4) identity
 *
 * This computes the fidelity between TK2(a, b, c) and the 2-qubit identity.
 *
 * a, b and c must be in the Weyl chamber, i.e. 1/2 >= a >= b >= |c|.
 *
 * This is computed using the formula
 *               Fidᵤ = (4 + Tr(Id₄U)) / 20 = (4 + Tr(U)) / 20
 * where U is the 4x4 matrix of TK2(a, b, c).
 *
 * Tr(U) can in turn be computed as
 *               Tr(U) = 4cos(a)cos(b)cos(c) − 4i sin(a)sin(b)sin(c).
 *
 * These are formulas B3 and B5 of https://arxiv.org/pdf/1811.12926.pdf. Refer
 * to that paper for more details.
 *
 *  @param a The XX interaction angle
 *  @param b The YY interaction angle
 *  @param c The ZZ interaction angle
 */
double trace_fidelity(double a, double b, double c);

}  // namespace tket
