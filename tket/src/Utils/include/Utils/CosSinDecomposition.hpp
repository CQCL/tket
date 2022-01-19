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

#include "EigenConfig.hpp"

namespace tket {

/**
 * Cosine-sine decomposition.
 *
 * The tuple (l0, l1, r0, r1, c, s) represents the decomposition
 *     [l0   ] [c -s] [r0   ]
 *     [   l1] [s  c] [   r1]
 * where l0, l1, r0 and r1 are unitaries of equal size, c and s are diagonal
 * matrices with non-negative entries, the diagonal entries of c are in
 * non-decreasimg order, and
 *     c^2 + s^2 = I
 */
typedef std::tuple<
    Eigen::MatrixXcd, Eigen::MatrixXcd, Eigen::MatrixXcd, Eigen::MatrixXcd,
    Eigen::MatrixXd, Eigen::MatrixXd>
    csd_t;

/**
 * Compute a cosine-sine decomposition of a unitary matrix.
 *
 * @param u unitary matrix to be decomposed
 * @return cosine-sine decomposition
 * @pre @p u has even dimensions
 */
csd_t CS_decomp(const Eigen::MatrixXcd &u);

}  // namespace tket
