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

#include "OpType/OpType.hpp"
#include "Utils/MatrixAnalysis.hpp"

namespace tket {
class Gate;

/** Functions to return the unitary matrix which a gate represents.
 *  All functions can throw GateUnitaryMatrixError exceptions.
 *  All functions will detect size mismatches (i.e., wrong number
 *  of parameters or qubits).
 */
struct GateUnitaryMatrix {
  /** The gate object knows how many qubits, and the (symbolic) parameters.
   *  Throws GateUnitaryMatrixError if not implemented, or upon error.
   *  Uses ILO-BE convention.
   */
  static Eigen::MatrixXcd get_unitary(const Gate& gate);

  /** Throws GateUnitaryMatrixError upon error, but does not check
   *  if the parameter values are all finite.
   *  Uses ILO-BE convention.
   */
  static Eigen::MatrixXcd get_unitary(
      OpType optype, unsigned number_of_qubits,
      const std::vector<double>& parameters);

  /** Return the unitary matrix of the gate, in sparse format, i.e.
   *  a collection of (i,j,z) tuples, meaning that U(i,j) = z.
   *  @param gate The gate. Throws GateUnitaryMatrixError upon error.
   *  @param abs_epsilon Used to decide whether an entry should be treated
   *          as zero. If std::abs(z) <= abs_epsilon then we treat z as
   *          zero exactly and so don't include it in the triplets.
   *  @return The triplets in the sparse representation.
   */
  static std::vector<TripletCd> get_unitary_triplets(
      const Gate& gate, double abs_epsilon = EPS);
};

}  // namespace tket
