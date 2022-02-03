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

#include "Utils/MatrixAnalysis.hpp"

namespace tket {
namespace tket_sim {
namespace internal {

struct GateNode {
  /** The triplets which make up the entries of the unitary matrix
   *  of a gate or box, if acting on qubits [0,1,2,...,k].
   */
  std::vector<TripletCd> triplets;

  /** The indices in the original top-level root circuit
   *  which this unitary acts upon.
   */
  std::vector<unsigned> qubit_indices;

  /** Lift the triplets to a full unitary matrix U acting on n qubits,
   *  then premultiply the given matrix by U.
   */
  void apply_full_unitary(
      Eigen::MatrixXcd& matr, unsigned full_number_of_qubits) const;
};

}  // namespace internal
}  // namespace tket_sim
}  // namespace tket
