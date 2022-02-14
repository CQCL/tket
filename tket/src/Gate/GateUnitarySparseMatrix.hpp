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
class Gate;
namespace internal {

/** For getting sparse unitary matrices directly for specific gates,
 *  without constructing the dense matrices.
 */
struct GateUnitarySparseMatrix {
  /** If the gate is an unknown type, returns an empty vector.
   *  (However, that only means that there is no specific sparse function.
   *  It may still be possible to get a dense unitary matrix
   *  from other functions).
   *
   *  Return the unitary matrix of the gate in sparse format, i.e.
   *  a collection of (i,j,z) tuples, meaning that U(i,j) = z.
   *
   *  @param gate unitary quantum gate
   *  @param abs_epsilon Used to decide whether an entry should be treated
   *          as zero. If std::abs(z) <= abs_epsilon then we treat z as
   *          zero exactly and so don't include it in the triplets.
   *  @return The triplets in the sparse representation.
   */
  static std::vector<TripletCd> get_unitary_triplets(
      const Gate& gate, double abs_epsilon = EPS);
};

}  // namespace internal
}  // namespace tket
