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

#include "PauliGraph/PauliGraph.hpp"

namespace tket {

struct cmp_tensors {
  bool operator()(
      const QubitPauliTensor &qps1, const QubitPauliTensor &qps2) const {
    return (qps1.string < qps2.string);
  }
};

/**
 * QubitOperator, defined to be useful for diagonalisation and
 * partitioning.
 */
typedef std::map<QubitPauliTensor, Expr, cmp_tensors> QubitOperator;

void insert_into_gadget_map(
    QubitOperator &gadget_map, const PauliGadgetProperties &pgp);

}  // namespace tket
