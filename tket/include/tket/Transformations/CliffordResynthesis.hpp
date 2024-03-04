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

namespace tket {

namespace Transforms {

/**
 * Resynthesise all Clifford subcircuits and simplify using Clifford rules.
 *
 * @param allow_swaps whether to allow introduction of implicit wire swaps
 * @return transform to perform Clifford resynthesis
 */
Transform clifford_resynthesis(bool allow_swaps = true);

}  // namespace Transforms

}  // namespace tket
