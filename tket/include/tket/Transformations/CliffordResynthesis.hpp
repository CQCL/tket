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

#include <optional>

#include "Transform.hpp"

namespace tket {

namespace Transforms {

/**
 * Resynthesise all Clifford subcircuits and simplify using Clifford rules.
 *
 * @param transform optional user-provided resynthesis method to apply to all
 *   Clifford subcircuits (a function taking a Clifford circuit as an argument
 *   and returning an equivalent circuit); if not provided, a default
 *   resynthesis method is applied
 * @param allow_swaps whether the rewriting may introduce wire swaps (only
 *   relevant to the default resynthesis method used when the `transform`
 *   argument is not provided)
 * @return transform to perform Clifford resynthesis
 */
Transform clifford_resynthesis(
    std::optional<std::function<Circuit(const Circuit&)>> transform =
        std::nullopt,
    bool allow_swaps = true);

}  // namespace Transforms

}  // namespace tket
