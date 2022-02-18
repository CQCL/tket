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

#include "Architecture/ArchitectureMapping.hpp"
#include "TokenSwapping/PartialTsaInterface.hpp"
#include "Utils/RNG.hpp"

namespace tket {
namespace tsa_internal {
namespace tests {

enum class RequiredTsaProgress { NONE, FULL, NONZERO };
enum class TokenOption {
  ALLOW_EMPTY_TOKEN_SWAP,
  DO_NOT_ALLOW_EMPTY_TOKEN_SWAP
};

/// Returns a summary string of the results, as well as doing the checks.
std::string run_tests(
    const ArchitectureMapping& arch_mapping,
    const std::vector<VertexMapping>& problems,
    RiverFlowPathFinder& path_finder, PartialTsaInterface& partial_tsa,
    RequiredTsaProgress progress,
    TokenOption token_option = TokenOption::DO_NOT_ALLOW_EMPTY_TOKEN_SWAP);

/// If no path finder is specified, will use the RiverFlowPathFinder
/// (which needs an RNG).
std::string run_tests(
    const ArchitectureMapping& arch_mapping,
    const std::vector<VertexMapping>& problems, RNG& rng,
    PartialTsaInterface& partial_tsa, RequiredTsaProgress progress,
    TokenOption token_option = TokenOption::DO_NOT_ALLOW_EMPTY_TOKEN_SWAP);

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
