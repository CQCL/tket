// Copyright 2019-2023 Cambridge Quantum Computing
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
#include <stdexcept>
#include <string>

namespace tket {
namespace WeightedSubgraphMonomorphism {

/** All edge weights are integer valued (no doubles anywhere),
 * so calculations should be fully portable and identical
 * IF no overflows occur. If an overflow does occur,
 * this exception may be thrown.
 */
struct IntegerOverflow : public std::runtime_error {
  IntegerOverflow(const std::string& str) : std::runtime_error(str) {}
};

/** Initialisation of the data to start solving should be fast;
 * if it is too slow (e.g., exceeds the WHOLE problem timeout),
 * something is probably wrong; throw.
 */
struct InitialisationTimeout : public std::runtime_error {
  InitialisationTimeout(const std::string& str) : std::runtime_error(str) {}
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
