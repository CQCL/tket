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

#include <map>
#include <tuple>
#include <utility>
#include <vector>

#include "Architecture/Architecture.hpp"

namespace tket {

/** This specifies desired source->target vertex mappings.
 *  Any nodes not occurring as a key might be moved by the algorithm.
 */
typedef std::map<Node, Node> NodeMapping;

/** Version 1.1, not too bad.
 *  @param architecture The raw object containing the graph.
 *  @param node_mapping The desired source->target node mapping.
 *  @return The required list of node pairs to swap.
 */
std::vector<std::pair<Node, Node>> get_swaps(
    const Architecture& architecture, const NodeMapping& node_mapping);

}  // namespace tket
