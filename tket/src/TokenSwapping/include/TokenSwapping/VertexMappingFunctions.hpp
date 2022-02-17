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
#include <optional>
#include <utility>

#include "SwapFunctions.hpp"

namespace tket {

/// The desired result of swapping is to move a token on each "key"
/// vertex to the "value" vertex.
typedef std::map<size_t, size_t> VertexMapping;

/** Are all tokens on their target vertices?
 *  @param vertex_mapping The desired mapping.
 *  @return Whether all tokens are on their target vertices.
 */
bool all_tokens_home(const VertexMapping& vertex_mapping);

/** Does nothing, except throwing if the mapping is invalid.
 *  @param vertex_mapping The desired mapping, to be checked.
 */
void check_mapping(const VertexMapping& vertex_mapping);

/** When you've already got another expendable VertexMapping object,
 *  it saves time to reuse instead of constructing a new one.
 *  @param vertex_mapping The desired mapping, to be checked.
 *  @param work_mapping A disposable object, will be overwritten.
 */
void check_mapping(
    const VertexMapping& vertex_mapping, VertexMapping& work_mapping);

/** We have a path [v(1), v(2), v(3), ..., v(N)].
 *  Calculate individual swaps along this path (i.e., using only
 *  Swap(v(i), v(i+1)) which we know are valid), which would swap the tokens
 *  (if any) on v(1), v(N), and perform the swaps.
 *  Only append nonempty swaps (i.e., where at least one token is moved).
 *  @param path The path (must be an actual possible path), whose start
 *      and end vertices are to be swapped (with all other vertices)
 *  @param vertex_mapping The source to target mapping, which will be updated.
 *  @param swap_list The list of swaps, which will be updated.
 */
void append_swaps_to_interchange_path_ends(
    const std::vector<size_t>& path, VertexMapping& vertex_mapping,
    SwapList& swap_list);

/** Given a source->target vertex mapping and a TARGET vertex, find the
 * corresponding source vertex. If the given target vertex does not appear in
 * the map, create it as a new fixed vertex, i.e. map[v] = v for the given
 * target vertex v.
 * @param source_to_target_map A source->target vertex mapping.
 * @param target_vertex A target vertex, to find in the map.
 * @return The source vertex corresponding to the target (possibly newly created
 * if the target was not present).
 */
size_t get_source_vertex(
    VertexMapping& source_to_target_map, size_t target_vertex);

/** We currently have a source->target mapping. Perform the vertex swap,
 * but if any vertex in the swap is not present in the map, add it to the map as
 * a new source vertex.
 * Note that, since we DON'T have a target->source map, we have to do an O(N)
 * search to find all target vertices.
 * @param source_to_target_map The map to update with the swap.
 * @param swap The swap to perform.
 */
void add_swap(VertexMapping& source_to_target_map, const Swap& swap);

}  // namespace tket
