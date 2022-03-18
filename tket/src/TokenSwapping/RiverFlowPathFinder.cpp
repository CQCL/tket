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

#include "RiverFlowPathFinder.hpp"

#include <sstream>
#include <stdexcept>

#include "SwapFunctions.hpp"
#include "Utils/Assert.hpp"

using std::vector;

namespace tket {
namespace tsa_internal {

struct RiverFlowPathFinder::Impl {
  DistancesInterface& distances_calculator;
  NeighboursInterface& neighbours_calculator;
  RNG& rng;

  typedef std::uint64_t EdgeCount;

  /** The key is an undirected edge; the value is the number of times
   *  that edge was already used in any requested path.
   *  (So, we favour flows in both directions).
   *  Overflow is basically impossible, but even if it did occur,
   *  it would not invalidate the results (it just means that some
   *  paths might change more than expected).
   */
  std::map<Swap, EdgeCount> edge_counts;

  struct ArrowData {
    size_t end_vertex;
    EdgeCount count;
  };

  /** A work vector. When we are trying to expand a path by one step,
   *  we need to list all those steps which are valid, i.e. reduce
   *  the distance to the target by one.
   */
  vector<ArrowData> candidate_moves;

  /// A work vector, will be built up
  vector<size_t> path;

  Impl(
      DistancesInterface& distances_interface,
      NeighboursInterface& neighbours_interface, RNG& random_generator)
      : distances_calculator(distances_interface),
        neighbours_calculator(neighbours_interface),
        rng(random_generator) {}

  void reset();

  /// Increases nonempty "path" towards the target vertex.
  void grow_path(size_t target_vertex, size_t required_path_size);

  /// Once "path" has been filled, update the counts (so that future paths
  /// through similar vertices are more likely to overlap).
  void update_data_with_path();
};

void RiverFlowPathFinder::Impl::reset() {
  for (auto& entry : edge_counts) {
    entry.second = 0;
  }
  rng.set_seed();
}

void RiverFlowPathFinder::Impl::grow_path(
    size_t target_vertex, size_t required_path_size) {
  TKET_ASSERT(path.size() < required_path_size);
  TKET_ASSERT(!path.empty());

  // We don't yet know how to move on, so we must choose a neighbour.
  // All candidates will have the same edge count.
  candidate_moves.clear();

  const auto remaining_distance = required_path_size - path.size();
  const auto& neighbours = neighbours_calculator(path.back());
  distances_calculator.register_neighbours(path.back(), neighbours);

  for (size_t neighbour : neighbours) {
    const auto neighbour_distance_to_target =
        distances_calculator(neighbour, target_vertex);

    if (neighbour_distance_to_target == remaining_distance - 1) {
      // Notice that nonexistent entries will be automatically set
      // to have count 0, by the C++ standard.
      const auto edge_count = edge_counts[get_swap(path.back(), neighbour)];
      if (!candidate_moves.empty()) {
        // We'll only add candidates with the same count or higher.
        if (candidate_moves[0].count > edge_count) {
          continue;
        }
        if (candidate_moves[0].count < edge_count) {
          candidate_moves.clear();
        }
      }
      candidate_moves.emplace_back();
      candidate_moves.back().end_vertex = neighbour;
      candidate_moves.back().count = edge_count;
      continue;
    }
    // GCOVR_EXCL_START
    TKET_ASSERT(
        neighbour_distance_to_target == remaining_distance ||
        neighbour_distance_to_target == remaining_distance + 1 ||
        AssertMessage() << "d(v_" << path.back() << ", v_" << target_vertex
                        << ")=" << remaining_distance << ". But v_"
                        << path.back() << " has neighbour v_" << neighbour
                        << ", at distance " << neighbour_distance_to_target
                        << " to the target v_" << target_vertex);
    // GCOVR_EXCL_STOP
  }
  // GCOVR_EXCL_START
  TKET_ASSERT(
      !candidate_moves.empty() ||
      AssertMessage() << "No neighbours of v_" << path.back()
                      << " at correct distance " << remaining_distance - 1
                      << " to target vertex v_" << target_vertex);
  // GCOVR_EXCL_STOP

  const auto& choice = rng.get_element(candidate_moves);
  path.push_back(choice.end_vertex);
}

void RiverFlowPathFinder::Impl::update_data_with_path() {
  for (size_t ii = 1; ii < path.size(); ++ii) {
    // Nonexistent counts automatically set to 0 initially
    ++edge_counts[get_swap(path[ii - 1], path[ii])];
  }
  distances_calculator.register_shortest_path(path);
}

RiverFlowPathFinder::RiverFlowPathFinder(
    DistancesInterface& distances_interface,
    NeighboursInterface& neighbours_interface, RNG& rng)
    : m_pimpl(std::make_unique<Impl>(
          distances_interface, neighbours_interface, rng)) {}

RiverFlowPathFinder::~RiverFlowPathFinder() {}

void RiverFlowPathFinder::reset() { m_pimpl->reset(); }

const vector<size_t>& RiverFlowPathFinder::operator()(
    size_t vertex1, size_t vertex2) {
  m_pimpl->path.clear();
  m_pimpl->path.push_back(vertex1);
  if (vertex1 == vertex2) {
    return m_pimpl->path;
  }

  // We must build up the path.
  // The number of vertices including the source and target.
  const size_t final_path_size =
      1 + m_pimpl->distances_calculator(vertex1, vertex2);

  for (size_t infinite_loop_guard = 10 * final_path_size;
       infinite_loop_guard != 0; --infinite_loop_guard) {
    m_pimpl->grow_path(vertex2, final_path_size);
    if (m_pimpl->path.size() == final_path_size) {
      TKET_ASSERT(m_pimpl->path.back() == vertex2);
      m_pimpl->update_data_with_path();
      return m_pimpl->path;
    }
  }
  throw std::runtime_error("get path - dropped out of loop");
}

void RiverFlowPathFinder::register_edge(size_t vertex1, size_t vertex2) {
  // Automatically zero if the edge doesn't exist.
  auto& edge_count = m_pimpl->edge_counts[get_swap(vertex1, vertex2)];
  ++edge_count;
}

}  // namespace tsa_internal
}  // namespace tket
