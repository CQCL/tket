// Copyright Quantinuum
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

#include "WeightedSquareGrid.hpp"

#include <catch2/catch_test_macros.hpp>
#include <numeric>
#include <tkrng/RNG.hpp>
#include <tkwsm/Common/GeneralUtils.hpp>

namespace tket {
namespace WeightedSubgraphMonomorphism {
namespace InitialPlacement {
namespace tests {

static unsigned get_width_from_number_of_edges(unsigned number_of_edges) {
  // We require width=height, where "width" means the "geometric" horizontal
  // distance across; thus each longest horizontal path actually has width+1
  // vertices.
  // Thus there are 2w(w+1) edges in total.
  // Solve the quadratic extremely crudely!
  REQUIRE(number_of_edges > 10);
  REQUIRE(number_of_edges < 100000);
  for (unsigned width = 1; width < 1000; ++width) {
    const unsigned recalc_number_of_edges = 2 * width * (width + 1);
    if (recalc_number_of_edges == number_of_edges) {
      return width;
    }
    REQUIRE(recalc_number_of_edges < number_of_edges);
  }
  // Will never be reached - but compilers probably won't work that out.
  return 0;
}

static std::vector<std::uint64_t> get_vertical_horizontal_patterns(
    unsigned number_to_try) {
  number_to_try = std::max<unsigned>(number_to_try, 10);

  std::vector<std::uint64_t> result(2);
  result.reserve(number_to_try);
  // Start off horizontal, until we can't any more.
  result[0] = 0;

  // Start off vertical, until we can't any more.
  set_maximum(result[1]);

  RNG rng;
  while (result.size() < number_to_try) {
    result.push_back(rng());
  }
  return result;
}

WeightedSquareGrid::WeightedSquareGrid(
    std::vector<WeightWSM> weights, unsigned number_of_primitive_gates_in_swap)
    : m_weights(std::move(weights)),
      m_width(get_width_from_number_of_edges(m_weights.size())),
      m_number_of_vertices((m_width + 1) * (m_width + 1)),
      m_vertical_horizontal_patterns(
          get_vertical_horizontal_patterns(10 * m_width)) {
  m_number_of_primitive_gates_in_swap = number_of_primitive_gates_in_swap;
  REQUIRE(m_number_of_primitive_gates_in_swap >= 1);
  REQUIRE(m_number_of_primitive_gates_in_swap <= 100);
  REQUIRE(m_width < 64);
  REQUIRE(get_vertex(0, 0) == 0);
  REQUIRE(get_vertex(m_width, 0) == m_width);
  REQUIRE(get_vertex(m_width, m_width) == m_number_of_vertices - 1);
  REQUIRE(get_vertex(0, m_width) + m_width == get_vertex(m_width, m_width));
}

GraphEdgeWeights WeightedSquareGrid::get_graph_data() const {
  GraphEdgeWeights result;
  for (unsigned xx = 0; xx <= m_width; ++xx) {
    for (unsigned yy = 0; yy <= m_width; ++yy) {
      const VertexWSM this_v = get_vertex(xx, yy);
      if (xx < m_width) {
        const VertexWSM horiz_next_v = get_vertex(xx + 1, yy);
        result[get_edge(this_v, horiz_next_v)] =
            m_weights.at(get_horizontal_weight_index(xx, yy));
      }
      if (yy < m_width) {
        const VertexWSM vert_next_v = get_vertex(xx, yy + 1);
        result[get_edge(this_v, vert_next_v)] =
            m_weights.at(get_vertical_weight_index(xx, yy));
      }
    }
  }
  REQUIRE(result.size() == m_weights.size());
  WeightWSM total_weight = 0;
  for (const auto& entry : result) {
    total_weight += entry.second;
  }
  REQUIRE(
      total_weight ==
      std::accumulate(m_weights.cbegin(), m_weights.cend(), WeightWSM(0)));
  return result;
}

VertexWSM WeightedSquareGrid::get_vertex(unsigned x, unsigned y) const {
  REQUIRE(x <= m_width);
  REQUIRE(y <= m_width);
  // The longest horizontal path for y=0 is [0 1 2 ... w], etc.
  return x + y * (m_width + 1);
}

unsigned WeightedSquareGrid::get_horizontal_weight_index(
    unsigned x, unsigned y) const {
  REQUIRE(x < m_width);
  REQUIRE(y <= m_width);
  return x + y * m_width;
}

unsigned WeightedSquareGrid::get_vertical_weight_index(
    unsigned x, unsigned y) const {
  REQUIRE(x <= m_width);
  REQUIRE(y < m_width);
  return get_horizontal_weight_index(m_width - 1, m_width) + 1 + y +
         x * m_width;
}

const PlacementCostModelInterface::Path& WeightedSquareGrid::get_path_to_use(
    VertexWSM vertex1, VertexWSM vertex2) const {
  REQUIRE(vertex1 != vertex2);
  const VertexWSM start_v = std::min(vertex1, vertex2);
  const VertexWSM end_v = std::max(vertex1, vertex2);

  const auto key = std::make_pair(start_v, end_v);
  const auto reversed_key = std::make_pair(end_v, start_v);

  // Insert empty vector if not already present.
  m_paths[key];

  auto& reversed_path = m_paths[reversed_key];

  // Necessary to do this here, not above, to ensure a valid reference!
  auto& path = m_paths[key];

  REQUIRE(path.size() == reversed_path.size());
  if (path.empty()) {
    // The path hasn't yet been calculated.
    // Use reversed_path as temporary storage.
    fill_path(path, reversed_path, start_v, end_v);
    reversed_path = path;
    reverse_path(reversed_path);
  }
  if (start_v == vertex1) {
    return path;
  }
  return reversed_path;
}

// This is all very crude; but it's a test, simplicity more important than
// efficiency!
void WeightedSquareGrid::fill_path(
    Path& path, Path& path_work_vector, VertexWSM vertex1,
    VertexWSM vertex2) const {
  const auto start_xy = get_xy(vertex1);
  const auto end_xy = get_xy(vertex2);
  const unsigned number_of_vertices_in_path =
      1 + std::abs(int(start_xy.first) - int(end_xy.first)) +
      std::abs(int(start_xy.second) - int(end_xy.second));

  WeightWSM best_weight;
  set_maximum(best_weight);
  const std::uint64_t mask_bit = 1;
  for (std::uint64_t pattern_uint : m_vertical_horizontal_patterns) {
    path_work_vector.clear();
    path_work_vector.emplace_back(vertex1, 0);
    WeightWSM total_weight = 0;

    // Incrementally add a vertex to the path.
    while (path_work_vector.back().first != vertex2 &&
           path_work_vector.size() < number_of_vertices_in_path) {
      bool go_horizontal = (mask_bit & pattern_uint) == 0;
      pattern_uint >>= 1;
      const auto current_xy = get_xy(path_work_vector.back().first);

      if (current_xy.first == end_xy.first) {
        go_horizontal = false;
      }
      if (current_xy.second == end_xy.second) {
        go_horizontal = true;
      }
      if (go_horizontal) {
        if (current_xy.first < end_xy.first) {
          const auto next_v =
              get_vertex(current_xy.first + 1, current_xy.second);
          path_work_vector.emplace_back(
              next_v, m_weights.at(get_horizontal_weight_index(
                          current_xy.first, current_xy.second)));
        } else {
          REQUIRE(current_xy.first > end_xy.first);
          const auto next_v =
              get_vertex(current_xy.first - 1, current_xy.second);
          path_work_vector.emplace_back(
              next_v, m_weights.at(get_horizontal_weight_index(
                          current_xy.first - 1, current_xy.second)));
        }
      } else {
        // Vertical motion.
        if (current_xy.second < end_xy.second) {
          const auto next_v =
              get_vertex(current_xy.first, current_xy.second + 1);
          path_work_vector.emplace_back(
              next_v, m_weights.at(get_vertical_weight_index(
                          current_xy.first, current_xy.second)));
        } else {
          REQUIRE(current_xy.second > end_xy.second);
          const auto next_v =
              get_vertex(current_xy.first, current_xy.second - 1);
          path_work_vector.emplace_back(
              next_v, m_weights.at(get_vertical_weight_index(
                          current_xy.first, current_xy.second - 1)));
        }
      }
      total_weight += path_work_vector.back().second;
      if (total_weight >= best_weight) {
        break;
      }
    }
    if (total_weight < best_weight &&
        path_work_vector.back().first == vertex2 &&
        path_work_vector.size() == number_of_vertices_in_path) {
      best_weight = total_weight;
      path = std::move(path_work_vector);
    }
  }
  // We now should have a valid path.
  REQUIRE(path.size() == number_of_vertices_in_path);
  REQUIRE(path[0].first == vertex1);
  REQUIRE(path.back().first == vertex2);

  REQUIRE(path[0].second == 0);
  REQUIRE(best_weight == get_total_weight(path));
}

std::pair<unsigned, unsigned> WeightedSquareGrid::get_xy(VertexWSM vv) const {
  const unsigned yy = vv / (m_width + 1);
  const unsigned xx = vv - yy * (m_width + 1);
  REQUIRE(vv == get_vertex(xx, yy));
  return std::make_pair(xx, yy);
}

}  // namespace tests
}  // namespace InitialPlacement
}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
