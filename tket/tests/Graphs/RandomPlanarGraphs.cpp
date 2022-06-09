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

#include "RandomPlanarGraphs.hpp"

#include <catch2/catch_test_macros.hpp>

#include "Utils/RNG.hpp"

using std::vector;

namespace tket {
namespace graphs {
namespace tests {

void RandomPlanarGraphs::reset(size_t width) {
  m_width = 5;
  m_width = std::max(width, m_width);

  m_remaining_gates.clear();
  m_remaining_gates.reserve(2 * m_width * (m_width - 1));

  Gate gate;
  // Think of increasing x as EAST, increasing y as NORTH.
  for (size_t xx = 0; xx < m_width; ++xx) {
    for (size_t yy = 0; yy < m_width; ++yy) {
      gate.vertex1 = yy * m_width + xx;

      if (xx + 1 != m_width) {
        // Join one square to the east.
        gate.vertex2 = gate.vertex1 + 1;
        m_remaining_gates.push_back(gate);
      }
      if (yy + 1 != m_width) {
        // Join one square to the north.
        gate.vertex2 = gate.vertex1 + m_width;
        m_remaining_gates.push_back(gate);
      }
    }
  }
  m_number_of_regions = m_width * m_width;
  m_region_ids.resize(m_number_of_regions);
  for (size_t ii = 0; ii < m_region_ids.size(); ++ii) {
    m_region_ids[ii] = ii;
  }
}

size_t RandomPlanarGraphs::merge_squares(RNG& rng) {
  if (m_remaining_gates.empty()) {
    return m_number_of_regions;
  }
  const Gate gate = rng.get_and_remove_element(m_remaining_gates);

  // The two 1x1 squares are now in the same region.
  const auto id1 = m_region_ids[gate.vertex1];
  const auto id2 = m_region_ids[gate.vertex2];
  if (id1 == id2) {
    return m_number_of_regions;
  }
  // The two region IDs now must be merged together.
  // Not the most efficient here!
  for (size_t& region_id : m_region_ids) {
    if (region_id == id2) {
      region_id = id1;
    }
  }
  --m_number_of_regions;
  return m_number_of_regions;
}

vector<vector<size_t>> RandomPlanarGraphs::get_region_data() const {
  m_region_data.resize(m_number_of_regions);
  for (auto& region_list : m_region_data) {
    region_list.clear();
  }
  m_old_id_to_new_id.clear();
  for (auto region_id : m_region_ids) {
    if (m_old_id_to_new_id.count(region_id) == 0) {
      const auto new_id = m_old_id_to_new_id.size();
      m_old_id_to_new_id[region_id] = new_id;
    }
  }
  REQUIRE(m_old_id_to_new_id.size() == m_number_of_regions);

  for (size_t xx = 0; xx < m_width; ++xx) {
    for (size_t yy = 0; yy < m_width; ++yy) {
      const size_t square = xx + yy * m_width;
      if (xx + 1 < m_width) {
        register_touching_regions(square, square + 1);
      }
      if (yy + 1 < m_width) {
        register_touching_regions(square, square + m_width);
      }
    }
  }
  return m_region_data;
}

void RandomPlanarGraphs::register_touching_regions(
    size_t square1, size_t square2) const {
  const auto id1 = m_old_id_to_new_id.at(m_region_ids.at(square1));
  const auto id2 = m_old_id_to_new_id.at(m_region_ids.at(square2));
  if (id1 != id2) {
    m_region_data[id1].push_back(id2);
  }
}

}  // namespace tests
}  // namespace graphs
}  // namespace tket
