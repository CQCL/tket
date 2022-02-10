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

#include <catch2/catch.hpp>
#include <random>

#include "../testutil.hpp"
#include "Architecture/Architecture.hpp"
#include "PlacementWithWSM/CalculatedPlacementMap.hpp"

// These tests are just copied from the old Placement tests.
// TODO make extensive tests!

namespace tket {
namespace test_Placement {

typedef std::map<Qubit, Node> QubitMapping;

static std::string get_qubit_map_str(const QubitMapping& map) {
  std::stringstream ss;
  for (const auto& entry : map) {
    ss << "(" << entry.first.repr() << " -> " << entry.second.repr() << ") ";
  }
  return ss.str();
}

static void test_wsm_get_placement_map(
    const Circuit& test_circ, const Architecture& test_arc,
    const std::string& map_encoding) {
  const CalculatedPlacementMap calc_map(test_circ, test_arc);
  const auto& test_m = calc_map.placement_map;
  const auto test_m_str = get_qubit_map_str(test_m);

  CHECK(test_m_str == map_encoding);

  // A complete placement should have been made.
  const qubit_vector_t all_qs = test_circ.all_qubits();
  CHECK(test_m.size() == all_qs.size());
  for (const auto& logical_qubit : all_qs) {
    CHECK(test_m.count(logical_qubit) != 0);
  }
  // For such small graphs, it should be very fast.
  CHECK(calc_map.full_placement_result.total_init_time_ms <= 10);
  CHECK(calc_map.full_placement_result.total_search_time_ms <= 10);
}

SCENARIO("Old get_placement tests, with WSM instead") {
  GIVEN("Old LinePlacement test data") {
    // Line 0-1-2-3
    const Architecture test_arc({{0, 1}, {1, 2}, {2, 3}});

    Circuit test_circ(4);
    // P-graph is a Y-shape, central vertex 1.
    add_2qb_gates(test_circ, OpType::CX, {{0, 1}, {2, 1}, {3, 1}});
    test_wsm_get_placement_map(
        test_circ, test_arc,
        // Manually checked: this is jointly, but not uniquely, best possible
        "(q[0] -> node[3]) (q[1] -> node[2]) (q[2] -> node[1]) (q[3] -> "
        "node[0]) ");
  }
  GIVEN("Old GraphPlacement and NoiseAwarePlacement test data") {
    // No obvious shape - just draw it and see!
    const Architecture test_arc(
        {{0, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 5}});

    Circuit test_circ(6);
    add_2qb_gates(
        test_circ, OpType::CX,
        {{0, 1}, {2, 1}, {3, 1}, {2, 5}, {3, 4}, {0, 5}});
    test_wsm_get_placement_map(
        test_circ, test_arc,
        // Plausible enough: all circuit P-edges except (q3,q4), (q0,q5)
        // are respected; and (q3,q4), (q0,q5) are only distance 2 apart in the
        // target architecture.
        "(q[0] -> node[3]) (q[1] -> node[1]) (q[2] -> node[2]) (q[3] -> "
        "node[4]) (q[4] -> node[0]) (q[5] -> node[5]) ");
  }
}

}  // namespace test_Placement
}  // namespace tket
