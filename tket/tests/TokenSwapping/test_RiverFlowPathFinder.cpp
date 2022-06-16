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

#include <array>
#include <catch2/catch_test_macros.hpp>

#include "Architecture/ArchitectureMapping.hpp"
#include "Architecture/DistancesFromArchitecture.hpp"
#include "Architecture/NeighboursFromArchitecture.hpp"
#include "TestUtils/ArchitectureEdgesReimplementation.hpp"
#include "TokenSwapping/RiverFlowPathFinder.hpp"
#include "Utils/RNG.hpp"

using std::vector;

namespace tket {
namespace tsa_internal {
namespace tests {

namespace {
// It is a cycle (ring) on vertices [0,1,2,..., N-1], with N ~ 0.
struct DistancesForCycle : public DistancesInterface {
  size_t number_of_vertices = 10;

  virtual size_t operator()(size_t v1, size_t v2) override {
    size_t distance1;
    if (v1 < v2) {
      distance1 = v2 - v1;
    } else {
      distance1 = v1 - v2;
    }
    const size_t distance2 = number_of_vertices - distance1;
    return std::min(distance1, distance2);
  }
};

class NeighboursForCycle : public NeighboursInterface {
 public:
  explicit NeighboursForCycle(size_t number_of_vertices)
      : m_number_of_vertices(number_of_vertices) {
    REQUIRE(number_of_vertices > 1);
    if (m_number_of_vertices == 2) {
      m_neighbours.resize(1);
    } else {
      m_neighbours.resize(2);
    }
  }

  virtual const vector<size_t>& operator()(size_t vertex) override {
    if (vertex >= m_number_of_vertices) {
      throw std::runtime_error("neighbours requested for invalid vertex");
    }
    m_neighbours[0] = (vertex + 1) % m_number_of_vertices;
    if (m_neighbours.size() > 1) {
      m_neighbours[1] =
          ((vertex + m_number_of_vertices) - 1) % m_number_of_vertices;
    }
    return m_neighbours;
  }

 private:
  size_t m_number_of_vertices;
  vector<size_t> m_neighbours;
};

struct TestResult {
  size_t total_number_of_path_calls = 0;
  size_t total_number_of_differing_extra_paths = 0;

  std::string str() const {
    std::stringstream ss;
    ss << "[ Number of path calls: " << total_number_of_path_calls
       << "  Extra paths: " << total_number_of_differing_extra_paths << " ]";
    return ss.str();
  }
};

}  // namespace

static void do_simple_path_test(
    const vector<size_t>& path, const Swap& endpoints) {
  REQUIRE(!path.empty());
  REQUIRE(path[0] == endpoints.first);
  REQUIRE(path.back() == endpoints.second);

  const std::set<size_t> vertices{path.cbegin(), path.cend()};
  REQUIRE(vertices.size() == path.size());
}

static void require_path_to_have_valid_edges(
    const vector<size_t>& path, NeighboursInterface& neighbours_interface) {
  std::array<Swap, 2> vertices;
  for (size_t ii = 0; ii + 1 < path.size(); ++ii) {
    vertices[0].first = path[ii];
    vertices[0].second = path[ii + 1];
    vertices[1].first = path[ii + 1];
    vertices[1].second = path[ii];
    for (const auto& pair : vertices) {
      const auto& neighbours = neighbours_interface(pair.first);
      bool is_neighbour = false;
      for (auto neigh : neighbours) {
        if (neigh == pair.second) {
          is_neighbour = true;
          break;
        }
      }
      REQUIRE(is_neighbour);
    }
  }
}

static void test(
    TestResult& result, RiverFlowPathFinder& path_finder,
    DistancesInterface& distance_calculator,
    NeighboursInterface& neighbours_calculator, size_t number_of_vertices,
    RNG& rng_for_test_data, size_t number_of_test_repeats = 10) {
  // We will check that calculated paths are mostly unchanged.
  std::map<Swap, vector<vector<size_t>>> calculated_paths;

  vector<Swap> possible_path_calls;
  possible_path_calls.reserve(number_of_vertices * number_of_vertices);
  for (size_t ii = 0; ii < number_of_vertices; ++ii) {
    for (size_t jj = 0; jj < number_of_vertices; ++jj) {
      possible_path_calls.emplace_back(ii, jj);
      calculated_paths[std::make_pair(ii, jj)];
    }
  }

  // The first time a path is calculated, its length will be checked using
  // the distance_calculator
  const auto get_path_size = [&calculated_paths, &distance_calculator](
                                 const Swap& end_vertices) -> size_t {
    if (end_vertices.first == end_vertices.second) {
      return 1;
    }
    const auto& existing_paths = calculated_paths[end_vertices];
    if (!existing_paths.empty()) {
      return existing_paths[0].size();
    }
    const auto& reversed_existing_paths = calculated_paths[std::make_pair(
        end_vertices.second, end_vertices.first)];

    if (!reversed_existing_paths.empty()) {
      return reversed_existing_paths[0].size();
    }
    return 1 + distance_calculator(end_vertices.first, end_vertices.second);
  };

  for (size_t counter = number_of_test_repeats; counter > 0; --counter) {
    rng_for_test_data.do_shuffle(possible_path_calls);
    result.total_number_of_path_calls += possible_path_calls.size();

    for (const Swap& end_vertices : possible_path_calls) {
      const auto& calc_path =
          path_finder(end_vertices.first, end_vertices.second);

      do_simple_path_test(calc_path, end_vertices);
      REQUIRE(calc_path.size() == get_path_size(end_vertices));

      auto& path_list = calculated_paths[end_vertices];
      bool found_path = false;
      for (auto& path : path_list) {
        if (path == calc_path) {
          found_path = true;
          break;
        }
      }
      if (!found_path) {
        if (!path_list.empty()) {
          ++result.total_number_of_differing_extra_paths;
        }
        path_list.emplace_back(calc_path);
        require_path_to_have_valid_edges(calc_path, neighbours_calculator);
      }
    }
  }
}

SCENARIO("Test path generation for cycles") {
  RNG rng_for_path_generation;
  RNG rng_for_test_data;
  DistancesForCycle distances;
  TestResult result;

  for (size_t number_of_vertices = 2; number_of_vertices <= 10;
       ++number_of_vertices) {
    INFO("number_of_vertices = " << number_of_vertices);
    distances.number_of_vertices = number_of_vertices;
    NeighboursForCycle neighbours(number_of_vertices);
    RiverFlowPathFinder path_finder(
        distances, neighbours, rng_for_path_generation);

    const auto current_differing_paths =
        result.total_number_of_differing_extra_paths;
    test(
        result, path_finder, distances, neighbours, number_of_vertices,
        rng_for_test_data);

    // Even cycles have non-unique paths, for polar opposite vertices;
    // odd cycles do not.
    if (number_of_vertices % 2 == 1) {
      // No extra paths were created.
      CHECK(
          current_differing_paths ==
          result.total_number_of_differing_extra_paths);
    }
  }
  REQUIRE(result.str() == "[ Number of path calls: 3840  Extra paths: 3 ]");
}

// Deliberately use the same RNG, so it's all mixed up;
// but we still expect not so many different paths.
static void test(
    TestResult& result, const ArchitectureMapping& arch_mapping, RNG& rng) {
  DistancesFromArchitecture distances(arch_mapping);
  NeighboursFromArchitecture neighbours(arch_mapping);
  RiverFlowPathFinder path_finder(distances, neighbours, rng);

  test(
      result, path_finder, distances, neighbours,
      arch_mapping.number_of_vertices(), rng);
}

SCENARIO("Path generation for ring graph") {
  RNG rng;
  TestResult result;
  const RingArch arch(7);
  const ArchitectureMapping arch_mapping(arch);
  test(result, arch_mapping, rng);
  REQUIRE(result.str() == "[ Number of path calls: 490  Extra paths: 0 ]");
}

SCENARIO("Path generation for square grids") {
  RNG rng;
  TestResult result;
  for (size_t ver = 2; ver <= 4; ver += 2) {
    for (size_t hor = 1; hor <= 5; hor += 2) {
      for (size_t layer = 1; layer <= 3; layer += 2) {
        const auto edges = get_square_grid_edges(ver, hor, layer);
        const Architecture arch(edges);
        const ArchitectureMapping arch_mapping(arch, edges);
        test(result, arch_mapping, rng);
      }
    }
  }
  REQUIRE(result.str() == "[ Number of path calls: 70000  Extra paths: 583 ]");
}

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
