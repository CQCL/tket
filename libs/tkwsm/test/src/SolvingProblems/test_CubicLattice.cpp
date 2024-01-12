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

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <random>
#include <tkwsm/Common/GeneralUtils.hpp>

#include "../TestUtils/CheckedSolution.hpp"
#include "../TestUtils/ResumedSolutionChecker.hpp"
#include "../TestUtils/TestSettings.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {
namespace {
// Represents the 3D points (x,y,z) with |x|,|y|,|z| <= k,
// x,y,z,k integers.
//
// We're going to try all self-embeddings. We know there are always 48,
// so we can easily exhaustively find the best solution.
class CubicLattice {
 public:
  explicit CubicLattice(unsigned k) {
    REQUIRE(k > 0);
    REQUIRE(k < 20);
    m_max_value = k;
    m_num_vertices = 2 * m_max_value + 1;
    m_num_vertices *= m_num_vertices * m_num_vertices;
  }

  void get_xyz(int& x, int& y, int& z, unsigned index) const {
    REQUIRE(index < m_num_vertices);
    // I = z' + Ky' + K^2.x'
    // sI = I/K = y' + Kx'
    const int index_shifted = index / (2 * m_max_value + 1);
    // x' = sI/K = (y' + Kx')/K
    const int x_dashed = index_shifted / (2 * m_max_value + 1);
    REQUIRE(x_dashed <= 2 * m_max_value);
    // y' = sI-kx'
    const int y_dashed = index_shifted - (2 * m_max_value + 1) * x_dashed;
    // z' = I - K(y'+Kx').
    const int z_dashed = index - (2 * m_max_value + 1) * index_shifted;
    x = x_dashed - m_max_value;
    y = y_dashed - m_max_value;
    z = z_dashed - m_max_value;
    REQUIRE(std::abs(x) <= m_max_value);
    REQUIRE(std::abs(y) <= m_max_value);
    REQUIRE(std::abs(z) <= m_max_value);
  }

  std::array<int, 3> get_xyz(unsigned index) const {
    std::array<int, 3> xyz;
    get_xyz(xyz[0], xyz[1], xyz[2], index);
    return xyz;
  }

  int get_max_value() const { return m_max_value; }

  unsigned get_index(int x, int y, int z) const {
    REQUIRE(std::abs(x) <= m_max_value);
    REQUIRE(std::abs(y) <= m_max_value);
    REQUIRE(std::abs(z) <= m_max_value);
    unsigned index = x + m_max_value;
    index *= (2 * m_max_value + 1);
    index += y + m_max_value;
    index *= (2 * m_max_value + 1);
    index += z + m_max_value;
    return index;
  }

  unsigned get_index(const std::array<int, 3>& point) const {
    return get_index(point[0], point[1], point[2]);
  }

  std::set<EdgeWSM> get_edges() const {
    std::set<EdgeWSM> edges;

    for (int x = -m_max_value; x <= m_max_value; ++x) {
      for (int y = -m_max_value; y <= m_max_value; ++y) {
        for (int z = -m_max_value; z <= m_max_value; ++z) {
          add_edge(x, y, z, x + 1, y, z, edges);
          add_edge(x, y, z, x, y + 1, z, edges);
          add_edge(x, y, z, x, y, z + 1, edges);
        }
      }
    }
    // How many edges? Use symmetry.
    const unsigned horiz_edges_in_each_plane =
        2 * m_max_value * (2 * m_max_value + 1);
    const unsigned horiz_planes = 2 * m_max_value + 1;
    const unsigned total_horiz_edges = horiz_edges_in_each_plane * horiz_planes;
    CHECK(edges.size() == 3 * total_horiz_edges);
    return edges;
  }

  unsigned get_n_vertices() const { return m_num_vertices; }

 private:
  int m_max_value;
  unsigned m_num_vertices;

  // Horribly inelegant!
  void add_edge(
      int x, int y, int z, int x2, int y2, int z2,
      std::set<EdgeWSM>& edges) const {
    if (!(std::abs(x2) <= m_max_value && std::abs(y2) <= m_max_value &&
          std::abs(z2) <= m_max_value)) {
      return;
    }
    edges.insert(get_edge(get_index(x, y, z), get_index(x2, y2, z2)));
  }
};
}  // namespace

static void add_random_weights(
    GraphEdgeWeights& edges_and_weights, std::mt19937_64& rng_64) {
  // Let's take weights 1,4,9,...,64, stored in 3 bits...
  std::uint64_t bits = 0;
  for (auto& entry : edges_and_weights) {
    // who cares about a bit of bias...
    if (bits == 0) {
      bits = rng_64();
    }
    const auto rnd_number = bits & 0x7;
    bits >>= 3;
    REQUIRE(rnd_number >= 0);
    REQUIRE(rnd_number <= 7);
    entry.second = (rnd_number + 1) * (rnd_number + 1);
  }
}

/*
See, of course,

https://en.wikipedia.org/wiki/Octahedral_symmetry

There are 48 possible orthogonal matrices mapping the cube {|x(i)| <= K}
to itself [thus, keeping (0,0,0) fixed].

They are simply the 3! = 6 permutations of (x,y,z),
combined with the 2^3 = 8  triples of +/- signs.

*/

// Assuming that the initial edge weights for the lattice are filled,
// perform all possible rotations/reflections on the cube
// to transform the points, and append them to the vector.
static void append_all_transformed_cubes_data(
    const CubicLattice& lattice,
    std::vector<GraphEdgeWeights>& list_of_weights_data) {
  REQUIRE(list_of_weights_data.size() == 1);

  // Start with the identity permutation.
  std::array<unsigned, 3> permutations;
  for (unsigned ii = 0; ii < 3; ++ii) {
    permutations[ii] = ii;
  }

  // The +/- signs. 1 means -, 0 means +.
  unsigned signs_code;

  const auto transform_point = [&permutations,
                                &signs_code](std::array<int, 3>& point) {
    const auto orig_copy = point;
    for (unsigned ii = 0; ii < 3; ++ii) {
      point[ii] = orig_copy[permutations[ii]];
      if ((signs_code & (1u << ii)) != 0) {
        point[ii] = -point[ii];
      }
    }
  };

  do {
    for (signs_code = 0; signs_code <= 7; ++signs_code) {
      if (list_of_weights_data.size() == 1 && signs_code == 0) {
        continue;
      }
      list_of_weights_data.emplace_back();
      const auto& original_data = list_of_weights_data[0];
      auto& new_data = list_of_weights_data.back();
      for (const auto& original_entry : original_data) {
        const auto& original_edge = original_entry.first;
        const auto& p1_index = original_edge.first;
        const auto& p2_index = original_edge.second;
        auto p1 = lattice.get_xyz(p1_index);
        auto p2 = lattice.get_xyz(p2_index);
        transform_point(p1);
        transform_point(p2);
        const auto new_index1 = lattice.get_index(p1);
        const auto new_index2 = lattice.get_index(p2);
        new_data[get_edge(new_index1, new_index2)] = original_entry.second;
      }
    }
  } while (std::next_permutation(permutations.begin(), permutations.end()));

  REQUIRE(list_of_weights_data.size() == 48);
  for (const auto& data : list_of_weights_data) {
    // They should all have the same keys.
    REQUIRE(data.size() == list_of_weights_data[0].size());
    for (const auto& entry : data) {
      REQUIRE(list_of_weights_data[0].count(entry.first) != 0);
    }
  }
}

// Already checked to be valid.
static WeightWSM get_optimal_self_embedding(
    const std::vector<GraphEdgeWeights>& list_of_weights_data) {
  WeightWSM best_weight;
  set_maximum(best_weight);
  for (const auto& transformed_data : list_of_weights_data) {
    WeightWSM total_weight = 0;
    for (const auto& entry : transformed_data) {
      total_weight += entry.second * list_of_weights_data[0].at(entry.first);
    }
    best_weight = std::min(best_weight, total_weight);
  }
  REQUIRE(!is_maximum(best_weight));
  return best_weight;
}

SCENARIO("Check cubic lattice indexing") {
  const CubicLattice lattice(1);
  std::array<int, 3> p;
  for (p[0] = -1; p[0] <= 1; p[0] += 1) {
    for (p[1] = -1; p[1] <= 1; p[1] += 1) {
      for (p[2] = -1; p[2] <= 1; p[2] += 1) {
        const auto index = lattice.get_index(p);
        const auto new_p = lattice.get_xyz(index);
        // INFO("index=" << index << ": (" << p[0] << ","
        //      << p[1] << "," << p[2] << ")");
        CHECK(p == new_p);
        CHECK(index < lattice.get_n_vertices());
      }
    }
  }
}

/* Performance NOTE:

A typical printout with older versions was:

cubic lattice: k=1, V=27, E=54, opt. soln 31870; time 1+8
cubic lattice: k=2, V=125, E=300, opt. soln 162234; time 20+111
cubic lattice: k=3, V=343, E=882, opt. soln 541713; time 142+473
cubic lattice: k=4, V=729, E=1944, opt. soln 1258558; time 546+1399
cubic lattice: k=5, V=1331, E=3630, opt. soln 2266498; time 1619+3378
@@@ fin. time 2328+5369

HOWEVER, it is now:

cubic lattice: k=1, V=27, E=54, opt. soln 31870; time 0+9; 25 iters; known
opt.val. 31870 cubic lattice: k=2, V=125, E=300, opt. soln 162234; time 3+126;
27 iters; known opt.val. 162234 cubic lattice: k=3, V=343, E=882, opt. soln
541713; time 22+1431; 41 iters; known opt.val. 541713 cubic lattice: k=4, V=729,
E=1944, opt. soln 1258558; time 99+8182; 167 iters; known opt.val. 1258558 cubic
lattice: k=5, V=1331, E=3630, opt. soln 2266498; time 364+10523; 17 iters; known
opt.val. 2266498; TIMED OUT
@@@ Cubic lattice fin. Time 488+20271

So, the newer version is SLOWER than the older version.
However, it is really ONLY the cubic lattices tests which are slower;
almost all other tests [including square grids] are faster.
It's unclear why cubic lattices should suffer like this
(e.g., square grids are also quite regular, homogeneous graphs);
it needs further investigation.
*/

void test_cubic_lattices(
    const std::vector<unsigned>& k_values, bool do_resumption_check) {
  const auto& os = TestSettings::get().os;

  std::mt19937_64 rng_64;
  std::vector<GraphEdgeWeights> list_of_weights_data;
  list_of_weights_data.reserve(48);

  CheckedSolution::Statistics stats("Cubic lattices");
  CheckedSolution::ProblemInformation info;
  const MainSolverParameters solver_params(10000);
  ResumedSolutionChecker resumption_checker;

  for (unsigned k_value : k_values) {
    const CubicLattice lattice(k_value);
    list_of_weights_data.clear();
    list_of_weights_data.resize(1);
    for (const auto& edge : lattice.get_edges()) {
      list_of_weights_data[0][edge];
    }
    add_random_weights(list_of_weights_data[0], rng_64);
    append_all_transformed_cubes_data(lattice, list_of_weights_data);
    const auto optimal_weight =
        get_optimal_self_embedding(list_of_weights_data);
    os << "\nk=" << lattice.get_max_value()
       << ", V=" << lattice.get_n_vertices()
       << ", E=" << list_of_weights_data[0].size() << ", opt. soln "
       << optimal_weight;

    const auto old_init_time = stats.total_init_time_ms;
    const auto old_search_time = stats.total_search_time_ms;

    info.known_optimal_solution = optimal_weight;

    const CheckedSolution checked_solution(
        list_of_weights_data[0], list_of_weights_data[0], info, solver_params,
        stats);

    if (do_resumption_check) {
      resumption_checker.check(
          checked_solution, list_of_weights_data[0], list_of_weights_data[0],
          solver_params);
    }
  }
  stats.finish();
  CHECK(stats.success_count == k_values.size());
  CHECK(stats.failure_count == 0);
  CHECK(stats.timeout_count == 0);
}

SCENARIO("Self-embed cubic lattices - quicker test") {
  const std::vector<unsigned> k_values{1, 2};
  test_cubic_lattices(k_values, true);
}

SCENARIO("Self-embed cubic lattices - slower test") {
  // k=3, V=343, E=882 is currently around 1.5 secs; slower than we'd like
  const std::vector<unsigned> k_values{3};
  test_cubic_lattices(k_values, false);
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
