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

#include <catch2/catch_test_macros.hpp>
#include <istream>
#include <map>
#include <random>
#include <string>
#include <tkwsm/Common/GeneralUtils.hpp>
#include <utility>

#include "../TestUtils/CheckedSolution.hpp"
#include "../TestUtils/TestSettings.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

namespace {
struct EdgeWeightAndRand {
  EdgeWSM edge;
  std::uint64_t rand_num;
};

struct TestResult {
  unsigned success_count = 0;
  unsigned failure_count = 0;
  unsigned timeout_count = 0;
  long long total_time_ms = 0;
  unsigned total_edges = 0;
  unsigned total_verts = 0;
};

struct TestParameters {
  unsigned timeout = 1000;

  // If nonzero, will check that the total time taken
  // is within these limits.
  unsigned expected_max_total_time_ms = 0;
  unsigned expected_min_total_time_ms = 0;

  // Extra checks against changing tests. Fill with the values, if known.
  unsigned total_number_of_vertices = 0;
  unsigned total_number_of_edges = 0;
};
}  // namespace

// The code is simply 3 numbers separated by spaces, like "30 1000 1111".
// They are: number of vertices; number of edges; the rng seed.
// Should be completely platform/compiler independent, as the behaviour
// of the Mersenne twister engine mt19937_64
// is guaranteed by the standard,
// even though distributions (conversion of raw 64 bits into,
// e.g., an approximately uniform int) are NOT.
// The weights will be assigned from less significant bits of the
// random numbers, so there is basically zero correlation with the sorting
// order (determined almost surely by much more significant bits).
static GraphEdgeWeights get_graph_data(
    const std::string& code, const std::vector<WeightWSM>& weights) {
  unsigned number_of_vertices = 0;
  unsigned number_of_edges = 0;
  unsigned seed = 0;
  std::istringstream iss(code);
  iss >> number_of_vertices;
  iss >> number_of_edges;
  iss >> seed;
  REQUIRE(number_of_vertices >= 5);
  REQUIRE(number_of_vertices <= 1000);
  REQUIRE(number_of_edges >= number_of_vertices);
  const unsigned expected_size =
      (number_of_vertices * (number_of_vertices - 1)) / 2;
  REQUIRE(number_of_edges <= expected_size);
  REQUIRE(seed <= 1000000);
  unsigned weights_mask;
  switch (weights.size()) {
    case 2:
      weights_mask = 1;
      break;
    case 4:
      weights_mask = 3;
      break;
    case 8:
      weights_mask = 7;
      break;
    default:
      REQUIRE(false);
      break;
  }
  for (auto ww : weights) {
    REQUIRE(ww > 0);
    REQUIRE(ww <= 1000);
  }
  std::mt19937_64 rng;
  rng.seed(seed);

  std::vector<EdgeWeightAndRand> data_vector;
  data_vector.reserve(expected_size);
  for (unsigned ii = 0; ii < number_of_vertices; ++ii) {
    for (unsigned jj = ii + 1; jj < number_of_vertices; ++jj) {
      data_vector.emplace_back();
      data_vector.back().edge = get_edge(ii, jj);
      // The raw 64 bits should be fully guaranteed by the standard...
      data_vector.back().rand_num = rng();
    }
  }
  REQUIRE(data_vector.size() == expected_size);
  std::sort(
      data_vector.begin(), data_vector.end(),
      [](const EdgeWeightAndRand& lhs, const EdgeWeightAndRand& rhs) -> bool {
        // Fully portable even with nonstable sort, as no duplicate elems
        return (lhs.rand_num < rhs.rand_num) ||
               (lhs.rand_num == rhs.rand_num && lhs.edge < rhs.edge);
      });

  GraphEdgeWeights result;
  for (unsigned ii = 0; ii < number_of_edges; ++ii) {
    result[data_vector[ii].edge] =
        weights.at((data_vector[ii].rand_num >> 2) & weights_mask);
  }
  return result;
}

// Copied from boost hash combine.
// Very annoyingly, boost hash combine is useless for us
// because it takes size_t, which varies across platforms
static void hash_combine(std::uint32_t& seed, std::uint32_t v) {
  seed ^= v + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

// Used to check that the final generated data hasn't changed.
// However, shortened by a few bits to ensure that it fits within
// a SIGNED int.
static std::uint32_t get_weights_hash(const GraphEdgeWeights& data) {
  std::uint32_t max_value;
  set_maximum(max_value);
  std::uint32_t result = 1;
  const auto add_value = [&result, max_value](std::intmax_t val) {
    REQUIRE(val <= max_value);
    hash_combine(result, val);
  };
  for (const auto& entry : data) {
    add_value(entry.first.first);
    add_value(entry.first.second);
    add_value(entry.second);
  }
  return (result >> 2) ^ (result & 3);
}

typedef std::vector<std::intmax_t> ResultsSummary;

namespace {
struct AllAgainstAllTester {
  std::vector<std::string> codes;
  std::vector<WeightWSM> weights;

  // expected_results should FIRST list the hashes of the graphs with weights,
  // THEN the scalar products S in order: S=0 means no solution,
  // S>0 means an optimal solution with S was found.
  // -1 means a timeout; we have never actually completed the solution,
  // so don't actually know the solution.
  std::vector<std::intmax_t> expected_results;

  TestResult test_all_against_all(const TestParameters& params) const {
    TestResult result;
    ResultsSummary calc_results;
    std::vector<GraphEdgeWeights> graphs;
    std::vector<unsigned> num_vertices;
    graphs.reserve(codes.size());
    for (const auto& code : codes) {
      graphs.emplace_back(get_graph_data(code, weights));
      result.total_edges += graphs.back().size();
      num_vertices.push_back(get_vertices(graphs.back()).size());
      result.total_verts += num_vertices.back();
      calc_results.push_back(get_weights_hash(graphs.back()));
    }
    const auto& os = TestSettings::get().os;
    if (params.total_number_of_edges != 0) {
      CHECK(params.total_number_of_edges == result.total_edges);
    }
    if (params.total_number_of_vertices != 0) {
      CHECK(params.total_number_of_vertices == result.total_verts);
    }
    std::stringstream ss;
    ss << "random graphs: " << graphs.size() << " graphs; "
       << result.total_edges << " edges; " << result.total_verts
       << " vertices; timeout " << params.timeout;

    CheckedSolution::Statistics statistics(ss.str());
    const MainSolverParameters solver_params(params.timeout);

    const auto update_calc_results =
        [&calc_results](const CheckedSolution& checked_solution) {
          if (!checked_solution.finished) {
            // A timeout.
            calc_results.push_back(-1);
            return;
          }
          // If no solution, it will be zero.
          calc_results.push_back(checked_solution.scalar_product);
        };

    for (unsigned ii = 0; ii < graphs.size(); ++ii) {
      for (unsigned jj = 0; jj < graphs.size(); ++jj) {
        if (calc_results.size() < expected_results.size() &&
            expected_results[calc_results.size()] == -1) {
          // It's known to be a timeout, so don't bother again.
          calc_results.push_back(-1);
          ++statistics.timeout_count;
          continue;
        }

        os << "\n#### embedding: G[" << ii << "]: (V=" << num_vertices[ii]
           << ",E=" << graphs[ii].size() << ") -> G[" << jj
           << "]: (V=" << num_vertices[jj] << ",E=" << graphs[jj].size() << ")";

        CheckedSolution::ProblemInformation info;
        if (ii == jj) {
          // Self embedding is always possible,
          // although we do not know the OPTIMAL solution.
          WeightWSM total_w = 0;
          for (const auto& entry : graphs[ii]) {
            total_w += entry.second * entry.second;
          }
          info.known_upper_bound = total_w;
          info.existence = CheckedSolution::ProblemInformation::
              SolutionsExistence::KNOWN_TO_BE_SOLUBLE;
        }

        const CheckedSolution checked_solution(
            graphs[ii], graphs[jj], info, solver_params, statistics);
        update_calc_results(checked_solution);
      }
    }
    result.total_time_ms =
        statistics.total_init_time_ms + statistics.total_search_time_ms;
    result.failure_count = statistics.failure_count;
    result.timeout_count = statistics.timeout_count;
    result.success_count = statistics.success_count;
    statistics.finish(
        CheckedSolution::Statistics::Expectation::ALL_SUCCESS_OR_TIMEOUT);

    CHECK(expected_results == calc_results);
    if (params.expected_max_total_time_ms != 0) {
      CHECK(result.total_time_ms <= params.expected_max_total_time_ms);
    }
    CHECK(result.total_time_ms >= params.expected_min_total_time_ms);
    return result;
  }
};
}  // namespace

SCENARIO("embedding random graphs - smaller graphs, small weights") {
  AllAgainstAllTester tester;
  tester.codes = {"5 8 111",  "5 9 12211", "6 10 13311", "7 10 222",
                  "7 15 333", "8 16 1111", "8 20 444",   "10 20 333"};
  tester.weights = {1, 2, 3, 8};
  tester.expected_results = {
      820581231, 797760108, 317578032, 996088179, 905537177, 505148537,
      63334049,  630164384, 87,        89,        0,         0,
      49,        67,        35,        45,        0,         222,
      0,         0,         116,       182,       99,        124,
      0,         0,         58,        0,         98,        0,
      49,        76,        0,         0,         0,         161,
      71,        99,        54,        63,        0,         0,
      0,         0,         279,       0,         155,       163,
      0,         0,         0,         0,         0,         425,
      0,         0,         0,         0,         0,         0,
      0,         0,         174,       0,         0,         0,
      0,         0,         0,         0,         0,         279};

  TestParameters params;
  params.timeout = 1000;
  params.total_number_of_edges = 108;
  params.total_number_of_vertices = 56;

  // Currently, ~10 ms.
  // However, Valgrind is much slower than normal runs.
  params.expected_max_total_time_ms = 1000 * 10;

  const auto result = tester.test_all_against_all(params);
  CHECK(result.success_count == 64);
  CHECK(result.failure_count == 0);
  CHECK(result.timeout_count == 0);
}

SCENARIO("embedding random graphs - medium graphs, small weights") {
  AllAgainstAllTester tester;
  tester.codes = {"10 20 1111", "10 30 2222",  "11 20 3333", "11 40 4444",
                  "15 30 5555", "16 50, 6666", "17 60 7777", "18 70 888"};
  tester.weights = {1, 2, 3, 8};
  tester.expected_results = {
      724217328, 705349590, 154711899, 916605139, 166486361, 875669872,
      325817875, 806972053, 411,       182,       0,         122,
      0,         259,       146,       128,       0,         616,
      0,         310,       0,         0,         0,         0,
      0,         0,         228,       100,       0,         192,
      97,        98,        0,         0,         0,         575,
      0,         0,         0,         0,         0,         0,
      0,         0,         590,       0,         278,       194,
      0,         0,         0,         0,         0,         1338,
      0,         0,         0,         0,         0,         0,
      0,         0,         1068,      0,         0,         0,
      0,         0,         0,         0,         0,         1257};
  TestParameters params;
  params.timeout = 10000;
  params.total_number_of_edges = 320;
  params.total_number_of_vertices = 108;

  // Currently, ~600 ms.
  params.expected_max_total_time_ms = 5000;
  params.expected_min_total_time_ms = 10;

  const auto result = tester.test_all_against_all(params);
  CHECK(result.success_count == 64);
  CHECK(result.failure_count == 0);
  CHECK(result.timeout_count == 0);
}

static AllAgainstAllTester get_large_graphs_small_weights_data() {
  AllAgainstAllTester tester;
  tester.codes = {"20 50 1111",  "22 80 2222",    "25 120 3333", "25 200 4444",
                  "30 200 5555", "32 300 6666",   "35 300 7777", "40 500 8888",
                  "50 500 9999", "55 1000 101010"};

  tester.weights = {1, 2, 3, 8};

  // Any value < -1 (e.g., -9999) means that we DID, once,
  // compute that there's NO solution, but it took a long time.
  tester.expected_results = {
      460517071,  255540664,  811304662, 581415081, 853453591, 367941120,
      1072813581, 1006422874, 309411091, 971368384, 1261,      0,
      619,        -1,         -1,        -1,        -1,        -1,
      -1,         -1,         0,         1732,      0,         -1,
      -9999,      -1,         -1,        -1,        -1,        -1,
      0,          0,          2132,      -1,        0,         -1,
      -1,         -1,         -1,        -1,        0,         0,
      0,          3463,       0,         0,         0,         -1,
      0,          -1,         0,         0,         0,         0,
      3955,       -1,         0,         -1,        0,         -1,
      0,          0,          0,         0,         0,         5758,
      0,          -1,         0,         -1,        0,         0,
      0,          0,          0,         0,         5869,      -1,
      0,          -1,         0,         0,         0,         0,
      0,          0,          0,         10612,     0,         -1,
      0,          0,          0,         0,         0,         0,
      0,          0,          9721,      -1,        0,         0,
      0,          0,          0,         0,         0,         0,
      0,          18810};
  return tester;
}

SCENARIO(
    "embedding random graphs - 2 nasty problems: large graphs, small weights") {
  auto tester = get_large_graphs_small_weights_data();
  for (auto& entry : tester.expected_results) {
    if (entry < -1) {
      // It's the no-solution problem, ~2 seconds.
      entry = 0;
      continue;
    }
    if (entry == 619 || entry > 1000000) {
      // Either the problem taking ~5 seconds, or a code number.
      continue;
    }
    // We'll just skip all other problems; pretend they're timeouts.
    entry = -1;
  }
  TestParameters params;
  params.timeout = 60000;
  params.total_number_of_edges = 3250;
  params.total_number_of_vertices = 334;

  params.expected_max_total_time_ms = 60000;
  params.expected_min_total_time_ms = 1;
  const auto result = tester.test_all_against_all(params);
  CHECK(result.success_count == 2);
  CHECK(result.failure_count == 0);
  CHECK(result.timeout_count == 98);
}

SCENARIO(
    "embedding random graphs - large graphs, small weights, shorter problems") {
  auto tester = get_large_graphs_small_weights_data();
  for (auto& entry : tester.expected_results) {
    if (entry < -1 || entry == 619) {
      // It's one of the longer problems; skip.
      entry = -1;
    }
  }
  TestParameters params;
  params.timeout = 1000;
  params.total_number_of_edges = 3250;
  params.total_number_of_vertices = 334;

  // test coverage takes longer than normal running.
  params.expected_max_total_time_ms = 20 * 1000;
  params.expected_min_total_time_ms = 10;
  const auto result = tester.test_all_against_all(params);
  CHECK(result.success_count == 68);
  CHECK(result.failure_count == 0);
  CHECK(result.timeout_count == 32);
}

static AllAgainstAllTester get_mixed_sizes_problems() {
  AllAgainstAllTester tester;
  tester.codes = {"5 7 111",     "6 14 222",    "10 20 1111",  "10 40 3333",
                  "20 50 4444",  "20 100 5555", "20 150 6666", "30 100 7777",
                  "30 200 8888", "30 400 9999"};

  tester.weights = {1, 2, 5, 20};

  // It so happens that all no-solution problems in this set are fairly quick.
  tester.expected_results = {
      911552196, 461091619, 772140787, 11588550,   1037162436, 766190752,
      951748272, 961275497, 870669976, 1033828678, 117,        91,
      126,       52,        0,         42,         36,         58,
      36,        -1,        0,         753,        0,          292,
      0,         470,       149,       0,          174,        -1,
      0,         0,         2219,      294,        0,          282,
      -1,        477,       181,       -1,         0,          0,
      0,         1304,      0,         0,          -1,         0,
      0,         -1,        0,         0,          0,          0,
      5321,      2461,      -1,        0,          -1,         -1,
      0,         0,         0,         0,          0,          12607,
      -1,        0,         0,         -1,         0,          0,
      0,         0,         0,         0,          15471,      0,
      0,         -1,        0,         0,          0,          0,
      0,         0,         0,         10600,      -1,         -1,
      0,         0,         0,         0,          0,          0,
      0,         0,         21893,     -1,         0,          0,
      0,         0,         0,         0,          0,          0,
      0,         31845};
  return tester;
}

static std::set<unsigned> get_medium_time_problem_values() {
  // Problems taking roughly ~100 ms to ~500ms.
  return {149, 294, 282, 181};
};

static unsigned get_long_problem_value() {
  // Just this one nasty problem takes almost all the time.
  return 2461;
}

SCENARIO(
    "embedding random graphs - mixed sizes and densities: short problems") {
  const auto long_problem_value = get_long_problem_value();
  const auto medium_problem_values = get_medium_time_problem_values();
  auto tester = get_mixed_sizes_problems();
  for (auto& entry : tester.expected_results) {
    if (entry == long_problem_value ||
        medium_problem_values.count(entry) != 0) {
      // It's one of the longer problems; skip.
      entry = -1;
    }
  }
  TestParameters params;
  params.timeout = 1000;
  params.total_number_of_edges = 1081;
  params.total_number_of_vertices = 181;

  // test coverage takes longer than normal running.
  params.expected_max_total_time_ms = 20 * 1000;

  params.expected_min_total_time_ms = 10;
  const auto result = tester.test_all_against_all(params);
  CHECK(result.success_count == 80);
  CHECK(result.failure_count == 0);
  CHECK(result.timeout_count == 20);
}

SCENARIO(
    "embedding random graphs - mixed sizes and densities: longer problems") {
  const auto long_problem_value = get_long_problem_value();
  const auto medium_problem_values = get_medium_time_problem_values();
  auto tester = get_mixed_sizes_problems();
  for (auto& entry : tester.expected_results) {
    if (entry != long_problem_value &&
        medium_problem_values.count(entry) == 0 && entry < 1000000) {
      // It's one of the short problems, and nort a codeword; skip.
      entry = -1;
    }
  }
  TestParameters params;
  params.timeout = 100000;
  params.total_number_of_edges = 1081;
  params.total_number_of_vertices = 181;
  params.expected_max_total_time_ms = 100000;
  params.expected_min_total_time_ms = 100;
  const auto result = tester.test_all_against_all(params);
  CHECK(result.success_count == 5);
  CHECK(result.failure_count == 0);
  CHECK(result.timeout_count == 95);
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
