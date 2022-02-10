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
#include <istream>
#include <map>
#include <string>
#include <utility>

#include "../TestUtils/CheckedSolution.hpp"
#include "../TestUtils/TestSettings.hpp"
#include "WeightSubgrMono/Common/GeneralUtils.hpp"

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
  bool test_trivially_impossible_embeddings = false;
  bool recalculate_known_timeouts = false;
  unsigned timeout = 1000;
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

// To save time, do not actually compute the solution for
// the (i,j) pair; instead, simply take the given value
// as if it had been computed.
typedef std::map<std::pair<unsigned, unsigned>, int> OverwriteValues;

// expected_results should FIRST list the hashes of the graphs with weights,
// THEN the scalar products S in order: S=0 means no solution,
// S>0 means an optimal solution with S was found, -1 means a timeout.
static TestResult test_all_against_all(
    const std::vector<std::string>& codes,
    const std::vector<WeightWSM>& weights,
    const ResultsSummary& expected_results, const TestParameters& params,
    const OverwriteValues& shortcut_values = {}) {
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
  os << "\n\n########### generated " << graphs.size()
     << " random graphs, total " << result.total_edges << " edges, "
     << result.total_verts
     << " vertices\n#### now testing all against all, timeout "
     << params.timeout;
  CheckedSolution::Statistics statistics;
  const MainSolver::Parameters solver_params(params.timeout);

  const auto update_calc_results =
      [&calc_results](const CheckedSolution& checked_solution) {
        if (!checked_solution.finished) {
          // A timeout.
          calc_results.push_back(-1);
          return;
        }
        if (checked_solution.complete_solution_weight) {
          calc_results.push_back(
              checked_solution.complete_solution_weight.value());
        } else {
          calc_results.push_back(0);
        }
      };

  for (unsigned ii = 0; ii < graphs.size(); ++ii) {
    for (unsigned jj = 0; jj < graphs.size(); ++jj) {
      const auto value_optional =
          get_optional_value(shortcut_values, std::make_pair(ii, jj));
      if (value_optional) {
        ++statistics.success_count;
        calc_results.push_back(value_optional.value());
        continue;
      }

      if ((num_vertices[ii] > num_vertices[jj] ||
           graphs[ii].size() > graphs[jj].size()) &&
          !params.test_trivially_impossible_embeddings) {
        ++statistics.success_count;
        calc_results.push_back(0);
        continue;
      }
      if (!params.recalculate_known_timeouts &&
          calc_results.size() < expected_results.size() &&
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
  os << "\n#### FIN: total time " << result.total_time_ms << " ms. ";
  CHECK(expected_results == calc_results);
  return result;
}

SCENARIO("embedding random graphs - smaller graphs, small weights") {
  const std::vector<std::string> codes{"5 8 111",  "5 9 12211", "6 10 13311",
                                       "7 10 222", "7 15 333",  "8 16 1111",
                                       "8 20 444", "10 20 333"};
  const std::vector<WeightWSM> weights{1, 2, 3, 8};
  const ResultsSummary expected_results{
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

  const auto result =
      test_all_against_all(codes, weights, expected_results, params);

  CHECK(result.total_time_ms < 10 * 6);
  CHECK(result.success_count == 64);
  CHECK(result.failure_count == 0);
  CHECK(result.timeout_count == 0);
  CHECK(result.total_edges == 108);
  CHECK(result.total_verts == 56);
}

SCENARIO("embedding random graphs - medium graphs, small weights") {
  const std::vector<std::string> codes{
      "10 20 1111", "10 30 2222",  "11 20 3333", "11 40 4444",
      "15 30 5555", "16 50, 6666", "17 60 7777", "18 70 888"};
  const std::vector<WeightWSM> weights{1, 2, 3, 8};
  const ResultsSummary expected_results{
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
  params.timeout = 5000;

  const auto result =
      test_all_against_all(codes, weights, expected_results, params);

  CHECK(result.total_time_ms < 10 * 728);
  CHECK(result.total_time_ms > 728 / 10);
  CHECK(result.success_count == 64);
  CHECK(result.failure_count == 0);
  CHECK(result.timeout_count == 0);
  CHECK(result.total_edges == 320);
  CHECK(result.total_verts == 108);
}

SCENARIO("embedding random graphs - large graphs, small weights") {
  const std::vector<std::string> codes{
      "20 50 1111",  "22 80 2222",    "25 120 3333", "25 200 4444",
      "30 200 5555", "32 300 6666",   "35 300 7777", "40 500 8888",
      "50 500 9999", "55 1000 101010"};
  const std::vector<WeightWSM> weights{1, 2, 3, 8};

  // 5000 ms timeout...
  const ResultsSummary expected_results{
      460517071,  255540664,  811304662, 581415081, 853453591, 367941120,
      1072813581, 1006422874, 309411091, 971368384, 1261,      0,
      619,        -1,         -1,        -1,        -1,        -1,
      -1,         -1,         0,         1732,      0,         -1,
      0,          -1,         -1,        -1,        -1,        -1,
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

  OverwriteValues overwrite_values;
  unsigned expected_time_ms = 8000;
  if (!TestSettings::get().run_slow_tests) {
    overwrite_values[std::make_pair<unsigned>(1, 4)] = 0;
    expected_time_ms -= 4000;
  }

  TestParameters params;
  params.timeout = 5000;

  const auto result = test_all_against_all(
      codes, weights, expected_results, params, overwrite_values);

  CHECK(result.total_time_ms < 10 * expected_time_ms);
  CHECK(result.total_time_ms > expected_time_ms / 10);
  CHECK(result.success_count == 70);
  CHECK(result.failure_count == 0);
  CHECK(result.timeout_count == 30);
  CHECK(result.total_edges == 3250);
  CHECK(result.total_verts == 334);
}

SCENARIO("embedding random graphs - mixed sizes and densities") {
  const std::vector<std::string> codes{
      "5 7 111",     "6 14 222",    "10 20 1111",  "10 40 3333",
      "20 50 4444",  "20 100 5555", "20 150 6666", "30 100 7777",
      "30 200 8888", "30 400 9999"};
  const std::vector<WeightWSM> weights{1, 2, 5, 20};

  // timeout 5000
  const ResultsSummary expected_results{
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

  TestParameters params;
  params.timeout = 5000;

  OverwriteValues overwrite_values;
  unsigned expected_time_ms = 7000;
  if (!TestSettings::get().run_slow_tests) {
    overwrite_values[std::make_pair<unsigned>(4, 5)] = 2461;
    expected_time_ms -= 3300;
  }
  const auto result = test_all_against_all(
      codes, weights, expected_results, params, overwrite_values);

  CHECK(result.total_time_ms < 10 * expected_time_ms);
  CHECK(result.total_time_ms > expected_time_ms / 10);
  CHECK(result.success_count == 85);
  CHECK(result.failure_count == 0);
  CHECK(result.timeout_count == 15);
  CHECK(result.total_edges == 1081);
  CHECK(result.total_verts == 181);
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
