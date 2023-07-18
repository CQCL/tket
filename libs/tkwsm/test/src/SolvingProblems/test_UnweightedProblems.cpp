// Copyright 2019-2023 Cambridge Quantum Computing
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

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <tkwsm/Common/GeneralUtils.hpp>
#include <utility>

#include "../TestUtils/CheckedSolution.hpp"
#include "../TestUtils/ResumedSolutionChecker.hpp"
#include "../TestUtils/TestSettings.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {
namespace {

// Try to embed graphs from the first sequence
// into graphs from the second sequence,
// recording the result in a string (for easy copy/paste).
struct EmbedGraphSequences {
  long long total_time_ms;

  // Simply use 0 for no embedding, 1 for an embedding,
  // * for timeout, and letters for errors.
  std::string result;

  EmbedGraphSequences(
      const std::vector<GraphEdgeWeights>& graph_sequence1,
      const std::vector<GraphEdgeWeights>& graph_sequence2, unsigned timeout_ms,
      const std::string& expected_result)
      : total_time_ms(0) {
    CheckedSolution::Statistics statistics(
        "unweighted problems; embedding graph sequences");

    MainSolverParameters solver_params(timeout_ms);
    solver_params.terminate_with_first_full_solution = true;

    CheckedSolution::ProblemInformation info;
    std::stringstream ss;
    unsigned result_index = 0;
    ResumedSolutionChecker resumption_checker;

    for (unsigned index1 = 0; index1 < graph_sequence1.size(); ++index1) {
      const auto& pattern_graph = graph_sequence1[index1];
      for (unsigned index2 = 0; index2 < graph_sequence2.size(); ++index2) {
        const auto& target_graph = graph_sequence2[index2];
        const bool timeout_expected = result_index < expected_result.size() &&
                                      expected_result.at(result_index) == '*';

        if (result_index % 8 == 0) {
          TestSettings::get().os << "\n### RI=" << result_index << ": ";
        }
        ++result_index;

        if (timeout_expected) {
          // To save time, don't bother trying to solve
          // known hard problems.
          ss << "*";
          continue;
        }

        const auto search_time_before = statistics.total_search_time_ms;
        const CheckedSolution checked_solution(
            pattern_graph, target_graph, info, solver_params, statistics);

        resumption_checker.check(
            checked_solution, pattern_graph, target_graph, solver_params);

        if (checked_solution.scalar_product == pattern_graph.size()) {
          ss << "1";
          continue;
        }
        if (checked_solution.scalar_product == 0) {
          if (checked_solution.finished) {
            ss << "0";
          } else {
            // Timed out.
            ss << "*";
          }
        } else {
          // Error: wrong scalar product!
          ss << "X";
        }
      }
    }
    result = ss.str();
    if (!expected_result.empty()) {
      CHECK(expected_result.size() == result.size());
      CHECK(result_index == result.size());
    }
    total_time_ms =
        statistics.total_init_time_ms + statistics.total_search_time_ms;

    statistics.finish();
  }
};

}  // namespace

typedef std::mt19937_64 RNG_64;

// Use 16 random bits as the sorting key,
// to get approx uniform distribution of permutations.
template <class T>
static void reorder(
    RNG_64& rng, std::vector<std::pair<std::uint_fast16_t, T>>& data) {
  std::uint64_t bits = 0;
  for (auto& entry : data) {
    if (bits == 0) {
      bits = rng();
    }
    entry.first = bits & 0xffff;
    bits >>= 16;
  }
  std::sort(data.begin(), data.end());
}

// Add edges gradually to a graph,
// to get a sequence of graphs each one of which embeds in the next,
// but randomly relabelling the vertices to make it harder for
// the solver.
static std::vector<GraphEdgeWeights> get_increasing_graph_sequence(
    unsigned number_of_vertices, unsigned num_entries, RNG_64& rng) {
  std::vector<std::pair<std::uint_fast16_t, EdgeWSM>> edges_data;
  edges_data.reserve((number_of_vertices * (number_of_vertices - 1)) / 2);
  std::vector<std::pair<std::uint_fast16_t, unsigned>> new_labels(
      number_of_vertices);

  for (unsigned ii = 0; ii < number_of_vertices; ++ii) {
    new_labels[ii].second = ii;
    for (unsigned jj = ii + 1; jj < number_of_vertices; ++jj) {
      edges_data.emplace_back();
      edges_data.back().second = get_edge(ii, jj);
    }
  }
  reorder(rng, edges_data);
  const unsigned num_edges_increment = edges_data.size() / (num_entries + 1);
  REQUIRE(num_edges_increment > 0);
  REQUIRE(num_edges_increment * num_entries < edges_data.size());

  // Now create the increasing graphs.
  std::vector<GraphEdgeWeights> graph_data;
  for (unsigned multiplier = 1; multiplier <= num_entries; ++multiplier) {
    const unsigned num_edges = num_edges_increment * multiplier;
    reorder(rng, new_labels);
    // Now add the edges.
    graph_data.emplace_back();
    for (unsigned ee = 0; ee < num_edges; ++ee) {
      const auto& edge = edges_data[ee].second;
      const auto new_edge = get_edge(
          new_labels[edge.first].second, new_labels[edge.second].second);

      // WeightWSM 1 for every edge.
      graph_data.back()[new_edge] = 1;
    }
    if (graph_data.size() >= num_entries) {
      break;
    }
  }
  return graph_data;
}

// A string like "111111110111..." records the results
// of trying to embed graph P(i) into T(j).
// The graphs come from increasing sequences,
// so there should be a cutoff point dividing 0 and 1.
static void check_monotonic_embedding_property(
    const std::string& str_result, unsigned n_target_graphs,
    bool same_sequence) {
  const unsigned n_pattern_graphs = str_result.size() / n_target_graphs;
  CHECK(n_pattern_graphs * n_target_graphs == str_result.size());
  if (same_sequence) {
    CHECK(n_pattern_graphs == n_target_graphs);
  }
  unsigned index = 0;
  unsigned embed_count = 0;
  unsigned nonembed_count = 0;

  // The pattern graphs and target graphs are both increasing.
  // In each target graph block, it should START at 0 and switch over to 1.

  for (unsigned p_index = 0; p_index < n_pattern_graphs; ++p_index) {
    const auto previous_embed_count = embed_count;
    const auto previous_nonembed_count = nonembed_count;
    embed_count = 0;
    nonembed_count = 0;

    for (unsigned t_index = 0; t_index < n_target_graphs; ++t_index) {
      auto& symbol = str_result.at(index);
      ++index;
      if (symbol == '1') {
        ++embed_count;
        if (same_sequence) {
          // If it happens to be the same increasing sequence
          // in the source and target, clearly this must hold.
          CHECK(t_index >= p_index);
        }
        continue;
      }
      if (symbol == '0') {
        ++nonembed_count;
        CHECK(embed_count == 0);
        if (same_sequence) {
          CHECK(t_index < p_index);
        }
        continue;
      }
    }

    if (embed_count + nonembed_count ==
            previous_embed_count + previous_nonembed_count &&
        embed_count + nonembed_count == n_target_graphs) {
      // No timeouts. The number of embeddings must be DECREASING,
      // because the pattern graphs are getting bigger.
      CHECK(embed_count <= previous_embed_count);
    }
  }
}

static const std::vector<std::string>& get_expected_results_ref() {
  static const std::vector<std::string> expected_results{
      "1111111101111111001111110001111100001111000001110000001100000001",
      "1111111111111111111111110111111100111111000111110001111100011111",
      "1111111111111111111111111111111101111111011111110011111100111111",
      "1111111111111111111111111111111111111111011111110011111100111111",
      "1111111111111111111111111111111111111111111111110111111101111111",
      "0000000100000000000000000000000000000000000000000000000000000000",
      "1111111101111111001111110001111100001111000001110000001100000001",
      "1111111101111111001111110000111100000111000000110000001100000001",
      "1111111101111111001111110001111100001111000011110000001100000011",
      "1111111111111111011111110011111100011111000111110000011100000011",
      "0000000000000000000000000000000000000000000000000000000000000000",
      "0011111100000000000000000000000000000000000000000000000000000000",
      "1111111101111111001111110001111100001111000001110000001100000001",
      "0111111100111111000111110000111100000111000000110000000100000001",
      "1111111101111111001111110001111100001111000001110000001100000001",
      "0000000000000000000000000000000000000000000000000000000000000000",
      "0000000000000000000000000000000000000000000000000000000000000000",
      "0011111100000000000000000000000000000000000000000000000000000000",
      "1111111101111111001111110001111100001111000001110000001100000001",
      "1111111100111111000011110000111100000111000000110000001100000001",
      "0000000000000000000000000000000000000000000000000000000000000000",
      "0000000000000000000000000000000000000000000000000000000000000000",
      "0000000000000000000000000000000000000000000000000000000000000000",
      "0000000000000000000000000000000000000000000000000000000000000000",
      "1111111101111111001111110001111100001111000001110000001100000001"};
  return expected_results;
}

static std::vector<std::vector<GraphEdgeWeights>>
get_list_of_increasing_graph_sequences() {
  std::vector<std::vector<GraphEdgeWeights>> list_of_increasing_graph_sequences;
  const unsigned num_entries = 8;
  unsigned number_of_vertices = 3;
  RNG_64 rng;
  for (int count = 0; count < 5; ++count) {
    number_of_vertices += 3;
    list_of_increasing_graph_sequences.emplace_back(
        get_increasing_graph_sequence(number_of_vertices, num_entries, rng));
  }
  {
    // A crude check that the test data hasn't changed,
    // and is identical across platforms.
    CHECK(rng() == 0x3c5c9fe803f69af3);
    const auto& final_list = list_of_increasing_graph_sequences.back();
    CHECK(final_list.size() == 8);
    const auto& final_graph = final_list.back();
    CHECK(final_graph.size() == 136);

    // Check an edge in the middle...
    const auto& middle_graph = final_list[final_list.size() / 2];
    CHECK(middle_graph.size() == 85);
    std::size_t counter = middle_graph.size() / 2;
    for (const auto& edge_weight_pair : middle_graph) {
      CHECK(edge_weight_pair.second == 1);
      if (counter != 0) {
        --counter;
        continue;
      }
      CHECK(edge_weight_pair.first.first == 5);
      CHECK(edge_weight_pair.first.second == 6);
      break;
    }
  }
  return list_of_increasing_graph_sequences;
}

static const std::vector<std::vector<GraphEdgeWeights>>&
get_list_of_increasing_graph_sequences_ref() {
  static const auto list_of_increasing_graph_sequences =
      get_list_of_increasing_graph_sequences();
  return list_of_increasing_graph_sequences;
}

static std::set<std::pair<unsigned, unsigned>> get_longer_ij_pair_tests() {
  return {{2, 3}, {2, 4}, {3, 3}, {3, 4}, {4, 4}};
}

static void test(bool short_test) {
  const unsigned num_entries = 8;
  const unsigned timeout_ms = short_test ? 1000 : 10000;
  const auto& expected_results = get_expected_results_ref();

  long long total_time_ms = 0;
  unsigned expected_str_index = 0;
  const auto& list_of_increasing_graph_sequences =
      get_list_of_increasing_graph_sequences_ref();
  const std::string short_test_str = short_test ? "SHORT" : "LONG";
  TestSettings::get().os << "\n\nRunning unweighted probs: " << short_test_str;
  const auto longer_test_pairs = get_longer_ij_pair_tests();
  unsigned pair_count = 0;

  for (unsigned ii = 0; ii < list_of_increasing_graph_sequences.size(); ++ii) {
    for (unsigned jj = 0; jj < list_of_increasing_graph_sequences.size();
         ++jj) {
      TestSettings::get().os << "\ni=" << ii << ", j=" << jj << " : ";
      const bool should_test =
          (longer_test_pairs.count(std::make_pair(ii, jj)) == 0) == short_test;

      if (!should_test) {
        TestSettings::get().os << "SKIPPED";
        ++expected_str_index;
        continue;
      }
      const EmbedGraphSequences embedding_tester(
          list_of_increasing_graph_sequences[ii],
          list_of_increasing_graph_sequences[jj], timeout_ms,
          expected_results.at(expected_str_index));

      CHECK(embedding_tester.result == expected_results.at(expected_str_index));
      ++expected_str_index;
      ++pair_count;
      TestSettings::get().os << "\n@@@@@ (" << ii << "," << jj << ") took time "
                             << embedding_tester.total_time_ms;
      total_time_ms += embedding_tester.total_time_ms;
      check_monotonic_embedding_property(
          embedding_tester.result, num_entries, ii == jj);
    }
  }
  CHECK(expected_str_index == expected_results.size());
  TestSettings::get().os << "\n::::END: all unweighted probs for " << pair_count
                         << " (i,j) pairs, " << short_test_str << " tests; "
                         << total_time_ms << " ms.\n";
}

SCENARIO("Increasing graph sequences: short tests") { test(true); }
SCENARIO("Increasing graph sequences: long tests") { test(false); }

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
