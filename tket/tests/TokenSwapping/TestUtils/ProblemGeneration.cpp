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

#include "ProblemGeneration.hpp"

#include <catch2/catch_test_macros.hpp>

#include "GetRandomSet.hpp"
#include "TokenSwapping/GeneralFunctions.hpp"

using std::vector;

namespace tket {
namespace tsa_internal {
namespace tests {

TSProblemParameters00::TSProblemParameters00()
    : token_density_percentage(10),
      min_number_of_tokens(1),
      max_number_of_tokens(10000) {}

VertexMapping TSProblemParameters00::get_problem(
    RNG& rng, unsigned number_of_vertices) const {
  unsigned number_of_tokens =
      (token_density_percentage * number_of_vertices) / 100;
  number_of_tokens = std::max(number_of_tokens, min_number_of_tokens);
  number_of_tokens = std::min(number_of_tokens, number_of_vertices);
  number_of_tokens = std::min(number_of_tokens, max_number_of_tokens);

  VertexMapping vertex_mapping;
  const auto tokens = get_random_set(rng, number_of_tokens, number_of_vertices);
  const auto targets_set =
      get_random_set(rng, number_of_tokens, number_of_vertices);
  REQUIRE(tokens.size() == number_of_tokens);
  REQUIRE(targets_set.size() == number_of_tokens);
  vector<size_t> targets{targets_set.cbegin(), targets_set.cend()};
  for (auto token : tokens) {
    vertex_mapping[token] = rng.get_and_remove_element(targets);
  }
  REQUIRE(targets.empty());
  REQUIRE(vertex_mapping.size() == number_of_tokens);
  return vertex_mapping;
}

ProblemGenerator00::ProblemGenerator00()
    : init_token_density_percentage(1), final_percentage(100), step(1) {}

vector<VertexMapping> ProblemGenerator00::get_problems(
    const std::string& arch_name, unsigned number_of_vertices, RNG& rng,
    // It will calculate a short summary string of the problems
    // and check against this string; this helps to detect
    // accidentally changed parameters/generation algorithms
    // leading to different tests.
    const std::string& expected_summary) const {
  REQUIRE(step > 0);

  TSProblemParameters00 params;
  vector<VertexMapping> vertex_mappings;

  // This will probably detect if the rng changes, or has different seed
  auto code = rng.get_size_t(255);
  unsigned tokens_count = 0;
  for (params.token_density_percentage = init_token_density_percentage;
       params.token_density_percentage <= final_percentage;
       params.token_density_percentage += step) {
    vertex_mappings.push_back(params.get_problem(rng, number_of_vertices));
    tokens_count += vertex_mappings.back().size();
  }
  code = (code << 8) + rng.get_size_t(255);
  std::stringstream ss;
  ss << "[" << arch_name << ": " << code << ": v" << number_of_vertices << " i"
     << init_token_density_percentage << " f" << final_percentage << " s"
     << step << ": " << vertex_mappings.size() << " problems; " << tokens_count
     << " tokens]";
  CHECK(ss.str() == expected_summary);
  return vertex_mappings;
}

RandomTreeGenerator00::RandomTreeGenerator00()
    : min_number_of_children(1),
      max_number_of_children(3),
      approx_number_of_vertices(10) {}

// Creates the edges of a random tree with vertex 0 being the root.
vector<std::pair<unsigned, unsigned>> RandomTreeGenerator00::get_tree_edges(
    RNG& rng) const {
  REQUIRE(max_number_of_children > min_number_of_children);
  REQUIRE(max_number_of_children > 1);
  REQUIRE(approx_number_of_vertices >= 3);
  // The vertices awaiting child nodes to be assigned.
  work_vector.resize(1);
  work_vector[0] = 0;

  vector<std::pair<unsigned, unsigned>> edges;
  for (auto infinite_loop_guard = 100 + 100 * approx_number_of_vertices;
       infinite_loop_guard > 0; --infinite_loop_guard) {
    const auto number_of_children =
        rng.get_size_t(min_number_of_children, max_number_of_children);
    const unsigned node = rng.get_and_remove_element(work_vector);
    for (unsigned ii = 0; ii < number_of_children; ++ii) {
      const unsigned new_vertex = edges.size() + 1;
      work_vector.push_back(new_vertex);
      edges.emplace_back(node, new_vertex);
      if (edges.size() + 1 >= approx_number_of_vertices) {
        return edges;
      }
    }
    if (work_vector.empty()) {
      return edges;
    }
  }
  REQUIRE(false);
  return edges;
}

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
