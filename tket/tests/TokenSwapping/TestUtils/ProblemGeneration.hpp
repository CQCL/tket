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

#pragma once

#include "Architecture/Architecture.hpp"
#include "TokenSwapping/VertexMappingFunctions.hpp"
#include "Utils/RNG.hpp"

namespace tket {
namespace tsa_internal {
namespace tests {

struct TSProblemParameters00 {
  // How many tokens are there, as a percentage of the number of vertices?
  // Will still work if above 100, just gets truncated to 100%.
  unsigned token_density_percentage;

  // For very small graphs, ensure a minimum number of tokens.
  unsigned min_number_of_tokens;
  unsigned max_number_of_tokens;

  TSProblemParameters00();

  // Using the above problem parameters
  VertexMapping get_problem(RNG& rng, unsigned number_of_vertices) const;
};

// Given an architecture, generate various test problems
// with varying numbers of tokens.
struct ProblemGenerator00 {
  unsigned init_token_density_percentage;
  unsigned final_percentage;
  unsigned step;

  ProblemGenerator00();

  std::vector<VertexMapping> get_problems(
      const std::string& arch_name, unsigned number_of_vertices, RNG& rng,
      // It will calculate a short summary string of the problems
      // and check against this string; this helps to detect
      // accidentally changed parameters/generation algorithms
      // leading to different tests.
      const std::string& expected_summary) const;
};

struct RandomTreeGenerator00 {
  // Every finite tree must have a leaf!
  // So, some vertices will end up being leaves (having no children),
  // even if the min is nonzero.
  unsigned min_number_of_children;
  unsigned max_number_of_children;
  unsigned approx_number_of_vertices;
  mutable std::vector<unsigned> work_vector;

  RandomTreeGenerator00();

  // Creates the edges of a random tree with vertices {0,1,2,...} with
  // vertex 0 being the root.
  // It might not find exactly the requested number of vertices.
  // Note that (number of vertices) == (number of edges+1), for a tree.
  std::vector<std::pair<unsigned, unsigned>> get_tree_edges(RNG& rng) const;
};

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
