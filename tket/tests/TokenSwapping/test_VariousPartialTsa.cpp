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

#include <catch2/catch_test_macros.hpp>

#include "TestUtils/ArchitectureEdgesReimplementation.hpp"
#include "TestUtils/DebugFunctions.hpp"
#include "TestUtils/PartialTsaTesting.hpp"
#include "TestUtils/ProblemGeneration.hpp"
#include "TokenSwapping/CyclesPartialTsa.hpp"
#include "TokenSwapping/RiverFlowPathFinder.hpp"
#include "TokenSwapping/TrivialTSA.hpp"
#include "Utils/RNG.hpp"

using std::vector;

namespace tket {
namespace tsa_internal {
namespace tests {

namespace {
struct Tester {
  vector<std::string> messages_full_trivial_tsa;
  vector<std::string> messages_partial_trivial_tsa;
  vector<std::string> messages_cycles_tsa_0;
  mutable RNG rng;
  mutable TrivialTSA trivial_tsa;
  mutable CyclesPartialTsa cycles_tsa;

  void run_test(
      const ArchitectureMapping& arch_mapping,
      const vector<VertexMapping>& problems, size_t index) const {
    trivial_tsa.set(TrivialTSA::Options::FULL_TSA);
    CHECK(
        run_tests(
            arch_mapping, problems, rng, trivial_tsa,
            RequiredTsaProgress::FULL) == messages_full_trivial_tsa[index]);

    trivial_tsa.set(TrivialTSA::Options::BREAK_AFTER_PROGRESS);
    CHECK(
        run_tests(
            arch_mapping, problems, rng, trivial_tsa,
            RequiredTsaProgress::NONZERO) ==
        messages_partial_trivial_tsa[index]);

    CHECK(
        run_tests(
            arch_mapping, problems, rng, cycles_tsa,
            RequiredTsaProgress::NONE) == messages_cycles_tsa_0[index]);
  }
};

}  // namespace
SCENARIO("Partial TSA: Rings") {
  const vector<std::string> problem_messages{
      "[Ring3: 51582: v3 i1 f100 s1: 100 problems; 135 tokens]",
      "[Ring4: 51481: v4 i1 f100 s1: 100 problems; 178 tokens]",
      "[Ring5: 51644: v5 i1 f100 s1: 100 problems; 224 tokens]",
      "[Ring6: 51528: v6 i1 f100 s1: 100 problems; 270 tokens]",
      "[Ring7: 51496: v7 i1 f100 s1: 100 problems; 318 tokens]",
      "[Ring30: 51633: v30 i1 f100 s1: 100 problems; 1473 tokens]"};

  Tester tester;
  tester.messages_full_trivial_tsa = {
      "[TSA=Trivial FULL PF=RiverFlow\n"
      "135 tokens; 69 total L; 55 swaps.\n"
      "L-decr %: min 100, max 100, av 100.\n"
      "Power %: min 50, max 100, av 82]",

      "[TSA=Trivial FULL PF=RiverFlow\n"
      "178 tokens; 156 total L; 144 swaps.\n"
      "L-decr %: min 100, max 100, av 100.\n"
      "Power %: min 33, max 100, av 69]",

      "[TSA=Trivial FULL PF=RiverFlow\n"
      "224 tokens; 260 total L; 273 swaps.\n"
      "L-decr %: min 100, max 100, av 100.\n"
      "Power %: min 33, max 100, av 59]",

      "[TSA=Trivial FULL PF=RiverFlow\n"
      "270 tokens; 405 total L; 464 swaps.\n"
      "L-decr %: min 100, max 100, av 100.\n"
      "Power %: min 30, max 100, av 52]",

      "[TSA=Trivial FULL PF=RiverFlow\n"
      "318 tokens; 511 total L; 596 swaps.\n"
      "L-decr %: min 100, max 100, av 100.\n"
      "Power %: min 30, max 100, av 49]",

      "[TSA=Trivial FULL PF=RiverFlow\n"
      "1473 tokens; 10908 total L; 16873 swaps.\n"
      "L-decr %: min 100, max 100, av 100.\n"
      "Power %: min 26, max 50, av 36]"};

  tester.messages_partial_trivial_tsa = {
      "[TSA=Trivial NONZERO PF=RiverFlow\n"
      "135 tokens; 69 total L; 49 swaps.\n"
      "L-decr %: min 50, max 100, av 97.\n"
      "Power %: min 50, max 100, av 82]",

      "[TSA=Trivial NONZERO PF=RiverFlow\n"
      "178 tokens; 156 total L; 101 swaps.\n"
      "L-decr %: min 20, max 100, av 80.\n"
      "Power %: min 16, max 100, av 67]",

      "[TSA=Trivial NONZERO PF=RiverFlow\n"
      "224 tokens; 260 total L; 129 swaps.\n"
      "L-decr %: min 12, max 100, av 61.\n"
      "Power %: min 16, max 100, av 58]",

      "[TSA=Trivial NONZERO PF=RiverFlow\n"
      "270 tokens; 405 total L; 186 swaps.\n"
      "L-decr %: min 7, max 100, av 49.\n"
      "Power %: min 8, max 100, av 52]",

      "[TSA=Trivial NONZERO PF=RiverFlow\n"
      "318 tokens; 511 total L; 196 swaps.\n"
      "L-decr %: min 7, max 100, av 39.\n"
      "Power %: min 5, max 100, av 50]",

      "[TSA=Trivial NONZERO PF=RiverFlow\n"
      "1473 tokens; 10908 total L; 273 swaps.\n"
      "L-decr %: min 0, max 50, av 2.\n"
      "Power %: min 1, max 100, av 46]"};

  tester.messages_cycles_tsa_0 = {
      "[TSA=Cycles PF=RiverFlow\n"
      "135 tokens; 69 total L; 55 swaps.\n"
      "L-decr %: min 100, max 100, av 100.\n"
      "Power %: min 50, max 100, av 82]",

      "[TSA=Cycles PF=RiverFlow\n"
      "178 tokens; 156 total L; 119 swaps.\n"
      "L-decr %: min 0, max 100, av 97.\n"
      "Power %: min 0, max 100, av 72]",

      "[TSA=Cycles PF=RiverFlow\n"
      "224 tokens; 260 total L; 194 swaps.\n"
      "L-decr %: min 0, max 100, av 94.\n"
      "Power %: min 0, max 100, av 65]",

      "[TSA=Cycles PF=RiverFlow\n"
      "270 tokens; 405 total L; 294 swaps.\n"
      "L-decr %: min 0, max 100, av 92.\n"
      "Power %: min 0, max 100, av 63]",

      "[TSA=Cycles PF=RiverFlow\n"
      "318 tokens; 511 total L; 357 swaps.\n"
      "L-decr %: min 0, max 100, av 89.\n"
      "Power %: min 0, max 100, av 62]",

      "[TSA=Cycles PF=RiverFlow\n"
      "1473 tokens; 10908 total L; 6344 swaps.\n"
      "L-decr %: min 42, max 100, av 79.\n"
      "Power %: min 50, max 86, av 61]"};

  std::string arch_name;
  const ProblemGenerator00 generator;

  for (size_t index = 0; index < problem_messages.size(); ++index) {
    auto num_vertices = index + 3;
    if (num_vertices == 8) {
      num_vertices = 30;
    }
    const RingArch arch(num_vertices);
    arch_name = "Ring" + std::to_string(num_vertices);

    // OK to reuse RNG, as it's reset before each problem.
    tester.rng.set_seed();
    const auto problems = generator.get_problems(
        arch_name, num_vertices, tester.rng, problem_messages[index]);

    const ArchitectureMapping arch_mapping(arch);
    tester.run_test(arch_mapping, problems, index);
  }
}

SCENARIO("Partial TSA: Square grid") {
  const vector<std::array<unsigned, 3>> grid_parameters = {
      {2, 3, 3}, {5, 5, 3}};
  const vector<std::string> problem_messages{
      "[Grid(2,3,3): 51683: v18 i1 f100 s1: 100 problems; 865 tokens]",
      "[Grid(5,5,3): 51573: v75 i1 f100 s1: 100 problems; 3751 tokens]"};

  Tester tester;
  tester.messages_full_trivial_tsa = {
      "[TSA=Trivial FULL PF=RiverFlow\n"
      "865 tokens; 1921 total L; 2592 swaps.\n"
      "L-decr %: min 100, max 100, av 100.\n"
      "Power %: min 31, max 100, av 41]",

      "[TSA=Trivial FULL PF=RiverFlow\n"
      "3751 tokens; 15297 total L; 23212 swaps.\n"
      "L-decr %: min 100, max 100, av 100.\n"
      "Power %: min 28, max 50, av 36]"};

  tester.messages_partial_trivial_tsa = {
      "[TSA=Trivial NONZERO PF=RiverFlow\n"
      "865 tokens; 1921 total L; 153 swaps.\n"
      "L-decr %: min 2, max 100, av 12.\n"
      "Power %: min 8, max 100, av 48]",

      "[TSA=Trivial NONZERO PF=RiverFlow\n"
      "3751 tokens; 15297 total L; 193 swaps.\n"
      "L-decr %: min 0, max 25, av 1.\n"
      "Power %: min 5, max 100, av 44]"};

  tester.messages_cycles_tsa_0 = {
      "[TSA=Cycles PF=RiverFlow\n"
      "865 tokens; 1921 total L; 1425 swaps.\n"
      "L-decr %: min 60, max 100, av 95.\n"
      "Power %: min 46, max 100, av 61]",

      "[TSA=Cycles PF=RiverFlow\n"
      "3751 tokens; 15297 total L; 11464 swaps.\n"
      "L-decr %: min 83, max 100, av 95.\n"
      "Power %: min 50, max 79, av 59]"};

  const ProblemGenerator00 generator;

  for (size_t index = 0; index < grid_parameters.size(); ++index) {
    const auto& parameters = grid_parameters[index];

    const auto edges =
        get_square_grid_edges(parameters[0], parameters[1], parameters[2]);
    const Architecture arch(edges);
    const ArchitectureMapping arch_mapping(arch, edges);

    std::stringstream ss;
    ss << "Grid(" << parameters[0] << "," << parameters[1] << ","
       << parameters[2] << ")";

    tester.rng.set_seed();
    const auto problems = generator.get_problems(
        ss.str(), arch.n_nodes(), tester.rng, problem_messages[index]);

    tester.run_test(arch_mapping, problems, index);
  }
}

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
