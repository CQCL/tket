#include <catch2/catch.hpp>

#include "TestUtils/FullTsaTesting.hpp"
#include "TestUtils/ProblemGeneration.hpp"
#include "TokenSwapping/HybridTsa00.hpp"

using std::vector;

namespace tket {
namespace tsa_internal {
namespace tests {

namespace {
struct FullTester {
  FullTsaTesting results;
  FullTsaTesting trivial_results;
  HybridTsa00 full_tsa;
  TrivialTSA trivial_tsa;
  RNG rng;
  ProblemGenerator00 generator;
  std::string test_name;

  void add_problems(
      const Architecture& arch, const std::string& arch_name,
      const std::string& problem_message) {
    rng.set_seed();
    const auto problems =
        generator.get_problems(arch_name, arch, rng, problem_message);

    // OK to reuse RNG, as it's reset before each problem.
    results.add_problems(arch, problems, test_name, rng, full_tsa);

    trivial_tsa.set(TrivialTSA::Options::FULL_TSA);
    trivial_results.add_problems(arch, problems, test_name, rng, trivial_tsa);
  }
};
}  // namespace

SCENARIO("Full TSA: stars") {
  const vector<std::string> problem_messages{
      "[Star3: 51481: v4 i1 f100 s1: 100 problems; 178 tokens]",
      "[Star5: 51528: v6 i1 f100 s1: 100 problems; 270 tokens]",
      "[Star10: 51662: v11 i1 f100 s1: 100 problems; 515 tokens]",
      "[Star20: 51494: v21 i1 f100 s1: 100 problems; 1015 tokens]"};
  const vector<size_t> num_spokes{3, 5, 10, 20};
  FullTester tester;
  tester.test_name = "Stars";
  std::string arch_name;
  vector<std::pair<unsigned, unsigned>> edges;

  for (size_t index = 0; index < problem_messages.size(); ++index) {
    arch_name = "Star" + std::to_string(num_spokes[index]);
    edges.clear();
    for (unsigned vv = 1; vv <= num_spokes[index]; ++vv) {
      edges.emplace_back(0, vv);
    }
    const Architecture arch(edges);
    tester.add_problems(arch, arch_name, problem_messages[index]);
  }
  CHECK(
      tester.results.str() ==
      "[Stars:HybridTSA_00: 400 probs; 1978 toks; 1623 tot.lb]\n"
      "[Total swaps: 2632 2588 2550 2539 2539 2550]\n"
      "[Winners: joint: 360 381 392 400 400 392  undisputed: 0 0 0 0 0 0]");

  CHECK(
      tester.trivial_results.str() ==
      "[Stars:Trivial: 400 probs; 1978 toks; 1623 tot.lb]\n"
      "[Total swaps: 3968 3804 3088 3088 3088 3088]\n"
      "[Winners: joint: 247 271 400 400 400 400  undisputed: 0 0 0 0 0 0]");
}

SCENARIO("Full TSA: wheels") {
  const vector<std::string> problem_messages{
      "[Wheel3: 51481: v4 i1 f100 s1: 100 problems; 178 tokens]",
      "[Wheel5: 51528: v6 i1 f100 s1: 100 problems; 270 tokens]",
      "[Wheel10: 51662: v11 i1 f100 s1: 100 problems; 515 tokens]",
      "[Wheel20: 51494: v21 i1 f100 s1: 100 problems; 1015 tokens]"};

  const vector<size_t> num_spokes{3, 5, 10, 20};
  FullTester tester;
  tester.test_name = "Wheels";
  std::string arch_name;
  vector<std::pair<unsigned, unsigned>> edges;

  for (size_t index = 0; index < problem_messages.size(); ++index) {
    arch_name = "Wheel" + std::to_string(num_spokes[index]);
    edges.clear();
    for (unsigned vv = 1; vv <= num_spokes[index]; ++vv) {
      edges.emplace_back(0, vv);
      if (vv == num_spokes[index]) {
        edges.emplace_back(vv, 1);
      } else {
        edges.emplace_back(vv, vv + 1);
      }
    }
    const Architecture arch(edges);
    tester.add_problems(arch, arch_name, problem_messages[index]);
  }
  CHECK(
      tester.results.str() ==
      "[Wheels:HybridTSA_00: 400 probs; 1978 toks; 1533 tot.lb]\n"
      "[Total swaps: 2482 2462 2430 2422 2422 2430]\n"
      "[Winners: joint: 374 384 395 400 400 395  undisputed: 0 0 0 0 0 0]");

  CHECK(
      tester.trivial_results.str() ==
      "[Wheels:Trivial: 400 probs; 1978 toks; 1533 tot.lb]\n"
      "[Total swaps: 3510 3410 2818 2818 2818 2818]\n"
      "[Winners: joint: 283 291 400 400 400 400  undisputed: 0 0 0 0 0 0]");
}

SCENARIO("Full TSA: Rings") {
  const vector<std::string> problem_messages{
      "[Ring3: 51582: v3 i1 f100 s1: 100 problems; 135 tokens]",
      "[Ring5: 51644: v5 i1 f100 s1: 100 problems; 224 tokens]",
      "[Ring10: 51634: v10 i1 f100 s1: 100 problems; 469 tokens]",
      "[Ring20: 51498: v20 i1 f100 s1: 100 problems; 974 tokens]"};
  const vector<size_t> num_vertices{3, 5, 10, 20};
  FullTester tester;
  tester.test_name = "Rings";
  std::string arch_name;

  for (size_t index = 0; index < problem_messages.size(); ++index) {
    const RingArch arch(num_vertices[index]);
    arch_name = "Ring" + std::to_string(num_vertices[index]);
    tester.add_problems(arch, arch_name, problem_messages[index]);
  }
  CHECK(
      tester.results.str() ==
      "[Rings:HybridTSA_00: 400 probs; 1802 toks; 3193 tot.lb]\n"
      "[Total swaps: 6302 5942 5118 5115 5113 5118]\n"
      "[Winners: joint: 292 328 399 399 400 399  undisputed: 0 0 0 0 1 0]");

  CHECK(
      tester.trivial_results.str() ==
      "[Rings:Trivial: 400 probs; 1802 toks; 3193 tot.lb]\n"
      "[Total swaps: 8922 8580 5104 5087 5079 5104]\n"
      "[Winners: joint: 231 252 394 397 400 394  undisputed: 0 0 0 0 3 0]");
}

SCENARIO("Full TSA: fully connected") {
  const vector<std::string> problem_messages{
      "[K3: 51582: v3 i1 f100 s1: 100 problems; 135 tokens]",
      "[K5: 51644: v5 i1 f100 s1: 100 problems; 224 tokens]",
      "[K10: 51634: v10 i1 f100 s1: 100 problems; 469 tokens]",
      "[K20: 51498: v20 i1 f100 s1: 100 problems; 974 tokens]"};

  const vector<size_t> num_vertices{3, 5, 10, 20};
  FullTester tester;
  tester.test_name = "FullyConn";
  std::string arch_name;

  for (size_t index = 0; index < problem_messages.size(); ++index) {
    const FullyConnected arch(num_vertices[index]);
    arch_name = "K" + std::to_string(num_vertices[index]);
    tester.add_problems(arch, arch_name, problem_messages[index]);
  }
  CHECK(
      tester.results.str() ==
      "[FullyConn:HybridTSA_00: 400 probs; 1802 toks; 867 tot.lb]\n"
      "[Total swaps: 1435 1435 1435 1435 1435 1435]\n"
      "[Winners: joint: 400 400 400 400 400 400  undisputed: 0 0 0 0 0 0]");

  CHECK(
      tester.trivial_results.str() ==
      "[FullyConn:Trivial: 400 probs; 1802 toks; 867 tot.lb]\n"
      "[Total swaps: 1435 1435 1435 1435 1435 1435]\n"
      "[Winners: joint: 400 400 400 400 400 400  undisputed: 0 0 0 0 0 0]");
}

SCENARIO("Full TSA: Square Grids") {
  const vector<std::array<unsigned, 3>> grid_parameters = {
      {2, 2, 2}, {3, 4, 4}};
  const vector<std::string> problem_messages{
      "[Grid(2,2,2): 51480: v8 i1 f100 s1: 100 problems; 368 tokens]",
      "[Grid(3,4,4): 51492: v48 i1 f100 s1: 100 problems; 2378 tokens]"};

  FullTester tester;
  tester.test_name = "Square grids";

  for (size_t index = 0; index < grid_parameters.size(); ++index) {
    const auto& parameters = grid_parameters[index];
    const SquareGrid arch(parameters[0], parameters[1], parameters[2]);
    std::stringstream ss;
    ss << "Grid(" << parameters[0] << "," << parameters[1] << ","
       << parameters[2] << ")";

    tester.add_problems(arch, ss.str(), problem_messages[index]);
  }
  CHECK(
      tester.results.str() ==
      "[Square grids:HybridTSA_00: 200 probs; 2746 toks; 4323 tot.lb]\n"
      "[Total swaps: 7083 7015 6863 6846 6842 6863]\n"
      "[Winners: joint: 148 163 188 198 200 188  undisputed: 0 0 0 0 2 0]");

  CHECK(
      tester.trivial_results.str() ==
      "[Square grids:Trivial: 200 probs; 2746 toks; 4323 tot.lb]\n"
      "[Total swaps: 12364 12208 9114 9039 8933 9114]\n"
      "[Winners: joint: 85 91 152 177 200 152  undisputed: 0 0 0 0 23 0]");
}

SCENARIO("Full TSA: Random trees") {
  RandomTreeGenerator00 tree_generator;
  FullTester tester;

  const vector<std::string> problem_messages{
      "[Tree0: 51644: v5 i1 f100 s1: 100 problems; 224 tokens]",
      "[Tree1: 51517: v16 i1 f100 s1: 100 problems; 766 tokens]",
      "[Tree2: 51481: v24 i1 f100 s1: 100 problems; 1168 tokens]"};
  tester.test_name = "Trees";
  std::string arch_name;

  for (size_t index = 0; index < problem_messages.size(); ++index) {
    tree_generator.min_number_of_children = index;
    tree_generator.max_number_of_children = 2 + 2 * index;
    tree_generator.approx_number_of_vertices =
        4 * tree_generator.max_number_of_children;

    const auto edges = tree_generator.get_tree_edges(tester.rng);
    const Architecture arch(edges);
    REQUIRE(arch.n_uids() == edges.size() + 1);
    arch_name = "Tree" + std::to_string(index);
    tester.add_problems(arch, arch_name, problem_messages[index]);
  }
  CHECK(
      tester.results.str() ==
      "[Trees:HybridTSA_00: 300 probs; 2158 toks; 2963 tot.lb]\n"
      "[Total swaps: 5216 5132 4844 4828 4817 4844]\n"
      "[Winners: joint: 227 251 286 296 300 286  undisputed: 0 0 0 0 4 0]");

  CHECK(
      tester.trivial_results.str() ==
      "[Trees:Trivial: 300 probs; 2158 toks; 2963 tot.lb]\n"
      "[Total swaps: 8128 7886 5592 5570 5563 5600]\n"
      "[Winners: joint: 128 148 282 297 300 280  undisputed: 0 0 0 0 3 0]");
}

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
