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
#include <map>
#include <utility>

#include "../TestUtils/CheckedSolution.hpp"
#include "../TestUtils/GraphGeneration.hpp"
#include "../TestUtils/ResumedSolutionChecker.hpp"
#include "../TestUtils/TestSettings.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

// First N elements: the 64-bit uint representing the graphs with weights.
// Next element is 0, to denote the end of the sequence.
// After that, the results of embedding graph i into j, listed in order;
// 0 meaning impossible.
typedef std::vector<std::uint64_t> EncodedSolvedProblems;

// KEY: is the problem name  VALUE: the collection of solved problems.
static std::map<std::string, EncodedSolvedProblems> get_data() {
  return {
      // clang-format off
    {"Density 50%",
      { 0x1093fb7292ecde4, 0x9372a0ee562901cc, 0x196df104e143cde2,
        0x4e1bc8532fd80f73, 0xadf9bf4ee6c8c7a0, 0x43ed2d52c30b3dd0,
        0x6cf3de54709208b3, 0, 847, 0, 0, 0, 0, 0, 0, 0, 530, 0, 0, 0, 0, 0,
        0, 0, 737, 0, 0, 0, 0, 0, 0, 0, 977, 0, 0, 0, 0, 0, 0, 0, 839, 0, 0,
        0, 0, 0, 0, 0, 881, 0, 0, 0, 0, 0, 0, 0, 832,
      }},

    {"Varying density1",
      { 0x602001850028000, 0x40000000010a0022, 0x480400101001020,
        0x39404005c30000, 0x18441680401004, 0x5a05411504000868,
        0xc2998c9048805a88, 0x51e0849148801350, 0xf799e5e09ab0fd07,
        0x8eb32ea49e57883c, 0, 75, 0, 0, 0, 0, 30, 57, 0, 35, 62, 0, 57, 0,
        61, 33, 18, 33, 18, 30, 33, 0, 0, 37, 56, 16, 13, 19, 16, 16, 22, 0,
        0, 0, 264, 51, 36, 65, 51, 48, 65, 0, 0, 0, 0, 56, 0, 50, 0, 45, 48,
        0, 0, 0, 0, 0, 91, 0, 0, 85, 137, 0, 0, 0, 0, 0, 0, 344, 0, 326, 0, 0,
        0, 0, 0, 0, 0, 0, 97, 200, 0, 0, 0, 0, 0, 0, 0, 0, 0, 716, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 672,
       }},

      {"Varying density2",
        { 0x8000208204000400, 0x100000c800, 0x8200040000309,
          0x800218ad88000080, 0x10000040124028, 0xc83208010000004c,
          0x9223095222193c20, 0xab993309a0c89072, 0x81b90406380c5ad1,
          0x3fa7570d8be0e861, 0, 66, 0, 0, 48, 0, 0, 21, 24, 21, 26, 29, 34,
          26, 26, 17, 29, 14, 14, 14, 14, 0, 0, 58, 53, 0, 0, 38, 26, 26, 26,
          0, 0, 0, 211, 0, 0, 0, 165, 90, 90, 0, 0, 0, 40, 52, 0, 22, 31, 19,
          22, 0, 0, 0, 0, 0, 243, 56, 100, 64, 83, 0, 0, 0, 0, 0, 0, 393, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 415, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 444,
          259, 0, 0, 0, 0, 0, 0, 0, 0, 0, 864,
       }},

      {"Varying density3",
        { 0x1024a88200041800, 0x20045066a208406, 0x805310ab1000401c,
          0xb50a425484146021, 0xc3546a06c0508080, 0x27808029a4222ff,
          0xd1f3e263335c894e, 0x66c8093a2ae09892, 0xf40f394a7123b317,
          0x687a2fd98bcc4f60, 0, 98, 0, 0, 53, 0, 135, 63, 73, 99, 90, 0, 135,
          0, 84, 0, 0, 152, 186, 0, 128, 0, 0, 313, 99, 0, 0, 150, 144, 120,
          144, 0, 0, 0, 204, 0, 0, 0, 0, 0, 148, 0, 0, 0, 0, 218, 0, 0, 197,
          0, 211, 0, 0, 0, 0, 0, 534, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 913, 0, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 472, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 897, 0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 879,
       }},
      // clang-format on
  };
}

/*
typical printout:
@@@@ Testing 'Density 50%' : timeout=200, got 7 probs, testing all-against-all.
@@@@ Testing 'Varying density1' : timeout=200, got 10 probs, testing
all-against-all.
@@@@ Testing 'Varying density2' : timeout=200, got 10 probs, testing
all-against-all.
@@@@ Testing 'Varying density3' : timeout=200, got 10 probs, testing
all-against-all.
 @@@ fixed small graphs fin. Total time 0+236
*/

SCENARIO("embedding all-against-all") {
  const auto solved_problems_map = get_data();
  EncodedSolvedProblems calc_problems;

  // These problems are small and easy.
  // <300ms TOTAL for the whole set, no timeouts anywhere near being hit.
  const unsigned timeout_ms = 1000;
  CheckedSolution::Statistics statistics("fixed small graphs");
  const MainSolverParameters solver_params(timeout_ms);
  const CheckedSolution::ProblemInformation info;
  ResumedSolutionChecker resumption_checker;

  // Holds the decoded graph data.
  std::vector<GraphEdgeWeights> gdata;
  const auto& os = TestSettings::get().os;

  // Go through all problem sets.
  for (const auto& entry : solved_problems_map) {
    os << "\n@@@@ Testing '" << entry.first << "' : timeout=" << timeout_ms;
    calc_problems.clear();

    // Decode the fixed problems and copy them into gdata.
    gdata.clear();
    for (auto code : entry.second) {
      calc_problems.push_back(code);
      if (calc_problems.back() == 0) {
        break;
      }
      const GraphGeneration::LimitedSizeGraphGeneral graph(code);
      gdata.push_back(graph.data);
    }
    REQUIRE(!gdata.empty());
    if (calc_problems.back() != 0) {
      calc_problems.push_back(0);
    }
    os << ", got " << gdata.size() << " probs, testing all-against-all.";

    // Pick out a specific test easily.
    // unsigned counter=0;

    // Now test all-against-all, and write the results into "calc_problems".
    for (unsigned ii = 0; ii < gdata.size(); ++ii) {
      REQUIRE(!gdata[ii].empty());

      // If an edge exists, it must have nonzero weight.
      for (const auto& entry : gdata[ii]) {
        REQUIRE(entry.second > 0);
      }

      for (unsigned jj = 0; jj < gdata.size(); ++jj) {
        //++counter;
        // os << "\n@@@@@@ COUNTER=" << counter << "; g[" << ii << "] -> g[" <<
        // jj << "]:\nPE: " << str(gdata[ii])
        //  << "\nTE: " << str(gdata[jj]) << "\n";

        const CheckedSolution checked_solution(
            gdata[ii], gdata[jj], info, solver_params, statistics);

        resumption_checker.check(
            checked_solution, gdata[ii], gdata[jj], solver_params);

        // Should be no timeouts!
        CHECK(checked_solution.finished);
        if (!checked_solution.assignments.empty()) {
          // All edge weights are positive.
          CHECK(checked_solution.scalar_product > 0);
          calc_problems.push_back(checked_solution.scalar_product);
          if (ii == jj) {
            // A self-embedding must exist. We can bound the solution.
            WeightWSM weight = 0;
            for (const auto& entry : gdata[ii]) {
              weight += entry.second * entry.second;
            }
            CHECK(checked_solution.scalar_product <= weight);
          }
          continue;
        }
        // no solution.
        calc_problems.push_back(0);
        // Self-embedding is always possible.
        CHECK(ii != jj);
      }
    }
    CHECK(entry.second == calc_problems);
  }
  statistics.finish();
  CHECK(statistics.success_count == 349);
  CHECK(statistics.failure_count == 0);
  CHECK(statistics.timeout_count == 0);
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
