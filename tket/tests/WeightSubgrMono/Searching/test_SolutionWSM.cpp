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
#include <cmath>
#include <map>
#include <utility>

#include "WeightSubgrMono/Common/GeneralUtils.hpp"
#include "WeightSubgrMono/Searching/SolutionWSM.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

SCENARIO("SolutionWSM : assignments") {
  SolutionWSM solution;
  GraphEdgeWeights edges_and_weights;
  edges_and_weights[get_edge(0, 1)] = 0;

  solution.complete = true;

  const auto check_expected_errors =
      [&solution, &edges_and_weights](
          const std::vector<std::string>& expected_substrings,
          const std::vector<std::string>& expected_absent_substrings = {}) {
        const auto message =
            solution.get_errors(edges_and_weights, edges_and_weights);
        for (const auto& expected : expected_substrings) {
          CHECK_THAT(message, Catch::Contains(expected));
        }
        for (const auto& absent : expected_absent_substrings) {
          CHECK_THAT(message, !Catch::Contains(absent));
        }
      };

  // missing p-vertices 0,1.
  check_expected_errors(
      {"P-edge (0,1) has unknown vertices",
       "number of used p vertices mismatch"});

  solution.assignments.emplace_back(0, 0);
  // missing p-vertex 1.
  check_expected_errors(
      {"P-edge (0,1) has unknown vertices",
       "number of used p vertices mismatch"});

  solution.assignments.emplace_back(1, 1);
  // finally correct! 0->0, 1->1.
  CHECK("" == solution.get_errors(edges_and_weights, edges_and_weights));

  solution.assignments.emplace_back(1, 1);
  // 1->1 is repeated. But no UNKNOWN vertices!
  check_expected_errors(
      {"Repeated assignments", "Duplicate value", "Sizes mismatch",
       "number of used p vertices mismatch"},
      {"has unknown vertices "});

  solution.assignments.pop_back();
  solution.assignments.emplace_back(1, 0);
  // 1->0 contradicts 1->1.
  check_expected_errors(
      {"Repeated assignments", "Duplicate value", "Sizes mismatch",
       "P vertices", "both map to", "number of used p vertices mismatch"});

  solution.assignments.pop_back();
  solution.assignments.emplace_back(2, 0);
  // 1->0 clashes with 2->0 (2 also being "unknown", but it gives up before
  // this).
  check_expected_errors(
      {"Duplicate value", "Sizes mismatch",
       "number of used p vertices mismatch"});

  solution.assignments.pop_back();
  solution.assignments.emplace_back(2, 2);
  // 2->2 is unknown. Sees too many p-vertices {0,1,2}.
  check_expected_errors(
      {"number of used p vertices mismatch"}, {"Repeated", "Duplicate"});
}

SCENARIO("SolutionWSM : int overflow") {
  GraphEdgeWeights p_edges_and_weights;
  const WeightWSM max_weight = std::numeric_limits<WeightWSM>::max();
  const auto edge = get_edge(0, 1);
  p_edges_and_weights[edge] = max_weight / WeightWSM(100);

  SolutionWSM solution;
  solution.assignments.emplace_back(0, 0);
  solution.assignments.emplace_back(1, 1);
  solution.complete = true;

  // (M/100)^2 definitely should overflow.
  CHECK(
      "\nOverflow: w(p-edge) * w(t-edge): "
      "184467440737095516*184467440737095516" ==
      solution.get_errors(p_edges_and_weights, p_edges_and_weights));

  GraphEdgeWeights t_edges_and_weights;
  for (WeightWSM x = 90; x <= 100; ++x) {
    // "M/100" is actually just under the real value of M/100,
    // so multiplying by 100, it will not overflow.
    t_edges_and_weights[edge] = x;
    const auto errors =
        solution.get_errors(p_edges_and_weights, t_edges_and_weights);
    CHECK_THAT(errors, !Catch::Contains("verflow"));
    CHECK_THAT(errors, Catch::Contains("Recalc/orig weights mismatch"));
  }
  t_edges_and_weights[edge] = 101;
  CHECK_THAT(
      solution.get_errors(p_edges_and_weights, t_edges_and_weights),
      Catch::Contains("Overflow: w(p-edge) * w(t-edge):"));

  // Now, one line into another, overflows only at the end.
  p_edges_and_weights.clear();
  t_edges_and_weights.clear();
  const unsigned vertices = 50;
  const WeightWSM big_weight = max_weight / vertices;

  for (unsigned ii = 0; ii <= vertices + 10; ++ii) {
    t_edges_and_weights[get_edge(ii, ii + 1)] = 1;
  }
  for (unsigned ii = 0; ii <= vertices + 10; ++ii) {
    if (ii >= 1) {
      solution.assignments.emplace_back(ii + 1, ii + 1);
    }
    // Assignments 0->0, 1->1, ..., (i+1)->(i+1),
    // and p-edges (0,1), (1,2), ..., (i,i+1), that is, i+1 of them.
    p_edges_and_weights[get_edge(ii, ii + 1)] = big_weight;
    if (ii < 3) {
      continue;
    }
    const auto message =
        solution.get_errors(p_edges_and_weights, t_edges_and_weights);
    if (ii + 1 <= vertices) {
      CHECK_THAT(message, Catch::Contains("Recalc/orig weights mismatch"));
      CHECK_THAT(message, !Catch::Contains("verflow"));
    } else {
      // Definitely should overflow now... (M/v)*(i+1) for i+1>v.
      CHECK_THAT(message, !Catch::Contains("mismatch"));
      CHECK_THAT(
          message, Catch::Contains("Overflow calculating total p-weight:"));
      CHECK_THAT(
          message, Catch::Contains("Overflow calculating total weight:"));
    }
  }
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
