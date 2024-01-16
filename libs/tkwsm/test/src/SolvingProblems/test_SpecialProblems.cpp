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

#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <tkrng/RNG.hpp>
#include <tkwsm/Common/GeneralUtils.hpp>
#include <utility>

#include "../TestUtils/CheckedSolution.hpp"
#include "../TestUtils/ResumedSolutionChecker.hpp"
#include "../TestUtils/TestSettings.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

/*
These are specially designed graphs where embeddings are "obvious" to a human
just from a picture, but non-obvious for an algorithm; thus they target
specific known weaknesses.

They were also designed to try to increase code coverage,
although most didn't; but keep them anyway because they're nice tests.
*/

static GraphEdgeWeights get_star_with_weights(unsigned number_of_spokes) {
  REQUIRE(number_of_spokes >= 2);
  GraphEdgeWeights result;
  for (unsigned ii = 1; ii <= number_of_spokes; ++ii) {
    result[get_edge(0, ii)] = ii;
  }
  return result;
}

// For embedding a pattern start into a target star,
// the optimal solution is unchanged as long as the target star has
// enough spokes (we just ignore the higher weight ones).
static WeightWSM get_optimal_solution_scalar_product(
    unsigned pattern_number_of_spokes) {
  REQUIRE(pattern_number_of_spokes >= 2);
  // If A, B are real sequences with A increasing, then
  //         min sum A[i].B[p(i)]
  // (where p ranges over all permutations) is solved by making
  // A,B have opposite order, i.e. B is decreasing.
  // (To maximise, we would make B increasing).
  // (Some simple algebra: consider what happens if B[i], B[j] are not in
  // the required order, and we swap them).
  // (Of course, one should read Hardy, Littlewood, Polya: "Inequalities"!)

  // Thus, we require the value 1.n + 2.(n-1) + 3.(n-2) + ... + n.1,
  // another fun exercise.
  return (pattern_number_of_spokes * (pattern_number_of_spokes + 1) *
          (pattern_number_of_spokes + 2)) /
         6;
}

// Superficially, these appear to be trivial problems; a human can easily
// write down a simple formula for the optimal solution.
// HOWEVER, unless an algorithm is
// clever enough to exclude many permutations based upon weight,
// AND uses a good heuristic to get a reasonable solution quickly,
// a naive algorithm might have to search through many permutations;
// hence these problems could quickly become very hard.
// Note that this is ENTIRELY about weight pruning;
// unweighted graph-theoretic considerations are useless.

/*
typical printouts:
for original vertex labels:
some shorter problems:

@@ Star[4] -> Star[5]: ; time 0+0; 12 iters; known opt.val. 20
@@ Star[4] -> Star[6]: ; time 0+0; 36 iters; known opt.val. 20
@@ Star[4] -> Star[7]: ; time 0+0; 110 iters; known opt.val. 20
@@ Star[4] -> Star[8]: ; time 0+0; 200 iters; known opt.val. 20
@@ Star[4] -> Star[9]: ; time 0+0; 119 iters; known opt.val. 20

some longer problems:

@@ Star[8] -> Star[9]: ; time 0+81; 16602 iters; known opt.val. 120
@@ Star[8] -> Star[10]: ; time 0+134; 40445 iters; known opt.val. 120
@@ Star[9] -> Star[9]: ; time 0+383; 51501 iters; known opt.val. 165
@@ Star[9] -> Star[10]: ; time 0+600; 117221 iters; known opt.val. 165

Thus, the algorithm is at least doing some sensible pruning.
However, with random relabelling the timings can change a lot:

@@ Star[4] -> Star[5]: ; time 0+0; 24 iters; known opt.val. 20
@@ Star[4] -> Star[6]: ; time 0+0; 136 iters; known opt.val. 20
@@ Star[4] -> Star[7]: ; time 0+0; 108 iters; known opt.val. 20
@@ Star[4] -> Star[8]: ; time 0+0; 291 iters; known opt.val. 20
@@ Star[4] -> Star[9]: ; time 0+1; 640 iters; known opt.val. 20

@@ Star[8] -> Star[9]: ; time 0+100; 20120 iters; known opt.val. 120
@@ Star[8] -> Star[10]: ; time 0+360; 102300 iters; known opt.val. 120
@@ Star[9] -> Star[9]: ; time 0+366; 51907 iters; known opt.val. 165
@@ Star[9] -> Star[10]: ; time 0+328; 68633 iters; known opt.val. 165
*/
SCENARIO("Solve WSM for star graphs") {
  const unsigned min_num_spokes = 2;
  const unsigned max_num_spokes = 6;
  std::stringstream ss;
  ss << "Star graphs: for E in [" << min_num_spokes << "," << max_num_spokes
     << "], all-against-all";
  CheckedSolution::Statistics stats(ss.str());
  CheckedSolution::ProblemInformation problem_info;
  const MainSolverParameters solver_params(10000);

  // Let's also try randomly changing the labels.
  // elem[0] will be for labels 0,1,2,... in the usual order;
  // elem[1] will have them relabelled.
  std::array<std::vector<GraphEdgeWeights>, 2> all_graph_data;
  {
    for (auto ii = min_num_spokes; ii <= max_num_spokes; ++ii) {
      all_graph_data[0].emplace_back(get_star_with_weights(ii));
    }
    all_graph_data[1].resize(all_graph_data[0].size());
    std::vector<VertexWSM> new_labels(max_num_spokes + 1);
    std::iota(new_labels.begin(), new_labels.end(), 0);
    RNG rng;

    for (unsigned ii = 0; ii < all_graph_data[0].size(); ++ii) {
      rng.do_shuffle(new_labels);
      for (const auto& original_entry : all_graph_data[0][ii]) {
        all_graph_data[1][ii][get_edge(
            new_labels[original_entry.first.first],
            new_labels[original_entry.first.second])] = original_entry.second;
      }
    }
  }
  ResumedSolutionChecker checker;

  for (unsigned flag = 0; flag <= 1; ++flag) {
    if (flag == 0) {
      TestSettings::get().os << "\n\n@@@@ original vertex labels";
    } else {
      TestSettings::get().os << "\n\n@@@@ random vertex labels";
    }
    const auto& gdata_list = all_graph_data[flag];
    for (const auto& pdata : gdata_list) {
      for (const auto& tdata : gdata_list) {
        const auto p_edges = pdata.size();
        const auto t_edges = tdata.size();

        if (p_edges <= t_edges) {
          TestSettings::get().os << "\n@@ Star[" << p_edges << "] -> Star["
                                 << t_edges << "]: ";
          problem_info.existence = CheckedSolution::ProblemInformation::
              SolutionsExistence::KNOWN_TO_BE_SOLUBLE;
          problem_info.known_optimal_solution =
              get_optimal_solution_scalar_product(p_edges);
        } else {
          problem_info.existence = CheckedSolution::ProblemInformation::
              SolutionsExistence::KNOWN_TO_BE_INSOLUBLE;
          problem_info.known_optimal_solution = {};
        }
        checker.check(
            CheckedSolution(pdata, tdata, problem_info, solver_params, stats),
            pdata, tdata, solver_params);
      }
    }
  }
  stats.finish();
}

// A graph of the form:   >-----<
// Although a human can immediately just see the answer,
// it's not so obvious how an algorithm can detect
// that this graph doesn't embed into a similar but wider one.
// The two ends could be very far apart,
// so simple distance counts will fail for wide enough arrows.
static GraphEdgeWeights get_double_arrow(unsigned width) {
  GraphEdgeWeights gdata;
  for (unsigned ii = 0; ii < width; ++ii) {
    gdata[get_edge(ii, ii + 1)] = 1 + (ii % 2);
  }
  gdata[get_edge(0, width + 1)] = 1;
  gdata[get_edge(0, width + 2)] = 2;
  gdata[get_edge(width, width + 3)] = 1;
  gdata[get_edge(width, width + 4)] = 2;
  return gdata;
}

SCENARIO("Embed double arrow graphs (stretched H graphs)") {
  std::vector<GraphEdgeWeights> gdata_list(10);
  for (unsigned ii = 0; ii < gdata_list.size(); ++ii) {
    gdata_list[ii] = get_double_arrow(ii * 10);
  }
  CheckedSolution::ProblemInformation problem_info;
  const MainSolverParameters solver_params(1000);

  CheckedSolution::Statistics stats("double arrows");
  for (const auto& p_gdata : gdata_list) {
    for (const auto& t_gdata : gdata_list) {
      problem_info.existence =
          (p_gdata.size() == t_gdata.size())
              ? CheckedSolution::ProblemInformation::SolutionsExistence::
                    KNOWN_TO_BE_SOLUBLE
              : CheckedSolution::ProblemInformation::SolutionsExistence::
                    KNOWN_TO_BE_INSOLUBLE;

      CheckedSolution(p_gdata, t_gdata, problem_info, solver_params, stats);
    }
  }
  stats.finish();
  CHECK(stats.failure_count == 0);
  CHECK(stats.timeout_count == 0);
}

namespace {
struct SpokeWithBobble {
  unsigned number_of_spokes;
  unsigned spoke_length;
  unsigned number_of_leaves_on_bobble;
  GraphEdgeWeights gdata;

  bool embeds_into_other(const SpokeWithBobble& other) const {
    return number_of_spokes <= other.number_of_spokes &&
           spoke_length == other.spoke_length &&
           number_of_leaves_on_bobble <= other.number_of_leaves_on_bobble;
  }

  unsigned num_edges() const {
    return number_of_spokes * (spoke_length + number_of_leaves_on_bobble);
  }

  void fill_gdata() {
    REQUIRE(spoke_length > 0);
    REQUIRE(number_of_spokes > 0);
    REQUIRE(number_of_leaves_on_bobble > 1);
    REQUIRE(gdata.empty());

    // Overkill, to ensure no vertex clashes...
    VertexWSM large_v = (10 + number_of_spokes) * (10 + spoke_length);

    for (unsigned ss = 0; ss < number_of_spokes; ++ss) {
      // spokes:  x--0--1--2   x--3--4--5   x--6--7--8 ... for some x >> 1.
      const unsigned first_v = ss * spoke_length;
      for (unsigned ll = 0; ll + 1 < spoke_length; ++ll) {
        gdata[get_edge(first_v + ll, first_v + ll + 1)] = 1 + ((ss + ll) % 3);
      }
      gdata[get_edge(large_v, first_v)] = 1;
    }

    // bobbles:  2--v  2--(v+1)  2--(v+2) ..., for some large v.
    for (unsigned ss = 0; ss < number_of_spokes; ++ss) {
      const VertexWSM bobble_center = (ss + 1) * spoke_length - 1;
      for (unsigned bb = 0; bb < number_of_leaves_on_bobble; ++bb) {
        ++large_v;
        gdata[get_edge(bobble_center, large_v)] = 1 + (ss % 2) + (bb % 3);
      }
    }
    // How many edges?
    REQUIRE(gdata.size() == num_edges());
    // These are, of course, TREES! So V=E+1.
    REQUIRE(get_vertices(gdata).size() == num_edges() + 1);
  }

  std::string str() const {
    std::stringstream ss;
    ss << "[" << number_of_spokes << "," << spoke_length << ","
       << number_of_leaves_on_bobble << "; " << gdata.size() << "]";
    return ss.str();
  }
};
}  // namespace

// Even though these are all trees, and it's "obvious" from drawing a picture,
// these actually can be quite testing even for V~20.
// Definitely should look at trying to improve these cases.
SCENARIO("Embed spokes_with_bobbles") {
  std::vector<SpokeWithBobble> graphs;
  SpokeWithBobble single_datum;
  for (single_datum.number_of_spokes = 1; single_datum.number_of_spokes < 4;
       ++single_datum.number_of_spokes) {
    for (single_datum.spoke_length = 2; single_datum.spoke_length < 6;
         single_datum.spoke_length += 2) {
      for (single_datum.number_of_leaves_on_bobble = 2;
           single_datum.number_of_leaves_on_bobble < 6;
           single_datum.number_of_leaves_on_bobble += 3) {
        graphs.emplace_back(single_datum);
        graphs.back().fill_gdata();
      }
    }
  }
  CheckedSolution::ProblemInformation problem_info;

  // 200 ms is fine for normal runs, but Valgrind etc. is slower.
  const MainSolverParameters solver_params(50 * 200);

  CheckedSolution::Statistics stats(
      std::string("spokes_with_bobbles : all-against-all"), graphs.size());

  unsigned counter = 0;

  const std::set<std::string> expected_problems_with_timeout{
      "[3,2,5; 21] -> [3,2,5; 21]", "[3,4,5; 27] -> [3,4,5; 27]"};
  std::set<std::string> calc_problems_with_timeout;
  std::string problem_string;

  for (const auto& p_graph : graphs) {
    for (const auto& t_graph : graphs) {
      const bool soluble = p_graph.embeds_into_other(t_graph);
      problem_info.existence =
          soluble ? CheckedSolution::ProblemInformation::SolutionsExistence::
                        KNOWN_TO_BE_SOLUBLE
                  : CheckedSolution::ProblemInformation::SolutionsExistence::
                        KNOWN_TO_BE_INSOLUBLE;

      problem_string = p_graph.str() + std::string(" -> ") + t_graph.str();

      if (soluble) {
        TestSettings::get().os << "\nN=" << counter << ": " << problem_string;
      }
      ++counter;
      if (expected_problems_with_timeout.count(problem_string) != 0) {
        TestSettings::get().os << " SKIP; takes >8 secs.\n";
        calc_problems_with_timeout.insert(problem_string);
        continue;
      }
      const auto timeout_count = stats.timeout_count;
      const CheckedSolution checked_solution(
          p_graph.gdata, t_graph.gdata, problem_info, solver_params, stats);

      if (timeout_count != stats.timeout_count) {
        calc_problems_with_timeout.insert(problem_string);
      }
      if (soluble) {
        TestSettings::get().os << "\n";
      }
    }
  }
  stats.finish();

  CHECK(stats.failure_count == 0);
  CHECK(stats.timeout_count == 0);
  CHECK(expected_problems_with_timeout == calc_problems_with_timeout);
}

namespace {
struct TreeParameters {
  unsigned number_of_vertices;
  GraphEdgeWeights gdata;

  void fill_gdata(RNG& rng) {
    REQUIRE(gdata.empty());
    REQUIRE(number_of_vertices > 2);
    gdata[get_edge(0, 1)] = 1;
    gdata[get_edge(1, 2)] = 2;
    for (std::size_t current_number_of_vertices = 3;
         current_number_of_vertices < number_of_vertices;
         ++current_number_of_vertices) {
      // Sprout a new edge from an existing vertex.
      const VertexWSM existing_v =
          rng.get_size_t(current_number_of_vertices - 1);
      gdata[get_edge(existing_v, current_number_of_vertices)] =
          1 + (current_number_of_vertices % 3);
    }
    // As always for a tree, V=E+1.
    REQUIRE(number_of_vertices == gdata.size() + 1);
  }
};
}  // namespace

// Searching for "subgraph isomorphism problem tree" shows many references;
// there definitely are clever polynomial time algorithms for the unweighted
// case. However, doesn't really help us much for the weighted case...
// ALTHOUGH, if the number of embeddings is small, we could simply
// enumerate all unweighted embeddings and choose the best one!
// Anyway, far too complicated to implement, for little benefit;
// in our applications, pattern graphs have no reason why they
// should be trees; and target graphs almost certainly will not be!
SCENARIO("Embed random trees") {
  RNG rng;
  std::vector<TreeParameters> tree_list(25);
  for (unsigned ii = 0; ii < tree_list.size(); ++ii) {
    tree_list[ii].number_of_vertices = 4 + (ii / 2);
    tree_list[ii].fill_gdata(rng);
  }
  // Impossible TV are quite rare; this is, so far,
  // the ONLY test where we found some!
  std::stringstream impossible_tv_stream;

  CheckedSolution::ProblemInformation problem_info;
  const MainSolverParameters solver_params(1000);

  CheckedSolution::Statistics stats(
      "Random trees; all-against-all", tree_list.size());

  for (unsigned ii = 0; ii < tree_list.size(); ++ii) {
    TestSettings::get().os << "\n";
    for (unsigned jj = 0; jj < tree_list.size(); ++jj) {
      if (ii == jj) {
        // The start of a run of maybe possible problems.
        TestSettings::get().os << "\n";
      }
      if (ii <= jj) {
        // It may be possible.
        TestSettings::get().os << "\nG[" << ii << "] -> G[" << jj
                               << "] (V:" << tree_list[ii].number_of_vertices
                               << "," << tree_list[jj].number_of_vertices
                               << ")";
      }
      if (ii == jj) {
        problem_info.existence = CheckedSolution::ProblemInformation::
            SolutionsExistence::KNOWN_TO_BE_SOLUBLE;
      } else {
        if (tree_list[ii].number_of_vertices <=
            tree_list[jj].number_of_vertices) {
          problem_info.existence =
              CheckedSolution::ProblemInformation::SolutionsExistence::UNKNOWN;
        } else {
          problem_info.existence = CheckedSolution::ProblemInformation::
              SolutionsExistence::KNOWN_TO_BE_INSOLUBLE;
        }
      }
      const CheckedSolution checked_solution(
          tree_list[ii].gdata, tree_list[jj].gdata, problem_info, solver_params,
          stats);

      if (!checked_solution.impossible_target_vertices.empty()) {
        impossible_tv_stream
            << "(" << ii << "," << jj
            << "):" << checked_solution.impossible_target_vertices.size()
            << " ";
      }
    }
  }
  stats.finish();

  CHECK(stats.failure_count == 0);
  CHECK(stats.timeout_count == 0);
  CHECK(impossible_tv_stream.str() == "(4,11):1 (4,16):1 ");
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
