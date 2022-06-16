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
#include <functional>
#include <set>
#include <sstream>

#include "TestUtils/DebugFunctions.hpp"
#include "TokenSwapping/SwapListOptimiser.hpp"
#include "Utils/RNG.hpp"

using std::vector;

namespace tket {
namespace tsa_internal {
namespace tests {

namespace {

// Only checks that swaps are correct, doesn't measure how good they are
class SwapCorrectnessTester {
 public:
  // Perform the raw swaps for comparison.
  void reset(const vector<Swap>& raw_swaps) {
    m_final_tracker.reset();
    for (const auto& swap : raw_swaps) {
      (void)m_final_tracker.do_vertex_swap(swap);
    }
    m_number_of_raw_swaps = raw_swaps.size();
  }

  void require_equal_permutations(const SwapList& swap_list) const {
    m_tracker_to_change.reset();
    size_t num_swaps = 0;
    for (auto id = swap_list.front_id(); id; id = swap_list.next(id.value())) {
      m_tracker_to_change.do_vertex_swap(swap_list.at(id.value()));
      ++num_swaps;
    }
    REQUIRE(num_swaps == swap_list.size());
    REQUIRE(m_tracker_to_change.equal_vertex_permutation_from_swaps(
        m_final_tracker));
    REQUIRE(m_number_of_raw_swaps >= num_swaps);
  }

 private:
  size_t m_number_of_raw_swaps = 0;
  DynamicTokenTracker m_final_tracker;
  mutable DynamicTokenTracker m_tracker_to_change;
};

// As well as correctness, also checks that optimisation passes
// do actually perform quite well.
class SwapTester {
 public:
  SwapTester() {
    m_optimisation_functions.reserve(5);

    m_optimisation_functions.emplace_back([](const vector<Swap>& raw_swaps,
                                             SwapList& list,
                                             SwapListOptimiser& optimiser) {
      for (const Swap& swap : raw_swaps) {
        optimiser.push_back(list, swap);
      }
    });
    m_optimisation_functions.emplace_back(
        [](const vector<Swap>&, SwapList& list, SwapListOptimiser& optimiser) {
          optimiser.optimise_pass_with_zero_travel(list);
        });
    m_optimisation_functions.emplace_back(
        [](const vector<Swap>&, SwapList& list, SwapListOptimiser& optimiser) {
          optimiser.optimise_pass_with_frontward_travel(list);
        });
    m_optimisation_functions.emplace_back(
        [](const vector<Swap>&, SwapList& list, SwapListOptimiser& optimiser) {
          optimiser.optimise_pass_with_token_tracking(list);
        });
    m_optimisation_functions.emplace_back(
        [](const vector<Swap>&, SwapList& list, SwapListOptimiser& optimiser) {
          optimiser.full_optimise(list);
        });
    reset_counters();
  }

  void reset_counters() {
    // Also includes the number of raw swaps,
    // and the number of tests.
    m_counts.resize(m_optimisation_functions.size() + 2);
    std::fill(m_counts.begin(), m_counts.end(), 0);
  }

  void test(const vector<Swap>& raw_swaps) {
    ++m_counts[0];
    m_counts[1] += raw_swaps.size();
    m_correctness_tester.reset(raw_swaps);

    for (size_t ii = 0; ii < m_optimisation_functions.size(); ++ii) {
      m_swap_list.clear();
      if (ii != 0) {
        for (const auto& swap : raw_swaps) {
          m_swap_list.push_back(swap);
        }
      }
      m_optimisation_functions[ii](raw_swaps, m_swap_list, m_optimiser);
      m_correctness_tester.require_equal_permutations(m_swap_list);
      m_counts[ii + 2] += m_swap_list.size();
    }
  }

  std::string get_final_result() const {
    std::stringstream ss;
    ss << "[ " << m_counts[0] << " tests; swap counts:";
    for (size_t ii = 1; ii < m_counts.size(); ++ii) {
      ss << " " << m_counts[ii] << " ";
    }
    ss << "]";
    return ss.str();
  }

 private:
  vector<
      std::function<void(const vector<Swap>&, SwapList&, SwapListOptimiser&)>>
      m_optimisation_functions;

  vector<size_t> m_counts;

  SwapList m_swap_list;
  SwapListOptimiser m_optimiser;
  SwapCorrectnessTester m_correctness_tester;
  size_t number_of_tests;
};
}  // namespace

SCENARIO("Random swaps are optimised") {
  RNG rng;
  SwapTester tester;
  vector<Swap> raw_swaps;
  const vector<size_t> num_vertices{5, 10, 20};

  // We will multiply the number of possible distinct swaps
  // by these numbers, then divide by 100, to determine how many swaps
  // to generate for the test.
  const vector<size_t> percentages{50, 100, 200, 500};

  // Not necessarily contiguous.
  std::set<size_t> vertices_set;

  for (size_t number_of_vertices : num_vertices) {
    const size_t possible_swaps =
        (number_of_vertices * (number_of_vertices - 1)) / 2;
    for (auto percent : percentages) {
      const size_t num_swaps = (possible_swaps * percent) / 100;
      vertices_set.clear();

      for (size_t ii = 0; ii < number_of_vertices; ++ii) {
        vertices_set.insert(ii);
      }
      const vector<size_t> vertices(vertices_set.cbegin(), vertices_set.cend());
      for (int test_counter = 0; test_counter < 1; ++test_counter) {
        INFO(
            "test_counter=" << test_counter << ", number_of_vertices="
                            << number_of_vertices << ", percent=" << percent);

        for (size_t jj = 0; jj < num_swaps; ++jj) {
          const auto v1 = rng.get_element(vertices);
          auto v2 = v1;
          while (v1 == v2) {
            v2 = rng.get_element(vertices);
          }
          raw_swaps.emplace_back(get_swap(v1, v2));
        }
        tester.test(raw_swaps);
      }
    }
  }
  CHECK(
      tester.get_final_result() ==
      "[ 12 tests; swap counts: 5636  5256  4976  4976  264  268 ]");
}

namespace {
// The above test just generates completely random swap sequences
// on N vertices. For a more realistic sequence, we try choosing them
// from a smaller list of possible swaps
// (thus, representing swaps on an incomplete graph).
// This might be more realistic.
struct EdgesGenerator {
  std::set<Swap> swaps_set;
  size_t approx_num_vertices = 5;
  size_t approx_num_edges = 10;
  size_t percentage_to_add_new_vertex = 50;

  vector<Swap> get_swaps(RNG& rng, size_t& actual_num_vertices) {
    actual_num_vertices = 2;
    swaps_set.clear();
    swaps_set.insert(get_swap(0, 1));

    for (size_t counter = 10 * approx_num_edges; counter > 0; --counter) {
      if (actual_num_vertices >= approx_num_vertices ||
          swaps_set.size() >= approx_num_edges) {
        break;
      }
      bool add_new_vertex = rng.check_percentage(percentage_to_add_new_vertex);
      if (!add_new_vertex) {
        const auto current_edges = swaps_set.size();
        for (int edge_attempt = 10; edge_attempt > 0; --edge_attempt) {
          const auto v1 = rng.get_size_t(actual_num_vertices - 1);
          const auto v2 = rng.get_size_t(actual_num_vertices - 1);
          if (v1 != v2) {
            swaps_set.insert(get_swap(v1, v2));
            if (current_edges != swaps_set.size()) {
              break;
            }
          }
        }
        if (current_edges != swaps_set.size()) {
          continue;
        }
        add_new_vertex = true;
      }
      if (add_new_vertex) {
        swaps_set.insert(get_swap(
            rng.get_size_t(actual_num_vertices - 1), actual_num_vertices));
        ++actual_num_vertices;
        continue;
      }
    }
    vector<Swap> result{swaps_set.cbegin(), swaps_set.cend()};
    return result;
  }
};

struct ManyTestsRunner {
  SwapTester tester;

  EdgesGenerator swaps_generator;
  vector<Swap> possible_swaps;
  size_t actual_num_vertices;
  vector<Swap> raw_swaps;

  void run(
      RNG& rng, const vector<size_t>& approx_num_vertices,
      const vector<size_t>& approx_num_edges_percentages,
      const vector<size_t>& swap_length_percentages,
      size_t num_tests_per_parameter_list) {
    for (auto approx_nv : approx_num_vertices) {
      swaps_generator.approx_num_vertices = approx_nv;
      for (auto approx_nep : approx_num_edges_percentages) {
        swaps_generator.approx_num_edges =
            approx_nv / 2 + (approx_nv * (approx_nv - 1) * approx_nep) / 200;
        for (size_t num_graphs = 0; num_graphs < 1; ++num_graphs) {
          possible_swaps = swaps_generator.get_swaps(rng, actual_num_vertices);
          for (auto slp : swap_length_percentages) {
            const size_t swap_list_length =
                1 + (possible_swaps.size() * slp) / 100;
            for (size_t test_counter = 0;
                 test_counter < num_tests_per_parameter_list; ++test_counter) {
              raw_swaps.clear();
              for (size_t nn = 0; nn < swap_list_length; ++nn) {
                raw_swaps.push_back(rng.get_element(possible_swaps));
              }
              tester.test(raw_swaps);
            }
          }
        }
      }
    }
  }
};
}  // namespace

SCENARIO("More realistic swap sequences") {
  RNG rng;
  const size_t num_tests_per_parameter_list = 10;

  // How many edges should we aim for, as a rough percentage of
  // the total number n(n-1)/2 of possibilities?
  const vector<size_t> approx_num_edges_percentages{5, 10, 20, 30, 40, 80};

  // How long should the swap length be, as a percentage of the
  // total possible number of swaps?
  const vector<size_t> swap_length_percentages{50, 100, 200};

  {
    const vector<size_t> approx_num_vertices{5, 8};
    ManyTestsRunner runner;
    runner.run(
        rng, approx_num_vertices, approx_num_edges_percentages,
        swap_length_percentages, num_tests_per_parameter_list);
    CHECK(
        runner.tester.get_final_result() ==
        "[ 360 tests; swap counts: 3160  2380  2104  2104  1396  1406 ]");
  }
  {
    const vector<size_t> approx_num_vertices{10, 12, 14};
    ManyTestsRunner runner;
    runner.run(
        rng, approx_num_vertices, approx_num_edges_percentages,
        swap_length_percentages, num_tests_per_parameter_list);
    CHECK(
        runner.tester.get_final_result() ==
        "[ 540 tests; swap counts: 10370  9048  7580  7580  5180  5216 ]");
  }
  {
    const vector<size_t> approx_num_vertices{30, 35, 40};
    ManyTestsRunner runner;
    runner.run(
        rng, approx_num_vertices, approx_num_edges_percentages,
        swap_length_percentages, num_tests_per_parameter_list);
    CHECK(
        runner.tester.get_final_result() ==
        "[ 540 tests; swap counts: 38900  37626  30944  30944  24714  "
        "24720 ]");
  }
}

// If we perform a sequence of swaps, then again in reverse order,
// (and thus, make a palindrome), it ALWAYS equals the identity permutation.
// (Of course, odd-length palindromes like "(0,1)" do NOT give the identity!)
// It seems "obvious" that zero-travel and frontwards-travel passes
// should optimise (even-length) palindromes to zero; but is it actually true?!
// Token-tracking passes definitely do NOT, but counterexamples are rare.
// (Even though token-tracking IRREDUCIBILITY can be shown to be
// STRICTLY STRONGER than zero-travel or frontwards-travel IRREDUCIBILITY!)
SCENARIO("Trivial swap list reversed order optimisation; pass comparisons") {
  vector<Swap> possible_swaps;
  const unsigned num_vertices = 4;

  for (unsigned ii = 0; ii < num_vertices; ++ii) {
    for (unsigned jj = ii + 1; jj < num_vertices; ++jj) {
      possible_swaps.push_back(get_swap(ii, jj));
    }
  }
  vector<Swap> raw_swaps;
  SwapList swaps;
  SwapListOptimiser optimiser;

  const auto push_back_swaps = [&raw_swaps, &swaps]() {
    swaps.fast_clear();
    for (auto swap : raw_swaps) {
      swaps.push_back(swap);
    }
  };

  const auto concatenate_reversed_swaps = [&raw_swaps, &swaps,
                                           &push_back_swaps]() {
    push_back_swaps();
    for (auto citer = raw_swaps.crbegin(); citer != raw_swaps.crend();
         ++citer) {
      swaps.push_back(*citer);
    }
  };

  size_t simple_travel_equals_token_tracking_count = 0;
  size_t simple_travel_beats_token_tracking_count = 0;
  size_t simple_travel_beaten_by_token_tracking_count = 0;
  size_t full_optimise_fully_reduces_palindrome = 0;
  size_t full_optimise_does_not_destroy_palindrome = 0;
  size_t token_tracking_pass_fully_reduces_palindrome = 0;
  size_t token_tracking_pass_does_not_destroy_palindrome = 0;

  RNG rng;

  for (int test_counter = 0; test_counter < 1000; ++test_counter) {
    if (raw_swaps.size() > 20) {
      raw_swaps.clear();
    }
    raw_swaps.push_back(rng.get_element(possible_swaps));

    concatenate_reversed_swaps();
    optimiser.optimise_pass_with_zero_travel(swaps);
    CHECK(swaps.size() == 0);

    concatenate_reversed_swaps();
    optimiser.optimise_pass_with_frontward_travel(swaps);
    CHECK(swaps.size() == 0);

    concatenate_reversed_swaps();
    optimiser.optimise_pass_with_token_tracking(swaps);
    if (swaps.size() == 0) {
      ++token_tracking_pass_fully_reduces_palindrome;
    } else {
      ++token_tracking_pass_does_not_destroy_palindrome;
    }

    concatenate_reversed_swaps();
    optimiser.full_optimise(swaps);
    if (swaps.size() == 0) {
      ++full_optimise_fully_reduces_palindrome;
    } else {
      ++full_optimise_does_not_destroy_palindrome;
    }

    push_back_swaps();
    optimiser.optimise_pass_with_zero_travel(swaps);
    const auto zero_travel_reduced_size = swaps.size();

    push_back_swaps();
    optimiser.optimise_pass_with_frontward_travel(swaps);
    const auto frontward_travel_reduced_size = swaps.size();
    CHECK(zero_travel_reduced_size == frontward_travel_reduced_size);

    push_back_swaps();
    optimiser.optimise_pass_with_token_tracking(swaps);

    const auto token_tracking_reduced_size = swaps.size();
    if (token_tracking_reduced_size == zero_travel_reduced_size) {
      ++simple_travel_equals_token_tracking_count;
    } else {
      if (token_tracking_reduced_size < zero_travel_reduced_size) {
        ++simple_travel_beaten_by_token_tracking_count;
      } else {
        ++simple_travel_beats_token_tracking_count;
      }
    }
  }
  CHECK(simple_travel_equals_token_tracking_count == 299);
  CHECK(simple_travel_beaten_by_token_tracking_count == 697);
  CHECK(simple_travel_beats_token_tracking_count == 4);
  CHECK(full_optimise_fully_reduces_palindrome == 1000);
  CHECK(full_optimise_does_not_destroy_palindrome == 0);
  CHECK(token_tracking_pass_fully_reduces_palindrome == 976);
  CHECK(token_tracking_pass_does_not_destroy_palindrome == 24);
}

SCENARIO("specific swap list optimisation counterexamples") {
  SwapList swaps;
  SwapListOptimiser optimiser;
  // Illustrates that general-travel irreducible does NOT imply token-tracking
  // irreducible. (Of course, we haven't IMPLEMENTED general-travel reduction,
  // but we can PROVE that general-travel irreducibility is equivalent to
  // zero-travel and frontwards-travel irreducibility).
  swaps.push_back(get_swap(0, 1));
  swaps.push_back(get_swap(0, 2));
  swaps.push_back(get_swap(0, 1));
  swaps.push_back(get_swap(0, 2));
  optimiser.optimise_pass_with_zero_travel(swaps);
  CHECK(swaps.size() == 4);
  optimiser.optimise_pass_with_frontward_travel(swaps);
  CHECK(swaps.size() == 4);
  optimiser.optimise_pass_with_token_tracking(swaps);
  CHECK(str(swaps) == " (0,2)  (0,1) ");

  // Are palindromes S + Reverse(S) ALWAYS optimised to an empty list by zero
  // travel or frontwards travel passes? Seems so, but how to prove it? (We know
  // that for IRREDUCIBILITY, zero-travel, frontwards-travel, general-travel
  // give equivalent concepts, and token-tracking gives a strictly stronger
  // pass, i.e. token-tracking irreducible => zero-travel irreducible, etc. but
  // NOT conversely. But we have no such results for sequence reduction, and
  // this counterexample illustrates that).
  const vector<Swap> swap_sequence_palindrome{
      {1, 2}, {1, 3}, {0, 2}, {1, 3}, {1, 3}, {2, 3}, {0, 1}, {1, 2},
      {0, 1}, {0, 2}, {1, 2}, {0, 3}, {0, 3}, {1, 2}, {0, 2}, {0, 1},
      {1, 2}, {0, 1}, {2, 3}, {1, 3}, {1, 3}, {0, 2}, {1, 3}, {1, 2}};
  REQUIRE(swap_sequence_palindrome.size() % 2 == 0);
  for (unsigned ii = 0; ii < swap_sequence_palindrome.size(); ++ii) {
    REQUIRE(
        swap_sequence_palindrome[ii] ==
        swap_sequence_palindrome[swap_sequence_palindrome.size() - 1 - ii]);
  }

  const auto push_back_swaps = [&swaps, &swap_sequence_palindrome]() {
    swaps.fast_clear();
    for (auto swap : swap_sequence_palindrome) {
      swaps.push_back(swap);
    }
  };

  push_back_swaps();
  optimiser.optimise_pass_with_frontward_travel(swaps);
  CHECK(swaps.size() == 0);

  push_back_swaps();
  optimiser.optimise_pass_with_zero_travel(swaps);
  CHECK(swaps.size() == 0);

  push_back_swaps();
  optimiser.optimise_pass_with_token_tracking(swaps);
  CHECK(str(swaps) == " (0,3)  (0,1)  (2,3)  (0,2)  (1,3)  (1,2) ");
}

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
