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

#include <algorithm>
#include <catch2/catch_test_macros.hpp>

#include "../Data/FixedCompleteSolutions.hpp"
#include "../Data/FixedSwapSequences.hpp"
#include "SwapSequenceReductionTester.hpp"

using std::vector;

namespace tket {
namespace tsa_internal {
namespace tests {

static void add_message(
    const SequenceReductionStats& stats, const std::string& extra_message,
    const SwapSequenceReductionTester::Options& options,
    vector<std::string>& calc_messages) {
  std::stringstream ss;
  ss << "[n=" << calc_messages.size() << ", " << extra_message
     << ": init segm optim? " << std::boolalpha
     << options.optimise_initial_segment_only << "]\n"
     << stats.str();
  calc_messages.push_back(ss.str());
}

static void check_final_messages(
    vector<std::string>& expected_messages,
    const vector<std::string>& calc_messages) {
  CHECK(expected_messages.size() == calc_messages.size());
  expected_messages.resize(calc_messages.size());
  for (unsigned ii = 0; ii < calc_messages.size(); ++ii) {
    CHECK(expected_messages[ii] == calc_messages[ii]);
  }
}

static void add_solutions(
    SwapSequenceReductionTester& tester,
    SwapSequenceReductionTester::Options& options, const unsigned skip_number,
    const vector<std::string>& seq_codes, SequenceReductionStats& stats) {
  for (unsigned ii = 0; ii < seq_codes.size(); ++ii) {
    if (ii % skip_number != 0) {
      continue;
    }
    const auto& code_str = seq_codes[ii];
    const DecodedProblemData problem_data(code_str);
    const auto reduced_size =
        tester.get_checked_solution_size(problem_data, options);
    stats.add_solution(problem_data.swaps.size(), reduced_size);
  }
}

static void run_reduction(
    SwapSequenceReductionTester& tester,
    SwapSequenceReductionTester::Options& options, const unsigned skip_number,
    const FixedSwapSequences& sequences, vector<std::string>& calc_messages) {
  for (int ii = 0; ii < 2; ++ii) {
    options.optimise_initial_segment_only = (ii % 2 == 0);
    {
      SequenceReductionStats full_tokens_stats;
      add_solutions(
          tester, options, skip_number, sequences.full, full_tokens_stats);
      add_solutions(
          tester, options, skip_number, sequences.full_with_errors,
          full_tokens_stats);
      add_message(full_tokens_stats, "Full tokens", options, calc_messages);
    }
    {
      SequenceReductionStats partial_tokens_stats;
      add_solutions(
          tester, options, skip_number, sequences.partial,
          partial_tokens_stats);
      add_solutions(
          tester, options, skip_number, sequences.partial_with_errors,
          partial_tokens_stats);
      add_message(
          partial_tokens_stats, "Partial tokens", options, calc_messages);
    }
  }
}

// Reduce the fixed swap sequences, with edge set implicitly defined
// by the swaps themselves.
// This long test take ~5 seconds on a 2021 Windows laptop.
SCENARIO("Fixed swap sequences reduction - long test", "[.long]") {
  vector<std::string> expected_messages{
      "[n=0, Full tokens: init segm optim? true]\n"
      "[478 equal probs (17115); 2 reduced probs (25 vs 29)]\n"
      "[Overall reduction 17140 vs 17144: 0%]",

      "[n=1, Partial tokens: init segm optim? true]\n"
      "[880 equal probs (25432); 16 reduced probs (385 vs 407)]\n"
      "[Overall reduction 25817 vs 25839: 0%]",

      "[n=2, Full tokens: init segm optim? false]\n"
      "[423 equal probs (14323); 57 reduced probs (2693 vs 2821)]\n"
      "[Overall reduction 17016 vs 17144: 0%]",

      "[n=3, Partial tokens: init segm optim? false]\n"
      "[658 equal probs (12376); 238 reduced probs (12962 vs 13463)]\n"
      "[Overall reduction 25338 vs 25839: 1%]"};
  const unsigned skip_number = 1;
  const FixedSwapSequences fixed_sequences;
  SwapSequenceReductionTester tester;
  SwapSequenceReductionTester::Options options;
  vector<std::string> calc_messages;

  run_reduction(tester, options, skip_number, fixed_sequences, calc_messages);
  check_final_messages(expected_messages, calc_messages);
}

// Reduce the fixed swap sequences, with edge set implicitly defined
// by the swaps themselves.
// This short test take ~0.4 seconds on a 2021 Windows laptop.
SCENARIO("Fixed swap sequences reduction") {
  vector<std::string> expected_messages{
      "[n=0, Full tokens: init segm optim? true]\n"
      "[25 equal probs (846); 0 reduced probs (0 vs 0)]\n"
      "[Overall reduction 846 vs 846: 0%]",

      "[n=1, Partial tokens: init segm optim? true]\n"
      "[46 equal probs (1348); 0 reduced probs (0 vs 0)]\n"
      "[Overall reduction 1348 vs 1348: 0%]",

      "[n=2, Full tokens: init segm optim? false]\n"
      "[24 equal probs (822); 1 reduced probs (22 vs 24)]\n"
      "[Overall reduction 844 vs 846: 0%]",

      "[n=3, Partial tokens: init segm optim? false]\n"
      "[34 equal probs (461); 12 reduced probs (844 vs 887)]\n"
      "[Overall reduction 1305 vs 1348: 3%]"};
  const unsigned skip_number = 20;
  const FixedSwapSequences fixed_sequences;
  SwapSequenceReductionTester tester;
  SwapSequenceReductionTester::Options options;
  vector<std::string> calc_messages;

  run_reduction(tester, options, skip_number, fixed_sequences, calc_messages);
  check_final_messages(expected_messages, calc_messages);
}

static void run_complete_problems(
    const unsigned skip_number, vector<std::string>& calc_messages) {
  SwapSequenceReductionTester::Options options;
  options.optimise_initial_segment_only = false;

  // Separate problems into small, medium, large.
  vector<SequenceReductionStats> stats(3);

  const FixedCompleteSolutions complete_solutions;
  SwapSequenceReductionTester tester;

  for (const auto& problem_entry : complete_solutions.solutions) {
    // First element encodes the edges.
    const DecodedArchitectureData arch_data(problem_entry.second[0]);
    for (unsigned ii = 1; ii < problem_entry.second.size(); ++ii) {
      if (ii % skip_number != 0) {
        continue;
      }
      const auto& problem_str = problem_entry.second[ii];
      const DecodedProblemData problem_data(
          problem_str, DecodedProblemData::RequireContiguousVertices::NO);

      // Small
      unsigned stats_index = 0;
      if (problem_str.size() > 25) {
        // Medium
        stats_index = 1;
      }
      if (problem_str.size() > 60) {
        // Large
        stats_index = 2;
      }
      const auto reduced_size =
          tester.get_checked_solution_size(problem_data, arch_data, options);
      stats[stats_index].add_solution(problem_data.swaps.size(), reduced_size);
    }
  }
  add_message(stats[0], "Small", options, calc_messages);
  add_message(stats[1], "Medium", options, calc_messages);
  add_message(stats[2], "Large", options, calc_messages);
}

// The actual problem input data: the graph may have extra edges
// not present in the returned solution.
// The long tests take ~10 seconds on a 2021 Windows laptop.
SCENARIO("Fixed complete problems - long test", "[.long]") {
  vector<std::string> expected_messages{
      "[n=0, Small: init segm optim? false]\n"
      "[249 equal probs (1353); 29 reduced probs (163 vs 204)]\n"
      "[Overall reduction 1516 vs 1557: 2%]",

      "[n=1, Medium: init segm optim? false]\n"
      "[167 equal probs (2650); 60 reduced probs (1107 vs 1234)]\n"
      "[Overall reduction 3757 vs 3884: 3%]",

      "[n=2, Large: init segm optim? false]\n"
      "[164 equal probs (12771); 408 reduced probs (43946 vs 45894)]\n"
      "[Overall reduction 56717 vs 58665: 3%]"};
  const unsigned skip_number = 1;
  vector<std::string> calc_messages;
  run_complete_problems(skip_number, calc_messages);
  check_final_messages(expected_messages, calc_messages);
}

// The actual problem input data: the graph may have extra edges
// not present in the returned solution.
// The shorter tests take ~0.4 seconds.
SCENARIO("Fixed complete problems") {
  vector<std::string> expected_messages{
      "[n=0, Small: init segm optim? false]\n"
      "[8 equal probs (48); 1 reduced probs (9 vs 10)]\n"
      "[Overall reduction 57 vs 58: 1%]",

      "[n=1, Medium: init segm optim? false]\n"
      "[8 equal probs (138); 1 reduced probs (23 vs 24)]\n"
      "[Overall reduction 161 vs 162: 0%]",

      "[n=2, Large: init segm optim? false]\n"
      "[10 equal probs (928); 16 reduced probs (1657 vs 1743)]\n"
      "[Overall reduction 2585 vs 2671: 3%]"};
  const unsigned skip_number = 20;
  vector<std::string> calc_messages;
  run_complete_problems(skip_number, calc_messages);
  check_final_messages(expected_messages, calc_messages);
}

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
