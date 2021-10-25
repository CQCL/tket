#include <algorithm>
#include <catch2/catch.hpp>

#include "../Data/FixedCompleteSolutions.hpp"
#include "../Data/FixedSwapSequences.hpp"
#include "SwapSequenceReductionTester.hpp"

;
using std::vector;

// NOTE: running all tests in this file currently takes ~19 seconds
// on an ordinary Windows laptop.

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

// Reduce the fixed swap sequences, with edge set implicitly defined
// by the swaps themselves.
SCENARIO("Fixed swap sequences reduction") {
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

  const FixedSwapSequences fixed_sequences;
  SwapSequenceReductionTester tester;
  SwapSequenceReductionTester::Options options;
  vector<std::string> calc_messages;

  const auto add_solutions = [&tester, &options](
                                 const vector<std::string>& seq_codes,
                                 SequenceReductionStats& stats) {
    for (const auto& code_str : seq_codes) {
      const DecodedProblemData problem_data(code_str);
      const auto reduced_size =
          tester.get_checked_solution_size(problem_data, options);
      stats.add_solution(problem_data.swaps.size(), reduced_size);
    }
  };

  for (int ii = 0; ii < 2; ++ii) {
    options.optimise_initial_segment_only = (ii % 2 == 0);
    {
      SequenceReductionStats full_tokens_stats;
      add_solutions(fixed_sequences.full, full_tokens_stats);
      add_solutions(fixed_sequences.full_with_errors, full_tokens_stats);
      add_message(full_tokens_stats, "Full tokens", options, calc_messages);
    }
    {
      SequenceReductionStats partial_tokens_stats;
      add_solutions(fixed_sequences.partial, partial_tokens_stats);
      add_solutions(fixed_sequences.partial_with_errors, partial_tokens_stats);
      add_message(
          partial_tokens_stats, "Partial tokens", options, calc_messages);
    }
  }
  check_final_messages(expected_messages, calc_messages);
}

// The actual problem input data: the graph may have extra edges
// not present in the returned solution.
SCENARIO("Fixed complete problems") {
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
  vector<std::string> calc_messages;
  add_message(stats[0], "Small", options, calc_messages);
  add_message(stats[1], "Medium", options, calc_messages);
  add_message(stats[2], "Large", options, calc_messages);
  check_final_messages(expected_messages, calc_messages);
}

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
