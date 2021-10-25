#include <catch2/catch.hpp>

#include "Data/FixedCompleteSolutions.hpp"
#include "Data/FixedSwapSequences.hpp"
#include "TestUtils/BestTsaTester.hpp"

// NOTE: currently, the tests in this file (solving ~2300 complete problems
// with the BestTSA, which includes full table lookup)
// take ~25 seconds on an ordinary Windows laptop.
//
/// TODO: The swap table optimiser currently tries to optimise many segments;
/// certainly it could be cut down, experimentation is needed
/// to find how much to cut it down, without degrading solution
/// quality too much.
//

;
using std::vector;

namespace tket {
namespace tsa_internal {
namespace tests {

namespace {
struct FixedSeqsStats {
  size_t equivalent_solns = 0;
  size_t equivalent_solns_swaps = 0;
  size_t better_solns = 0;
  size_t better_solns_swaps = 0;
  size_t better_solns_known_swaps = 0;
  size_t better_solns_total_swap_diff = 0;
  size_t better_solns_percent_decr_total = 0;
  size_t worse_solns = 0;
  size_t worse_solns_swaps = 0;
  size_t worse_solns_known_swaps = 0;
  size_t worse_solns_total_swap_diff = 0;
  size_t worse_solns_percent_incr_total = 0;

  void add(size_t known_size, size_t calc_size) {
    if (known_size == calc_size) {
      ++equivalent_solns;
      equivalent_solns_swaps += known_size;
      return;
    }
    if (calc_size < known_size) {
      ++better_solns;
      better_solns_swaps += calc_size;
      better_solns_known_swaps += known_size;
      const auto decr = known_size - calc_size;
      better_solns_total_swap_diff += decr;
      better_solns_percent_decr_total += (decr * 100) / known_size;
      return;
    }
    ++worse_solns;
    worse_solns_swaps += calc_size;
    worse_solns_known_swaps += known_size;
    const auto incr = calc_size - known_size;
    worse_solns_total_swap_diff += incr;
    worse_solns_percent_incr_total += (incr * 100) / known_size;
  }

  std::string str() const {
    std::stringstream ss;
    size_t good_soln_av_decr = 0;
    if (better_solns > 0) {
      good_soln_av_decr = better_solns_percent_decr_total / better_solns;
    }
    size_t bad_soln_av_incr = 0;
    if (worse_solns > 0) {
      bad_soln_av_incr = worse_solns_percent_incr_total / worse_solns;
    }

    ss << "[" << equivalent_solns << " equal (" << equivalent_solns_swaps
       << "); " << better_solns << " BETTER (" << better_solns_swaps << " vs "
       << better_solns_known_swaps << "): av " << good_soln_av_decr
       << "% decr\n"
       << worse_solns << " WORSE (" << worse_solns_swaps << " vs "
       << worse_solns_known_swaps << "): av " << bad_soln_av_incr << "% incr]";
    return ss.str();
  }
};
}  // namespace

static void check_overall_percentage_improvement(
    unsigned total_number_of_problems, unsigned total_calc_swaps,
    unsigned total_orig_swaps, double expected_percentage) {
  const double actual_decrease =
      100.0 - (100.0 * total_calc_swaps) / (double)total_orig_swaps;
  if (std::abs(actual_decrease - expected_percentage) < 1e-4) {
    return;
  }
  INFO(
      "Solved " << total_number_of_problems
                << " problems; known solutions have total swaps "
                << total_orig_swaps << ". We calculated " << total_calc_swaps
                << ", giving percentage decrease " << actual_decrease
                << ". But we expected " << expected_percentage);
  CHECK(false);
}

namespace {
struct Summary {
  std::string str;
  unsigned total_calc_swaps;
  unsigned total_orig_swaps;
  unsigned total_number_of_problems;

  Summary(
      const vector<std::string>& encoded_swap_sequences, BestTsaTester& tester)
      : total_calc_swaps(0), total_orig_swaps(0), total_number_of_problems(0) {
    FixedSeqsStats stats;
    for (const auto& code_str : encoded_swap_sequences) {
      const DecodedProblemData data(code_str);
      const auto known_size = data.swaps.size();
      REQUIRE(known_size > 0);
      try {
        const auto calc_soln_size = tester.get_checked_solution_size(data);
        stats.add(known_size, calc_soln_size);
        total_calc_swaps += calc_soln_size;
        total_orig_swaps += known_size;
        ++total_number_of_problems;
      } catch (const std::exception& e) {
        INFO(
            "Swap seq encoding string '"
            << code_str << "'\n...encoded " << data.swaps.size() << " swaps, "
            << data.vertex_mapping.size() << " tokens on "
            << data.number_of_vertices
            << " vertices. Gave error: " << e.what());
        REQUIRE(false);
      }
    }
    str = stats.str();
  }

  void check_overall_improvement(double expected_percentage) const {
    check_overall_percentage_improvement(
        total_number_of_problems, total_calc_swaps, total_orig_swaps,
        expected_percentage);
  }
};
}  // namespace

SCENARIO("Best TSA : solve problems from fixed swap sequences") {
  const FixedSwapSequences sequences;
  BestTsaTester tester;

  const Summary full_seqs_summary(sequences.full, tester);
  CHECK(full_seqs_summary.total_number_of_problems == 453);
  CHECK(
      full_seqs_summary.str ==
      "[248 equal (6088); 104 BETTER (4645 vs 4979): av 7% decr\n"
      "101 WORSE (5893 vs 5451): av 8% incr]");

  // The fixed swap sequences have been optimised quite a lot already,
  // so are probably quite close to optimal.
  full_seqs_summary.check_overall_improvement(-0.653832);

  const Summary partial_seqs_summary(sequences.partial, tester);
  CHECK(partial_seqs_summary.total_number_of_problems == 755);
  CHECK(
      partial_seqs_summary.str ==
      "[455 equal (6487); 165 BETTER (7044 vs 7457): av 7% decr\n"
      "135 WORSE (9124 vs 8604): av 6% incr]");

  partial_seqs_summary.check_overall_improvement(-0.474543);
}

// Now we want to solve complete problems; this is one of
// our most important tests. It is a bit silly
// to put problems with 5 vertices and problems with
// 50 vertices in the same test. Therefore, we crudely sort by length of
// encoding string, which is roughly "problem size",
// and distribute the final statistics amongst a number of categories
// based upon problem size.
namespace {
class StatisticsGrouper {
 public:
  StatisticsGrouper(
      unsigned number_of_messages,
      const vector<unsigned>& sorted_problem_sizes) {
    REQUIRE(number_of_messages >= 3);
    REQUIRE(sorted_problem_sizes.size() >= 5 * number_of_messages);
    REQUIRE(sorted_problem_sizes[0] >= 5);
    m_stats.resize(number_of_messages);
    m_problem_size_boundaries.resize(number_of_messages);
    const unsigned step = sorted_problem_sizes.size() / number_of_messages;
    for (unsigned ii = 0; ii + 1 < number_of_messages; ++ii) {
      m_problem_size_boundaries[ii] = sorted_problem_sizes[(ii + 1) * step];
    }
    m_problem_size_boundaries.back() = sorted_problem_sizes.back() + 1;
  }

  void add(
      const std::string& problem_str,
      const DecodedArchitectureData& arch_data) {
    unsigned allowed_index = 0;
    for (unsigned index = 0; index < m_problem_size_boundaries.size();
         ++index) {
      if (problem_str.size() <= m_problem_size_boundaries[index]) {
        allowed_index = index;
        break;
      }
    }
    // Now we know which category it's in, so do the calculation
    auto& stats = m_stats[allowed_index];
    const DecodedProblemData data(
        problem_str, DecodedProblemData::RequireContiguousVertices::NO);
    const auto known_size = data.swaps.size();
    REQUIRE(known_size > 0);
    try {
      const auto calc_soln_size =
          m_tester.get_checked_solution_size(data, arch_data);
      stats.add(known_size, calc_soln_size);
      m_total_calc_swaps += calc_soln_size;
      m_total_orig_swaps += known_size;
      ++m_total_number_of_problems;
    } catch (const std::exception& e) {
      INFO(
          "Swap seq encoding string '" << problem_str << "'\n...encoded "
                                       << data.swaps.size()
                                       << " swaps, error: " << e.what());
      CHECK(false);
    }
  }

  vector<std::string> get_final_messages() const {
    vector<std::string> messages(m_stats.size());
    for (unsigned ii = 0; ii < m_stats.size(); ++ii) {
      messages[ii] = m_stats[ii].str();
    }
    return messages;
  }

  void check_overall_improvement(double expected_percentage) const {
    check_overall_percentage_improvement(
        m_total_number_of_problems, m_total_calc_swaps, m_total_orig_swaps,
        expected_percentage);
  }

 private:
  unsigned m_total_calc_swaps = 0;
  unsigned m_total_orig_swaps = 0;
  unsigned m_total_number_of_problems = 0;
  BestTsaTester m_tester;
  vector<FixedSeqsStats> m_stats;
  vector<unsigned> m_problem_size_boundaries;
};
}  // namespace

SCENARIO("Best TSA : solve complete problems") {
  const FixedCompleteSolutions complete_solutions;

  // For a good test, very different problems should not be amalgamated
  // in the statistics. Thus we determine the different categories using length
  // of encoding string, which presumably roughly corresponds to "problem size"
  // and problem hardness.
  const vector<std::string> expected_messages{
      "[210 equal (1018); 19 BETTER (84 vs 111): av 24% decr\n"
      "2 WORSE (19 vs 15): av 26% incr]",

      "[145 equal (1822); 39 BETTER (451 vs 525): av 13% decr\n"
      "17 WORSE (269 vs 242): av 11% incr]",

      "[58 equal (1619); 122 BETTER (3465 vs 3832): av 9% decr\n"
      "34 WORSE (1321 vs 1232): av 6% incr]",

      "[18 equal (1382); 114 BETTER (8322 vs 8856): av 5% decr\n"
      "83 WORSE (6875 vs 6457): av 5% incr]",

      "[8 equal (1470); 164 BETTER (25183 vs 27141): av 6% decr\n"
      "44 WORSE (8722 vs 8384): av 3% incr]"};

  vector<unsigned> problem_sizes;
  for (const auto& entry : complete_solutions.solutions) {
    REQUIRE(entry.second.size() >= 2);
    // The first string encodes the edges in that architecture,
    // rather than a problem.
    for (unsigned ii = 1; ii < entry.second.size(); ++ii) {
      problem_sizes.push_back(entry.second[ii].size());
    }
  }
  std::sort(problem_sizes.begin(), problem_sizes.end());
  StatisticsGrouper grouper(expected_messages.size(), problem_sizes);

  // Now go through the problems, let the grouper object collate the stats
  // appropriately
  for (const auto& entry : complete_solutions.solutions) {
    const DecodedArchitectureData arch_data(entry.second[0]);
    for (unsigned ii = 1; ii < entry.second.size(); ++ii) {
      grouper.add(entry.second[ii], arch_data);
    }
  }
  const auto calc_messages = grouper.get_final_messages();
  REQUIRE(calc_messages.size() == expected_messages.size());
  for (unsigned ii = 0; ii < calc_messages.size(); ++ii) {
    INFO("for message[" << ii << "]: ");
    CHECK(calc_messages[ii] == expected_messages[ii]);
  }
  // A positive result is good; the fixed complete problems are DIRECTLY
  // comparing our TSA with the solver used to generate them.
  grouper.check_overall_improvement(3.25087);
}

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
