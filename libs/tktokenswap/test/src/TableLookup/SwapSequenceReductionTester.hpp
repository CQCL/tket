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

#pragma once

#include <string>
#include <tktokenswap/SwapListOptimiser.hpp>
#include <tktokenswap/SwapListTableOptimiser.hpp>

#include "../TestUtils/DecodedProblemData.hpp"

namespace tket {
namespace tsa_internal {
namespace tests {

/** Directly test the results of table reductions on fixed swap sequences. */
class SwapSequenceReductionTester {
 public:
  struct Options {
    bool optimise_initial_segment_only;
  };

  // Reduces the sequence of swaps, checks it, and returns the size.
  std::size_t get_checked_solution_size(
      const DecodedProblemData& problem_data,
      const DecodedArchitectureData& architecture_data, const Options& options);

  std::size_t get_checked_solution_size(
      const DecodedProblemData& problem_data, const Options& options);

 private:
  SwapListOptimiser m_general_optimiser;
  // SwapList m_raw_swap_list;
};

struct SequenceReductionStats {
  std::size_t problems;
  std::size_t reduced_problems;
  std::size_t total_original_swaps;

  // This only includes problems where the number of swaps strictly decreased
  // after table reduction.
  std::size_t total_original_swaps_for_reduced_problems;

  // This is the sum of "reduced_swaps" passed in, over all problems (including
  // those where there was no decrease).
  std::size_t total_reduced_swaps;

  SequenceReductionStats();

  void add_solution(std::size_t original_swaps, std::size_t reduced_swaps);

  std::string str() const;
};

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
