// Copyright 2019-2023 Cambridge Quantum Computing
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

#include <cstdint>
#include <limits>
#include <string>

namespace tket {
namespace tsa_internal {
namespace tests {

struct MinMaxAv {
  std::size_t min = std::numeric_limits<std::size_t>::max();
  std::size_t max = 0;
  std::size_t total = 0;

  void add(std::size_t result);
};

struct PartialTsaStatistics {
  std::size_t number_of_problems = 0;
  std::size_t total_of_L = 0;
  std::size_t total_of_L_decreases = 0;
  std::size_t total_number_of_tokens = 0;
  std::size_t total_number_of_swaps = 0;

  MinMaxAv l_decrease_percentages;

  // The "power" of a swap sequence (with given token configuration)
  // is defined to be  (decrease in L)/(number of swaps).
  // Thus, it's always between 0 and 2 (if all swaps make progress).
  // However, we multiply by 50, to make the power between 0 and 100%.
  MinMaxAv powers;

  void add_problem_result(
      std::size_t initial_L, std::size_t final_L, std::size_t tokens,
      std::size_t swaps);

  std::string str(std::size_t number_of_problems) const;
};

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
