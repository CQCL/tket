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

#pragma once

namespace tket {
namespace tsa_internal {
namespace tests {

/** If we want to use the same adjustable parameters across all
 * TokenSwapping tests simultaneously, put them here.
 */
struct TSGlobalTestParameters {
  /** Running all the token swapping tests can take ~30 seconds
   * on an ordinary laptop. Set this to false in order to test
   * a smaller set.
   */
  bool run_long_tests;

  // TSGlobalTestParameters() : run_long_tests(true) {}
  TSGlobalTestParameters() : run_long_tests(false) {}
};

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
