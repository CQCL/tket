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
#include <tkwsm/EndToEndWrappers/MainSolver.hpp>

namespace tket {
namespace WeightedSubgraphMonomorphism {

struct CheckedSolution;

// Take a solution which has already run,
// then re-solve, breaking into several steps by setting max_iterations.
// The end result should be identical.
class ResumedSolutionChecker {
 public:
  // Check that we end up with an identical solution
  // even when we stop and resume multiple times.
  void check(
      const CheckedSolution& solution, const GraphEdgeWeights& pdata,
      const GraphEdgeWeights& tdata, MainSolverParameters solver_params);

 private:
  unsigned m_min_number_of_iterations = 10;
  unsigned m_max_number_of_iterations = 10000;
  unsigned m_number_of_chunks = 5;
  unsigned m_remaining_problems = 10;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
