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

#include <string>
#include <vector>

namespace tket {
namespace tsa_internal {
namespace tests {

/** Random problems are easy to generate, but good ones maybe not so easy.
 *  Another way to generate possibly better problems is to take, first,
 *  a random problem and solve it with a reasonable TSA.
 *  We then record the sequence of swaps generated,
 *  and use the sequence to GENERATE a problem
 *  (taking the graph to have only edges which arise from mentioned vertex
 *  swaps and no others, and taking the desired vertex mapping to be
 *  exactly that which arises from performing the swaps).
 *  We then solve this new problem with our TSA and compare the number of swaps.
 *
 *  This has the benefit of providing a solution to compare against
 *  (the original swaps), which is also perhaps quite hard to improve upon,
 *  because at least one other TSA did not do any better.
 *
 *  Of course it is not a direct comparison of our TSA with others, because
 *
 *    (1): We obtained these swap sequences by removing unused edges in a
 * returned solution. This actually changes the problem, so it is possible that
 *        the returned solution would change if presented with this new problem.
 *        (Although it would be most elegant mathematically if this did not
 * occur, it seems hard to enforce it in an algorithm. There is not much benefit
 * in doing so, so it seems unlikely that it would arise "by chance". Even if it
 * did, proving that such a property did hold would be hard).
 *
 *    (2): Vertex relabelling also changes the problem, even though it is
 * "isomorphic". It seems very unlikely that an algorithm would always return
 *         isomorphic solutions to isomorphic problems. (Even if it were an
 * optimal algorithm, the solutions may not be unique even up to isomorphism).
 *
 *    (3): The TSA may be non-deterministic, due to RNGs. (Our algorithm is
 * deterministic, however, since we deliberately set all RNGs to a default seed
 * before use).
 *
 *    (4): The returned swaps have already been run through SwapListOptimiser
 * passes to reduce them.
 */
struct FixedSwapSequences {
  /* Encoding swap sequences as strings rather than inside a vector<Swap>
   *  should give smaller C++ and .obj files.
   *  For convenience, the vertex numbers in each problem should be
   * {0,1,2,...,n} with no gaps. Also for convenience, sorted by string length;
   * shorter strings are usually "simpler" problems.
   *
   * "Full" sequences came from problems where every vertex had a token.
   * "Partial" sequences came from problems where only some vertics had a token.
   * Thus, the vertices which did initially have tokens are also specified;
   * for a fair test, this is essential as it may enable reductions
   * which would be invalid in the "full" case.
   *
   * Note that some sequences currently give errors with the best TSA.
   * It is due to disconnected architectures, which can cause errors
   * (although not always). This is a bug which should be fixed, although
   * every architecture we use in practice should be connected.
   */

  std::vector<std::string> full;
  std::vector<std::string> partial;
  std::vector<std::string> full_with_errors;
  std::vector<std::string> partial_with_errors;

  /** Upon construction, the fixed sequences will all be set. */
  FixedSwapSequences();
};

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
