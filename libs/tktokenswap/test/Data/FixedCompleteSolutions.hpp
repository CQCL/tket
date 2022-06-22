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

#include <map>
#include <string>
#include <vector>

namespace tket {
namespace tsa_internal {
namespace tests {

/** These store complete solutions to fixed token swapping problems.
 *  This is similar to FixedSwapSequences.
 *  However, it is different in several ways:
 *
 *  (1) The solutions have not been processed further (thus, we
 *    do not expect the swap sequences to be irreducible).
 *    In particular, we did not relabel the vertices of the solutions,
 *    so they will not be contiguous.
 *
 *  (2) The full set of edges passed into the original solver is preserved
 *    (thus, we expect more variety in possible solutions; there may be
 *    more shortcuts making use of different edges).
 *    In particular, all architectures are connected, so there should be
 *    NO errors when running our TSA.
 *
 *  (3) Several real architectures are included.
 *
 *  I have tried to include a reasonable range of architectures
 *  and problem sizes.
 *
 *  Thus, this allows a direct comparison between our TSA
 *  and the one used to generate these solutions, and hopefully will show
 *  improvements more clearly over time.
 *  These are also hopefully more realistic problems.
 *  However, as noted also in FixedSwapSequences, we must remember that:
 *
 *  (a) relabelling vertices will, in most cases, give different solutions
 *  [even though the problems are "isomorphic"]; this is just an unavoidable
 *  consequence of the token swapping problem being hard and, presumably,
 *  often having many "nonisomorphic" optimal solutions [although this hasn't
 *  been precisely defined]. Thus, we can never REALLY do a direct comparison
 *  because we're always going to get small differences just by "chance",
 *  depending upon our vertex labelling;
 *
 *  (b) Many algorithms involve an RNG and hence do not give the same solution
 *  each time (although, our TSAs are careful always to reset the RNG seed,
 *  so should actually be deterministic).
 */
struct FixedCompleteSolutions {
  // KEY: the architecture name
  // VALUE: the problems, encoded as strings; the first element
  //  encodes the complete collection of edges (which cannot be deduced from the
  //  solution swaps because, of course, some edges might be unused). The
  //  remaining elements are the calculated solutions to actual problems, with
  //  the same encoding as in FixedSwapSequences. Thus the tokens are given, but
  //  the vertex mapping is not, since it can be deduced from the swaps as
  //  usual.
  std::map<std::string, std::vector<std::string>> solutions;

  // Fill in all the problem data upon construction.
  FixedCompleteSolutions();
};

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
