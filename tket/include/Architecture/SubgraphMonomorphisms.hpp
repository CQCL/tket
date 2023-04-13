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

#pragma once
#include "ArchitectureMapping.hpp"

namespace tket {

/** A simple wrapper around the Weighted Subgraph Monomorphism code. */
struct SubgraphMonomorphisms {
  struct Parameters {
    unsigned long timeout_ms = 10000;
    unsigned max_number_of_mappings = 1;
  };

  SubgraphMonomorphisms(
      const ArchitectureMapping& pattern_arch_mapping,
      const ArchitectureMapping& target_arch_mapping,
      const Parameters& parameters);

  /** Element[i] is the target vertex number which the pattern vertex number i
   * is mapped to. */
  typedef std::vector<size_t> Mapping;

  /** Only complete, valid mappings are included; they are all distinct from
   * each other. This should be fully deterministic and platform independent
   * (but only if the vertex numbers don't change). To obtain the original
   * qubits etc. from the vertex numbers, use the appropriate
   * ArchitectureMapping. WARN: the ordering of mappings, and the mappings found
   * if not run to completion, are highly liable to change. WARN: relabelling
   * (permuting) the vertex numbers will usually change the ordering
   * (unsurprisingly, since even if it used an inbuilt subgraph isomorphism
   * detector - which it doesn't - there's no canonical ordering), and time
   * taken.
   */
  std::vector<Mapping> mappings;

  /** The actual time taken, in milliseconds. */
  unsigned long long time_taken_ms;
};

}  // namespace tket
