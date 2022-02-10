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
#include "FullPlacementResult.hpp"
#include "Utils/UnitID.hpp"

namespace tket {

class Architecture;
class Circuit;

struct CalculatedPlacementMap {
  /** For testing, it's helpful to have the full internal result. */
  FullPlacementResult full_placement_result;
  std::map<Qubit, Node> placement_map;

  /** Extra algorithmic parameters to configure the placement. */
  struct Parameters {
    /** The timeout in milliseconds. */
    unsigned timeout_ms = 10000;
  };

  CalculatedPlacementMap(
      const Circuit& circ, const Architecture& arch,
      const Parameters& parameters);
};

}  // namespace tket
