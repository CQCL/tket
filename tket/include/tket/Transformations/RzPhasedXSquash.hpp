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

#include "PQPSquash.hpp"
#include "tket/Gate/GatePtr.hpp"

namespace tket {

namespace Transforms {

/**
 * @brief Squash Rx, Rz gates into PhasedX and Rz gates
 *
 */
class RzPhasedXSquasher : public PQPSquasher {
 public:
  RzPhasedXSquasher(bool reversed = false);
  std::pair<Circuit, Gate_ptr> flush(
      std::optional<Pauli> commutation_colour = std::nullopt) const;
};

}  // namespace Transforms

}  // namespace tket