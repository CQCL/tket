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

#include "Transform.hpp"

namespace tket {

namespace Transforms {

// extends PhaseGadgets by a qubit on identifying a pair of CXs around it
// Expects: CX, PhaseGadget, and any single-qubit gates
// Produces: CX, PhaseGadget, and any single-qubit gates
Transform smash_CX_PhaseGadgets();

// tries to match up ports on adjacent PhaseGadgets to enable maximal
// annihilation after synthesis (ignoring any intervening gates) Expects: Any
// gates Produces: The same gate set
Transform align_PhaseGadgets();

}  // namespace Transforms

}  // namespace tket
