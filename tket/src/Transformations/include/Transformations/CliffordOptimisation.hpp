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

// mutli-qubit patterns that decrease the CX count
// inserting swaps can sometimes cause errors elsewhere (e.g. routing), so
// they can be turned off
Transform multiq_clifford_replacement(bool allow_swaps = false);

// copies Z through the target of a CX and X through the control
Transform copy_pi_through_CX();

// These apply some of the Clifford rules in the paper "Optimising Clifford
// Circuits with Quantomatic" All of these expect and produce CX, Z, X, S, V,
// and other single qubit gates

Transform singleq_clifford_sweep();

}  // namespace Transforms

}  // namespace tket
