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

#include <memory>

#include "Gate.hpp"
#include "Ops/OpPtr.hpp"

namespace tket {

typedef std::shared_ptr<const Gate> Gate_ptr;

/**
 * Cast a general `Op` (of gate type) to a `Gate`.
 *
 * @throws NotValid if @p op is not a gate.
 */
Gate_ptr as_gate_ptr(Op_ptr op);

}  // namespace tket
