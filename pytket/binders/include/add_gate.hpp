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
#include <nanobind/nanobind.h>

#include <variant>
#include <vector>

#include "nanobind-stl.hpp"
#include "tket/Circuit/Circuit.hpp"
#include "typecast.hpp"
namespace nb = nanobind;

namespace tket {

Circuit *add_gate_method_any(
    Circuit *circ, const Op_ptr &op,
    const std::variant<std::vector<unsigned>, std::vector<UnitID>> &args,
    const nb::kwargs &kwargs);

}  // namespace tket
