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

#include <vector>

#include "Utils/Json.hpp"

namespace tket {

/** Type of a wire in a circuit or input to an op */
enum class EdgeType {
  /**
   * A wire carrying quantum information, corresponding to some allocated @ref
   * Qubit. Since these are persistent, every node in the DAG (except for an
   * input or output) has the same number of Quantum input and output wires.
   */
  Quantum,

  /**
   * A wire carrying classical information, corresponding to some allocated
   * @ref Bit. Since these are persistent, every node in the DAG (except for
   * an input or output) has the same number of Classical input and output
   * wires.
   */
  Classical,

  /**
   * A wire carrying a bit of classical information from a classical output
   * port of one op to a classical input port of another, not corresponding to
   * any allocated @ref Bit.
   */
  Boolean
};

NLOHMANN_JSON_SERIALIZE_ENUM(
    EdgeType, {
                  {EdgeType::Quantum, "Q"},
                  {EdgeType::Classical, "C"},
                  {EdgeType::Boolean, "B"},
              });

typedef std::vector<EdgeType> op_signature_t;

}  // namespace tket
