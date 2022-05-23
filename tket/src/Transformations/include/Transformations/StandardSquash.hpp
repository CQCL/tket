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

#include "Gate/GatePtr.hpp"
#include "Gate/Rotation.hpp"
#include "OpType/OpTypeFunctions.hpp"
#include "SingleQubitSquash.hpp"
#include "Utils/Expression.hpp"

namespace tket {

namespace Transforms {

/**
 * @brief Implements the AbstractSquasher interface for SingleQubitSquash
 *
 * The StandardSquasher squashes chains of single qubit gates to
 * the circuit given by the tk1_replacment function passed as parameter.
 *
 * At the moment, it does not commute anything through multi-qubit gates.
 */
class StandardSquasher : public AbstractSquasher {
 public:
  using Func = std::function<Circuit(const Expr &, const Expr &, const Expr &)>;

  StandardSquasher(const OpTypeSet &singleqs, const Func &tk1_replacement);

  bool accepts(Gate_ptr gp) const override;

  void append(Gate_ptr gp) override;

  std::pair<Circuit, Gate_ptr> flush(
      std::optional<Pauli> = std::nullopt) const override;

  void clear() override;

  std::unique_ptr<AbstractSquasher> clone() const override;

 private:
  const OpTypeSet singleqs_;
  const Func squash_fn_;
  Rotation combined_;
  Expr phase_;
};

}  // namespace Transforms

}  // namespace tket
