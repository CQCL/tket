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

#include "Gate/GatePtr.hpp"
#include "OpType/OpType.hpp"
#include "SingleQubitSquash.hpp"

namespace tket {

namespace Transforms {

/**
 * @brief Implements the AbstractSquasher interface for SingleQubitSquash
 *
 * The PQP Squasher squashes chains of single qubit gates to minimal sequences
 * of any two type of rotations (Rx, Ry or Rz), using Euler angle
 * decompositions.
 *
 * Decompositions can either be strict or smart. Strict decompositions
 * will always decompose to P-Q-P, whereas smart will decompose to P-Q-P
 * or Q-P-Q and will always try to commute P or Q through multi-qubit gates.
 * For smart decomposition, swapping P and Q gives the same result.
 *
 * Note that even in strict mode, if a rotation in the PQP decomposition
 * has angle 0, it will be omitted.
 *
 * The squash is made backwards, so that rotations get pushed towards the front.
 * This was chosen to be compatible with the `commute_through_multis` pass,
 * which also proceeds backwards. There are several reasons to prefer commuting
 * rotations to the front rather than the back
 *  - Noise heuristic: CX gate distribute the errors out, so it makes sense to
 *    delay them as much as possible
 *  - For contextual optimisation, it also makes sense to move Clifford
 *    operations (i.e. CX) to the back, as this can be then removed
 *  - Some initial benchmarking by Seyon seemed to show that commuting to the
 *    front performed better, but this might be an artefact of the benchmarking
 *    circuits used.
 * Note finally there is also an argument for commuting non-Clifford rotations
 * towards the back, as this makes generating automatic assertions easier (not
 * implemented at the time of writing).
 *
 * The PQPSquasher shouldn't have to know if the squash is made forwards
 * or backwards. Here, however, `reversed` is needed in `fixup_angles` for
 * consistency, so that forward and backward passes squash to the same normal
 * form
 */
class PQPSquasher : public AbstractSquasher {
 public:
  PQPSquasher(
      OpType p, OpType q, bool smart_squash = true, bool reversed = false);

  void clear() override;

  bool accepts(Gate_ptr gp) const override;

  void append(Gate_ptr gp) override;

  std::pair<Circuit, Gate_ptr> flush(
      std::optional<Pauli> commutation_colour = std::nullopt) const override;

  std::unique_ptr<AbstractSquasher> clone() const override;

 private:
  OpType p_;
  OpType q_;
  bool smart_squash_;
  bool reversed_;
  std::vector<Gate_ptr> rotation_chain;
};

}  // namespace Transforms

}  // namespace tket