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

#include "tket/Transformations/StandardSquash.hpp"

#include <memory>

#include "tket/Circuit/DAGDefs.hpp"
#include "tket/Gate/Rotation.hpp"
#include "tket/OpType/OpTypeInfo.hpp"
#include "tket/Transformations/BasicOptimisation.hpp"
#include "tket/Transformations/SingleQubitSquash.hpp"
#include "tket/Utils/Expression.hpp"

namespace tket {

namespace Transforms {

StandardSquasher::StandardSquasher(
    const OpTypeSet &singleqs, const Func &tk1_replacement)
    : singleqs_(singleqs),
      squash_fn_(tk1_replacement),
      combined_(),
      phase_(0.) {
  for (OpType ot : singleqs_) {
    if (!is_single_qubit_type(ot))
      throw BadOpType(
          "OpType given to standard_squash is not a single qubit gate", ot);
  }
}

bool StandardSquasher::accepts(OpType optype) const {
  return (singleqs_.find(optype) != singleqs_.end()) &&
         !is_projective_type(optype);
}

void StandardSquasher::append(Gate_ptr gp) {
  std::vector<Expr> angs = gp->get_tk1_angles();
  combined_.apply(Rotation(OpType::Rz, angs.at(2)));
  combined_.apply(Rotation(OpType::Rx, angs.at(1)));
  combined_.apply(Rotation(OpType::Rz, angs.at(0)));
  phase_ += angs.at(3);
}

std::pair<Circuit, Gate_ptr> StandardSquasher::flush(
    std::optional<Pauli>) const {
  auto [a, b, c] = combined_.to_pqp(OpType::Rz, OpType::Rx);
  Circuit replacement = squash_fn_(c, b, a);
  BGL_FORALL_VERTICES(rv, replacement.dag, DAG) {
    OpType v_type = replacement.get_OpType_from_Vertex(rv);
    if (!is_boundary_q_type(v_type) &&
        singleqs_.find(v_type) == singleqs_.end()) {
      throw BadOpType(
          "tk1_replacement given to standard_squash "
          "does not preserve gate set",
          v_type);
    }
  }
  replacement.add_phase(phase_);
  return {replacement, nullptr};
}

void StandardSquasher::clear() {
  combined_ = Rotation();
  phase_ = 0.;
}

std::unique_ptr<AbstractSquasher> StandardSquasher::clone() const {
  return std::make_unique<StandardSquasher>(*this);
}

static bool standard_squash(
    Circuit &circ, const OpTypeSet &singleqs,
    const std::function<Circuit(const Expr &, const Expr &, const Expr &)>
        &tk1_replacement,
    bool always_squash_symbols) {
  auto squasher = std::make_unique<StandardSquasher>(singleqs, tk1_replacement);
  return SingleQubitSquash(
             std::move(squasher), circ, false, always_squash_symbols)
      .squash();
}

Transform squash_factory(
    const OpTypeSet &singleqs,
    const std::function<Circuit(const Expr &, const Expr &, const Expr &)>
        &tk1_replacement,
    bool always_squash_symbols) {
  return Transform([=](Circuit &circ) {
    return standard_squash(
        circ, singleqs, tk1_replacement, always_squash_symbols);
  });
}

}  // namespace Transforms

}  // namespace tket
