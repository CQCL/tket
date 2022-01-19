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

#include "Circuit/DAGDefs.hpp"
#include "Gate/Rotation.hpp"
#include "SingleQubitSquash.hpp"

namespace tket {

/**
 * @brief Implements the Squasher interface for SingleQubitSquash
 *
 * The StandardSquasher squashes chains of single qubit gates to
 * the circuit given by the tk1_replacment function passed as parameter.
 *
 * At the moment, it does not commute anything through multi-qubit gates.
 */
class StandardSquasher {
 private:
  using Func = std::function<Circuit(const Expr &, const Expr &, const Expr &)>;

 public:
  StandardSquasher(const OpTypeSet &singleqs, const Func &tk1_replacement)
      : singleqs(singleqs), squash_fn(tk1_replacement), combined() {
    for (OpType ot : singleqs) {
      if (!is_single_qubit_type(ot))
        throw NotValid(
            "OpType given to standard_squash is not a single qubit gate");
    }
  }

  bool accepts(OpType type) const {
    return (singleqs.find(type) != singleqs.end()) && !is_projective_type(type);
  }

  void clear() { combined = Rotation(); }

  void append(Gate_ptr gate) {
    std::vector<Expr> angs = gate->get_tk1_angles();
    combined.apply(Rotation(OpType::Rz, angs.at(2)));
    combined.apply(Rotation(OpType::Rx, angs.at(1)));
    combined.apply(Rotation(OpType::Rz, angs.at(0)));
  }

  std::pair<Circuit, Gate_ptr> flush(
      std::optional<Pauli> = std::nullopt) const {
    auto [a, b, c] = combined.to_pqp(OpType::Rz, OpType::Rx);
    Circuit replacement = squash_fn(c, b, a);
    BGL_FORALL_VERTICES(rv, replacement.dag, DAG) {
      OpType v_type = replacement.get_OpType_from_Vertex(rv);
      if (!is_boundary_q_type(v_type) &&
          singleqs.find(v_type) == singleqs.end()) {
        throw NotValid(
            "tk1_replacement given to standard_squash "
            "does not preserve gate set");
      }
    }
    return {replacement, nullptr};
  }

 private:
  const OpTypeSet &singleqs;
  const Func &squash_fn;
  Rotation combined;
};

static bool standard_squash(
    Circuit &circ, const OpTypeSet &singleqs,
    const std::function<Circuit(const Expr &, const Expr &, const Expr &)>
        &tk1_replacement) {
  StandardSquasher squasher(singleqs, tk1_replacement);
  return SingleQubitSquash(squasher, false).squash(circ);
}

Transform Transform::squash_factory(
    const OpTypeSet &singleqs,
    const std::function<Circuit(const Expr &, const Expr &, const Expr &)>
        &tk1_replacement) {
  return Transform([=](Circuit &circ) {
    return standard_squash(circ, singleqs, tk1_replacement);
  });
}

}  // namespace tket