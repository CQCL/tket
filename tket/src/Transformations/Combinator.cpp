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

#include "Transform.hpp"

namespace tket {

Transform operator>>(const Transform &lhs, const Transform &rhs) {
  std::vector<Transform> l = {lhs, rhs};
  return Transform::sequence(l);
}

Transform Transform::sequence(std::vector<Transform> &tvec) {
  return Transform([=](Circuit &circ) {
    bool success = false;
    for (std::vector<Transform>::const_iterator it = tvec.begin();
         it != tvec.end(); ++it) {
      success = it->apply(circ) || success;
    }
    return success;
  });
}

Transform Transform::repeat(const Transform &trans) {
  return Transform([=](Circuit &circ) {
    bool success = false;
    while (trans.apply(circ)) success = true;
    return success;
  });
}

Transform Transform::repeat_with_metric(
    const Transform &trans, const Metric &eval) {
  return Transform([=](Circuit &circ) {
    bool success = false;
    int currentVal = eval(circ);
    Circuit *currentCircuit = &circ;
    Circuit newCircuit = circ;
    trans.apply(newCircuit);
    int newVal = eval(newCircuit);
    while (newVal < currentVal) {
      currentCircuit = &newCircuit;
      currentVal = newVal;
      success = true;
      trans.apply(newCircuit);
      newVal = eval(newCircuit);
    }
    if (&circ != currentCircuit) circ = *currentCircuit;
    return success;
  });
}

Transform Transform::repeat_while(
    const Transform &cond, const Transform &body) {
  return Transform([=](Circuit &circ) {
    bool success = false;
    while (cond.apply(circ)) {
      success = true;
      body.apply(circ);
    }
    return success;
  });
}

}  // namespace tket
