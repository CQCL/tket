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

#include "ZX/Rewrite.hpp"

namespace tket {

namespace zx {

Rewrite::Rewrite(const RewriteFun &fun) : apply(fun) {}

Rewrite Rewrite::sequence(const std::vector<Rewrite> &rvec) {
  return Rewrite([=](ZXDiagram &diag) {
    bool success = false;
    for (std::vector<Rewrite>::const_iterator it = rvec.begin();
         it != rvec.end(); ++it) {
      success = it->apply(diag) || success;
    }
    return success;
  });
}

Rewrite Rewrite::repeat(const Rewrite &rw) {
  return Rewrite([=](ZXDiagram &diag) {
    bool success = false;
    while (rw.apply(diag)) success = true;
    return success;
  });
}

Rewrite Rewrite::repeat_with_metric(
    const Rewrite &rw, const Rewrite::Metric &eval) {
  return Rewrite([=](ZXDiagram &diag) {
    bool success = false;
    unsigned currentVal = eval(diag);
    ZXDiagram newDiag = diag;
    rw.apply(newDiag);
    unsigned newVal = eval(newDiag);
    while (newVal < currentVal) {
      currentVal = newVal;
      success = true;
      rw.apply(newDiag);
      newVal = eval(newDiag);
    }
    if (success) {
      diag = newDiag;
    }
    return success;
  });
}

Rewrite Rewrite::repeat_while(const Rewrite &cond, const Rewrite &body) {
  return Rewrite([=](ZXDiagram &diag) {
    bool success = false;
    while (cond.apply(diag)) {
      success = true;
      body.apply(diag);
    }
    return success;
  });
}

}  // namespace zx

}  // namespace tket
