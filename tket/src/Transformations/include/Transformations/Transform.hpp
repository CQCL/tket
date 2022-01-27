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

#include <functional>

#include "Circuit/Circuit.hpp"

namespace tket {

class Transform {
 public:
  typedef std::function<bool(Circuit&)> Transformation;
  typedef std::function<unsigned(const Circuit&)> Metric;

  // the actual transformation to be applied
  // performs transformation in place and returns true iff made some change
  Transformation apply;  // this would ideally be `const`, but that deletes the
                         // copy assignment operator for Transform.

  explicit Transform(const Transformation& trans) : apply(trans) {}

  friend Transform operator>>(const Transform& lhs, const Transform& rhs);
};

namespace Transforms {

// identity Transform (does nothing to Circuit)
inline const Transform id = Transform([](const Circuit&) {
  return false;
});  // returns `false` as it does not change the Circuit in any way

}  // namespace Transforms

}  // namespace tket
