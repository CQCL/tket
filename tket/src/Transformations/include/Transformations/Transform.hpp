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
#include <memory>

#include "Circuit/Circuit.hpp"

namespace tket {

typedef struct {
  unit_bimap_t initial;
  unit_bimap_t final;
} unit_bimaps_t;

class Transform {
 public:
  typedef std::function<bool(Circuit&, std::shared_ptr<unit_bimaps_t>)>
      Transformation;
  typedef std::function<bool(Circuit&)> SimpleTransformation;
  typedef std::function<unsigned(const Circuit&)> Metric;

  Transformation apply_fn;

  explicit Transform(const Transformation& trans) : apply_fn(trans) {}

  explicit Transform(const SimpleTransformation& trans)
      : apply_fn([=](Circuit& circ, std::shared_ptr<unit_bimaps_t>) {
          return trans(circ);
        }) {}

  bool apply(Circuit& circ) const { return apply_fn(circ, nullptr); }

  friend Transform operator>>(const Transform& lhs, const Transform& rhs);
};

namespace Transforms {

// identity Transform (does nothing to Circuit)
inline const Transform id =
    Transform([](Circuit&, std::shared_ptr<unit_bimaps_t>) {
      return false;
    });  // returns `false` as it does not change the Circuit in any way

}  // namespace Transforms

}  // namespace tket
