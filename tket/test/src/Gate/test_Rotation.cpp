// Copyright 2019-2023 Cambridge Quantum Computing
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

#include <catch2/catch_test_macros.hpp>
#include <vector>

#include "../testutil.hpp"
#include "tket/Circuit/Circuit.hpp"
#include "tket/Gate/Rotation.hpp"
#include "tket/OpType/OpType.hpp"
#include "tket/Utils/Expression.hpp"

namespace tket {

SCENARIO("xyx decomposition") {
  const std::vector<std::tuple<Expr, Expr, Expr>> angles = {
      {0.2, 0.3, 0.4},
      {0.4, 0.8, 1.4},
      {0.4, 0.8, 0.6},
  };
  for (const std::tuple<Expr, Expr, Expr> &abc : angles) {
    Expr a = std::get<0>(abc);
    Expr b = std::get<1>(abc);
    Expr c = std::get<2>(abc);
    Rotation rxa(OpType::Rx, a);
    Rotation ryb(OpType::Ry, b);
    Rotation rxc(OpType::Rx, c);
    Rotation r(rxa);
    r.apply(ryb);
    r.apply(rxc);
    std::tuple<Expr, Expr, Expr> abc1 = r.to_pqp(OpType::Rx, OpType::Ry);
    Expr a1 = std::get<0>(abc1);
    Expr b1 = std::get<1>(abc1);
    Expr c1 = std::get<2>(abc1);
    Circuit circ(1);
    circ.add_op<unsigned>(OpType::Rx, a, {0});
    circ.add_op<unsigned>(OpType::Ry, b, {0});
    circ.add_op<unsigned>(OpType::Rx, c, {0});
    Circuit circ1(1);
    circ1.add_op<unsigned>(OpType::Rx, a1, {0});
    circ1.add_op<unsigned>(OpType::Ry, b1, {0});
    circ1.add_op<unsigned>(OpType::Rx, c1, {0});
    REQUIRE(test_unitary_comparison(circ, circ1));
  }
}

}  // namespace tket