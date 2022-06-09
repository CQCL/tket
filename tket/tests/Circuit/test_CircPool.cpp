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

#include <catch2/catch_test_macros.hpp>

#include "Circuit/CircPool.hpp"
#include "Simulation/CircuitSimulator.hpp"

namespace tket {
namespace test_CircPool {

SCENARIO("CircPool identities are correct") {
  Circuit orig, res;

  GIVEN("tk1_to_tk1") {
    orig = Circuit(1);
    orig.add_op<unsigned>(OpType::TK1, {0.2, 0.3, 0.4}, {0});
    res = CircPool::tk1_to_tk1(0.2, 0.3, 0.4);
  }
  GIVEN("CCX") {
    orig = Circuit(3);
    orig.add_op<unsigned>(OpType::CCX, {0, 1, 2});
    res = CircPool::CCX();
  }
  GIVEN("BRIDGE") {
    orig = Circuit(3);
    orig.add_op<unsigned>(OpType::BRIDGE, {0, 1, 2});
    res = CircPool::BRIDGE();
  }
  GIVEN("H_CZ_H") {
    orig = Circuit(2);
    orig.add_op<unsigned>(OpType::CX, {0, 1});
    res = CircPool::H_CZ_H();
  }

  auto u_orig = tket_sim::get_unitary(orig);
  auto u_res = tket_sim::get_unitary(res);
  REQUIRE(u_res.isApprox(u_orig));
}

}  // namespace test_CircPool
}  // namespace tket
