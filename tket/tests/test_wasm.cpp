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

#include <catch2/catch.hpp>

#include "Circuit/CircUtils.hpp"
#include "Circuit/Circuit.hpp"
#include "Converters/PhasePoly.hpp"
#include "Eigen/src/Core/Matrix.h"
#include "Ops/ClassicalOps.hpp"
#include "Simulation/CircuitSimulator.hpp"
#include "testutil.hpp"

namespace tket {
namespace test_Boxes {

SCENARIO("generating circ with wasm") {
  std::string wasm_file = "string/with/path/to/wasm/file";

  std::string wasm_func = "stringNameOfWASMFunc";

  GIVEN("wasmop creation") {
    WASMOp wop = WASMOp(1, wasm_func, wasm_file);

    REQUIRE(wop.get_n() == 1);
    REQUIRE(wop.get_func_name() == wasm_func);
    REQUIRE(wop.get_file_path() == wasm_file);
  }
  GIVEN("wasmop to json") {
    const std::shared_ptr<WASMOp> wop_ptr =
        std::make_shared<WASMOp>(1, wasm_func, wasm_file);

    nlohmann::json j = wop_ptr->serialize();

    auto wopj = WASMOp::deserialize(j);

    REQUIRE(wop_ptr->is_equal(*wopj));
  }
  GIVEN("add wasmop to circ") {
    Circuit u(2, 2);

    const std::shared_ptr<WASMOp> wop_ptr =
        std::make_shared<WASMOp>(1, wasm_func, wasm_file);
    u.add_op<unsigned>(wop_ptr, {0});

    const std::shared_ptr<WASMOp> wop_ptr_2 =
        std::make_shared<WASMOp>(2, wasm_func, wasm_file);

    u.add_op<unsigned>(wop_ptr_2, {0, 1});
  }
  GIVEN("compare wasmop") {
    WASMOp wop = WASMOp(1, wasm_func, wasm_file);

    WASMOp wop_2 = WASMOp(2, wasm_func, wasm_file);

    REQUIRE(!wop.is_equal(wop_2));
  }
  GIVEN("compare wasmop II") {
    WASMOp wop = WASMOp(1, wasm_func, wasm_file);

    WASMOp wop_2 = WASMOp(1, wasm_file, wasm_func);

    REQUIRE(!wop.is_equal(wop_2));
  }
  GIVEN("compare wasmop III") {
    WASMOp wop = WASMOp(1, wasm_func, wasm_file);

    WASMOp wop_2 = WASMOp(1, wasm_func, wasm_func);

    REQUIRE(!wop.is_equal(wop_2));
  }
  GIVEN("compare wasmop IV") {
    WASMOp wop = WASMOp(1, wasm_func, wasm_file);

    WASMOp wop_2 = WASMOp(1, wasm_func, wasm_file);

    REQUIRE(wop.is_equal(wop_2));
  }
  GIVEN("wasmop is_extern") {
    const std::shared_ptr<WASMOp> wop_ptr =
        std::make_shared<WASMOp>(1, wasm_func, wasm_file);

    REQUIRE(wop_ptr->is_extern());
  }
}

}  // namespace test_Boxes
}  // namespace tket
