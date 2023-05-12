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

#include <Eigen/Core>
#include "testutil.hpp"
#include "tket/Circuit/CircUtils.hpp"
#include "tket/Circuit/Circuit.hpp"
#include "tket/Converters/PhasePoly.hpp"
#include "tket/Ops/ClassicalOps.hpp"
#include "tket/Simulation/CircuitSimulator.hpp"

namespace tket {
namespace test_Boxes {

SCENARIO("generating circ with wasm") {
  std::string wasm_file = "string/with/path/to/wasm/file";

  std::string wasm_func = "stringNameOfWASMFunc";

  std::vector<unsigned> uv = {2, 1};

  std::vector<unsigned> uv_2 = {1};

  std::vector<unsigned> uv_3 = {};

  GIVEN("wasmop creation") {
    WASMOp wop = WASMOp(4, 1, uv, uv_2, wasm_func, wasm_file);

    REQUIRE(wop.get_n_i32() == 3);
    REQUIRE(wop.get_func_name() == wasm_func);
    REQUIRE(wop.get_wasm_file_uid() == wasm_file);
  }
  GIVEN("wasmop to json") {
    const std::shared_ptr<WASMOp> wop_ptr =
        std::make_shared<WASMOp>(4, 1, uv, uv_2, wasm_func, wasm_file);

    nlohmann::json j = wop_ptr->serialize();

    auto wopj = WASMOp::deserialize(j);

    REQUIRE(wop_ptr->is_equal(*wopj));
  }
  GIVEN("add wasmop to circ") {
    Circuit u(6, 6);

    const std::shared_ptr<WASMOp> wop_ptr =
        std::make_shared<WASMOp>(1, 1, uv_2, uv_3, wasm_func, wasm_file);

    std::vector<Bit> bits = {Bit(0)};

    std::vector<UnitID> new_args;

    new_args.push_back(Bit(0));

    WasmState wuid = WasmState(0);

    new_args.push_back(wuid);

    u.add_op<UnitID>(wop_ptr, new_args);

    const std::shared_ptr<WASMOp> wop_ptr_2 =
        std::make_shared<WASMOp>(6, 1, uv, uv, wasm_func, wasm_file);

    u.add_op<UnitID>(
        wop_ptr_2,
        {Bit(0),          //
         Bit(1),          //
         Bit(2),          //
         Bit(3),          //
         Bit(4),          //
         Bit(5),          //
         WasmState(0)});  // needs 6 bits + WasmState
  }
  GIVEN("compare wasmop") {
    WASMOp wop = WASMOp(4, 1, uv, uv_2, wasm_func, wasm_file);

    WASMOp wop_2 = WASMOp(2, 1, uv_2, uv_2, wasm_func, wasm_file);

    REQUIRE(!wop.is_equal(wop_2));
  }
  GIVEN("compare wasmop II") {
    WASMOp wop = WASMOp(6, 1, uv, uv, wasm_func, wasm_file);

    WASMOp wop_2 = WASMOp(4, 1, uv, uv_2, wasm_func, wasm_file);

    REQUIRE(!wop.is_equal(wop_2));
  }
  GIVEN("compare wasmop III") {
    WASMOp wop = WASMOp(4, 1, uv, uv_2, wasm_func, wasm_file);

    WASMOp wop_2 = WASMOp(4, 1, uv, uv_2, wasm_func, wasm_func);

    REQUIRE(!wop.is_equal(wop_2));
  }
  GIVEN("compare wasmop IV") {
    WASMOp wop = WASMOp(4, 1, uv, uv_2, wasm_func, wasm_file);

    WASMOp wop_2 = WASMOp(4, 1, uv, uv_2, wasm_file, wasm_file);

    REQUIRE(!wop.is_equal(wop_2));
  }
  GIVEN("compare wasmop V") {
    WASMOp wop = WASMOp(4, 1, uv, uv_2, wasm_func, wasm_file);

    WASMOp wop_2 = WASMOp(4, 1, uv, uv_2, wasm_func, wasm_file);

    REQUIRE(wop.is_equal(wop_2));
  }
  GIVEN("wasmop is_extern") {
    const std::shared_ptr<WASMOp> wop_ptr =
        std::make_shared<WASMOp>(4, 1, uv, uv_2, wasm_func, wasm_file);

    REQUIRE(wop_ptr->is_extern());
  }
  GIVEN("wasmop add circuit") {
    const std::shared_ptr<WASMOp> wop_ptr =
        std::make_shared<WASMOp>(4, 1, uv, uv_2, wasm_func, wasm_file);

    REQUIRE(wop_ptr->is_extern());
  }
  GIVEN("wasmop add circuit II") {
    const std::shared_ptr<WASMOp> wop_ptr =
        std::make_shared<WASMOp>(4, 1, uv, uv_2, wasm_func, wasm_file);

    REQUIRE(wop_ptr->is_extern());
  }
  GIVEN("wasmop add circuit III") {
    Circuit u(1, 1);

    const std::shared_ptr<WASMOp> wop_ptr =
        std::make_shared<WASMOp>(2, 1, uv_2, uv_2, wasm_func, wasm_file);

    const std::shared_ptr<WASMOp> wop_ptr_2 =
        std::make_shared<WASMOp>(1, 1, uv_2, uv_3, wasm_func, wasm_file);

    u.add_op<UnitID>(wop_ptr, {Bit(0), Bit(0), WasmState(0)});
    u.add_op<UnitID>(wop_ptr_2, {Bit(0), WasmState(0)});

    u.assert_valid();

    REQUIRE(u.depth() == 2);
  }
  GIVEN("wasmop add circuit IV") {
    Circuit u(1, 1);

    const std::shared_ptr<WASMOp> wop_ptr =
        std::make_shared<WASMOp>(2, 1, uv_2, uv_2, wasm_func, wasm_file);

    const std::shared_ptr<WASMOp> wop_ptr_2 =
        std::make_shared<WASMOp>(1, 1, uv_2, uv_3, wasm_func, wasm_file);

    REQUIRE_THROWS(u.add_op<unsigned>(wop_ptr, {0}));
    REQUIRE_THROWS(u.add_op<unsigned>(wop_ptr_2, {0, 0, 0}));

    u.assert_valid();
    REQUIRE(u.depth() == 0);
  }
  GIVEN("wasmop add circuit V") {
    Circuit u(1, 1);

    const std::shared_ptr<WASMOp> wop_ptr =
        std::make_shared<WASMOp>(2, 1, uv_2, uv_2, wasm_func, wasm_file);

    const std::shared_ptr<WASMOp> wop_ptr_2 =
        std::make_shared<WASMOp>(1, 3, uv_2, uv_3, wasm_func, wasm_file);

    u.add_op<UnitID>(wop_ptr, {Bit(0), Bit(0), WasmState(0)});
    u.add_op<UnitID>(
        wop_ptr_2, {Bit(0), WasmState(0), WasmState(1), WasmState(2)});

    u.assert_valid();
    REQUIRE(u.depth() == 2);
    REQUIRE(u.w_inputs().size() == 3);
    REQUIRE(u.w_outputs().size() == 3);
  }
}

SCENARIO("test wasm uid") {
  GIVEN("wasm uid") { WasmState wuid = WasmState(); }
  GIVEN("wasm uid - compare") {
    WasmState wuid = WasmState();
    WasmState wuid_2 = WasmState();
    REQUIRE(wuid == wuid_2);
  }
  GIVEN("wasm uid - compare 2") {
    WasmState wuid = WasmState(1);
    WasmState wuid_2 = WasmState(3);
    WasmState wuid_3 = WasmState(3);

    REQUIRE(wuid != wuid_2);
    REQUIRE(wuid_3 == wuid_2);
    REQUIRE(wuid != wuid_3);
  }
  GIVEN("wasm uid - create bit from wasm") {
    WasmState wuid = WasmState();
    REQUIRE_THROWS(Bit(wuid));
  }
}

}  // namespace test_Boxes
}  // namespace tket
