// Copyright 2019-2024 Cambridge Quantum Computing
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
#include <catch2/generators/catch_generators.hpp>

#include "tket/Circuit/Circuit.hpp"
#include "tket/Transformations/BasicOptimisation.hpp"

namespace tket {
namespace test_RedundancyRemoval {

struct TestGate {
  OpType opType;
  std::optional<Expr> param;
  std::vector<unsigned> args;
  TestGate(OpType _opType, std::vector<unsigned> _args)
      : opType(_opType), args(_args) {}
  TestGate(
      OpType _opType, std::optional<Expr> _param, std::vector<unsigned> _args)
      : opType(_opType), param(_param), args(_args) {}
};

struct TestCase {
  std::string name;
  TestGate gate1;
  TestGate gate2;
  bool gatesShouldCancel;
};

SCENARIO("Transforms::remove_redundancies removes typical redundancies") {
  Circuit original_circuit(5);
  Circuit test_circuit(original_circuit);
  auto test_case = GENERATE(
      TestCase{
          "noops", TestGate{OpType::noop, {0}}, TestGate{OpType::noop, {1}},
          true},
      TestCase{"H-H", TestGate{OpType::H, {0}}, TestGate{OpType::H, {0}}, true},
      TestCase{
          "Hs on different qubits", TestGate{OpType::H, {0}},
          TestGate{OpType::H, {1}}, false},
      TestCase{"X-X", TestGate{OpType::X, {0}}, TestGate{OpType::X, {0}}, true},
      TestCase{"Y-Y", TestGate{OpType::Y, {0}}, TestGate{OpType::Y, {0}}, true},
      TestCase{"Z-Z", TestGate{OpType::Z, {0}}, TestGate{OpType::Z, {0}}, true},
      TestCase{
          "S-Sdg", TestGate{OpType::S, {0}}, TestGate{OpType::Sdg, {0}}, true},
      TestCase{
          "Sdg-S", TestGate{OpType::Sdg, {0}}, TestGate{OpType::S, {0}}, true},
      TestCase{
          "T-Tdg", TestGate{OpType::T, {0}}, TestGate{OpType::Tdg, {0}}, true},
      TestCase{
          "Tdg-T", TestGate{OpType::Tdg, {0}}, TestGate{OpType::T, {0}}, true},
      TestCase{
          "V-Vdg", TestGate{OpType::V, {0}}, TestGate{OpType::Vdg, {0}}, true},
      TestCase{
          "Vdg-V", TestGate{OpType::Vdg, {0}}, TestGate{OpType::V, {0}}, true},
      TestCase{
          "U1-U1*", TestGate{OpType::U1, 0.5, {0}},
          TestGate{OpType::U1, -0.5, {0}}, true},
      TestCase{
          "Rz-Rz*", TestGate{OpType::Rz, 0.5, {0}},
          TestGate{OpType::Rz, -0.5, {0}}, true},
      TestCase{
          "Rx-Rx*", TestGate{OpType::Rx, 0.5, {0}},
          TestGate{OpType::Rx, -0.5, {0}}, true},
      TestCase{
          "Ry-Ry*", TestGate{OpType::Ry, 0.5, {0}},
          TestGate{OpType::Ry, -0.5, {0}}, true},
      TestCase{
          "SWAPS with matching ports", TestGate{OpType::SWAP, {0, 1}},
          TestGate{OpType::SWAP, {0, 1}}, true},
      TestCase{
          "SWAPS with swapped ports", TestGate{OpType::SWAP, {0, 1}},
          TestGate{OpType::SWAP, {1, 0}}, true},
      TestCase{
          "Cancelling CHs (same port order)", TestGate{OpType::CH, {0, 1}},
          TestGate{OpType::CH, {0, 1}}, true},
      TestCase{
          "Non-cancelling CHs (swapped port order)",
          TestGate{OpType::CH, {0, 1}}, TestGate{OpType::CH, {1, 0}}, false},
      TestCase{
          "Cancelling CXs (same port order)", TestGate{OpType::CX, {0, 1}},
          TestGate{OpType::CX, {0, 1}}, true},
      TestCase{
          "Non-cancelling CXs (swapped port order)",
          TestGate{OpType::CX, {0, 1}}, TestGate{OpType::CX, {1, 0}}, false},
      TestCase{
          "Cancelling CYs (same port order)", TestGate{OpType::CY, {0, 1}},
          TestGate{OpType::CY, {0, 1}}, true},
      TestCase{
          "Non-cancelling CYs (swapped port order)",
          TestGate{OpType::CY, {0, 1}}, TestGate{OpType::CY, {1, 0}}, false},
      TestCase{
          "Cancelling CZs (same port order)", TestGate{OpType::CZ, {0, 1}},
          TestGate{OpType::CZ, {0, 1}}, true},
      TestCase{
          "Cancelling CZs (swapped port order)", TestGate{OpType::CZ, {0, 1}},
          TestGate{OpType::CZ, {1, 0}}, true},
      TestCase{
          "Cancelling XXPhase s (same port order)",
          TestGate{OpType::XXPhase, 0.5, {0, 1}},
          TestGate{OpType::XXPhase, -0.5, {0, 1}}, true},
      TestCase{
          "Cancelling XXPhases (swapped port order)",
          TestGate{OpType::XXPhase, 0.5, {0, 1}},
          TestGate{OpType::XXPhase, -0.5, {1, 0}}, true},
      TestCase{
          "Cancelling YYPhases (same port order)",
          TestGate{OpType::YYPhase, 0.5, {0, 1}},
          TestGate{OpType::YYPhase, -0.5, {0, 1}}, true},
      TestCase{
          "Cancelling YYPhases (swapped port order)",
          TestGate{OpType::YYPhase, 0.5, {0, 1}},
          TestGate{OpType::YYPhase, -0.5, {1, 0}}, true},
      TestCase{
          "Cancelling ZZPhases (same port order)",
          TestGate{OpType::ZZPhase, 0.5, {0, 1}},
          TestGate{OpType::ZZPhase, -0.5, {0, 1}}, true},
      TestCase{
          "Cancelling ZZPhases (swapped port order)",
          TestGate{OpType::ZZPhase, 0.5, {0, 1}},
          TestGate{OpType::ZZPhase, -0.5, {1, 0}}, true},
      TestCase{
          "Cancelling CV=CVdg (same port order)", TestGate{OpType::CV, {0, 1}},
          TestGate{OpType::CVdg, {0, 1}}, true},
      TestCase{
          "Cancelling CVdg=CV (same port order)",
          TestGate{OpType::CVdg, {0, 1}}, TestGate{OpType::CV, {0, 1}}, true},
      TestCase{
          "Non-cancelling CV=CVdg (swapped port order)",
          TestGate{OpType::CV, {0, 1}}, TestGate{OpType::CVdg, {1, 0}}, false},
      TestCase{
          "Cancelling CSX=CSxdg (same port order)",
          TestGate{OpType::CSX, {0, 1}}, TestGate{OpType::CSXdg, {0, 1}}, true},
      TestCase{
          "Cancelling CSXdg=CSX (same port order)",
          TestGate{OpType::CSXdg, {0, 1}}, TestGate{OpType::CSX, {0, 1}}, true},
      TestCase{
          "Non-cancelling CSX=CSXdg (swapped port order)",
          TestGate{OpType::CSX, {0, 1}}, TestGate{OpType::CSXdg, {1, 0}},
          false},
      TestCase{
          "Cancelling CSdg=CS (same port order)",
          TestGate{OpType::CSdg, {0, 1}}, TestGate{OpType::CS, {0, 1}}, true},
      TestCase{
          "Cancelling CSdg=CS (swapped port order)",
          TestGate{OpType::CSdg, {0, 1}}, TestGate{OpType::CS, {1, 0}}, true},
      TestCase{
          "Cancelling CCXs (same port order)", TestGate{OpType::CCX, {0, 1, 2}},
          TestGate{OpType::CCX, {0, 1, 2}}, true},
      TestCase{
          "Cancelling CCXs (control ports swapped)",
          TestGate{OpType::CCX, {1, 0, 2}}, TestGate{OpType::CCX, {0, 1, 2}},
          true},
      TestCase{
          "Non-cancelling CCXs (X port swapped)",
          TestGate{OpType::CCX, {0, 1, 2}}, TestGate{OpType::CCX, {0, 2, 1}},
          false},
      TestCase{
          "Cancelling CSWAPs (same port order)",
          TestGate{OpType::CSWAP, {0, 1, 2}},
          TestGate{OpType::CSWAP, {0, 1, 2}}, true},
      TestCase{
          "Cancelling CSWAPs (swap ports swapped)",
          TestGate{OpType::CSWAP, {0, 2, 1}},
          TestGate{OpType::CSWAP, {0, 1, 2}}, true},
      TestCase{
          "Non-cancelling CSWAPs (control port swapped)",
          TestGate{OpType::CSWAP, {0, 1, 2}},
          TestGate{OpType::CSWAP, {1, 0, 2}}, false},
      TestCase{
          "Cancelling ECRs (same port order)", TestGate{OpType::ECR, {0, 1}},
          TestGate{OpType::ECR, {0, 1}}, true},
      TestCase{
          "Non-cancelling ECRs (swapped port order)",
          TestGate{OpType::ECR, {1, 0}}, TestGate{OpType::ECR, {0, 1}}, false},
      TestCase{
          "Cancelling BRIDGEs (same port order)",
          TestGate{OpType::BRIDGE, {0, 1, 2}},
          TestGate{OpType::BRIDGE, {0, 1, 2}}, true},
      TestCase{
          "Non-cancelling BRIDGEs (swapped port order)",
          TestGate{OpType::BRIDGE, {0, 1, 2}},
          TestGate{OpType::BRIDGE, {0, 2, 1}}, false},
      TestCase{
          "Cancelling XXPhase3s (same port order)",
          TestGate{OpType::XXPhase3, 0.5, {0, 1, 2}},
          TestGate{OpType::XXPhase3, -0.5, {0, 1, 2}}, true},
      TestCase{
          "Cancelling XXPhase3s (swapped port order)",
          TestGate{OpType::XXPhase3, 0.5, {0, 1, 2}},
          TestGate{OpType::XXPhase3, -0.5, {0, 2, 1}}, true},
      TestCase{
          "Cancelling CnXs n=4 (same port order)",
          TestGate{OpType::CnX, {0, 1, 2, 3, 4}},
          TestGate{OpType::CnX, {0, 1, 2, 3, 4}}, true},
      TestCase{
          "Cancelling CnXs n=4 (swapped control ports)",
          TestGate{OpType::CnX, {0, 2, 1, 3, 4}},
          TestGate{OpType::CnX, {2, 3, 0, 1, 4}}, true},
      TestCase{
          "Cancelling CnXs n=4 (swapped X port)",
          TestGate{OpType::CnX, {0, 2, 1, 3, 4}},
          TestGate{OpType::CnX, {2, 4, 0, 1, 3}}, false},
      TestCase{
          "Cancelling CnYs n=4 (same port order)",
          TestGate{OpType::CnY, {0, 1, 2, 3, 4}},
          TestGate{OpType::CnY, {0, 1, 2, 3, 4}}, true},
      TestCase{
          "Cancelling CnYs n=4 (swapped control ports)",
          TestGate{OpType::CnY, {0, 2, 1, 3, 4}},
          TestGate{OpType::CnY, {3, 2, 0, 1, 4}}, true},
      TestCase{
          "Cancelling CnYs n=4 (swapped Y port)",
          TestGate{OpType::CnY, {0, 2, 1, 3, 4}},
          TestGate{OpType::CnY, {3, 2, 4, 1, 0}}, false},
      TestCase{
          "Cancelling CnZs n=4 (same port order)",
          TestGate{OpType::CnZ, {0, 1, 2, 3, 4}},
          TestGate{OpType::CnZ, {0, 1, 2, 3, 4}}, true},
      TestCase{
          "Cancelling CnZs n=4 (swapped ports, also Z port)",
          TestGate{OpType::CnZ, {0, 2, 1, 3, 4}},
          TestGate{OpType::CnZ, {2, 4, 0, 1, 3}}, true},
      TestCase{
          "Cancelling PhaseGadgets (same port order)",
          TestGate{OpType::PhaseGadget, 0.5, {0, 1, 2, 3, 4}},
          TestGate{OpType::PhaseGadget, -0.5, {0, 1, 2, 3, 4}}, true},
      TestCase{
          "Cancelling PhaseGadgets (ports permuted)",
          TestGate{OpType::PhaseGadget, 0.5, {0, 1, 2, 3, 4}},
          TestGate{OpType::PhaseGadget, -0.5, {2, 3, 4, 1, 0}}, true},
      TestCase{
          "Cancelling CU1s (same port order)",
          TestGate{OpType::CU1, 0.5, {0, 1}},
          TestGate{OpType::CU1, -0.5, {0, 1}}, true},
      TestCase{
          "Cancelling CU1s (swapped port order)",
          TestGate{OpType::CU1, 0.5, {0, 1}},
          TestGate{OpType::CU1, -0.5, {1, 0}}, true});
  const auto& test_name = test_case.name;
  const auto& gate1 = test_case.gate1;
  const auto& gate2 = test_case.gate2;
  GIVEN(test_name) {
    if (gate1.param) {
      test_circuit.add_op<unsigned>(
          gate1.opType, gate1.param.value(), gate1.args);
      test_circuit.add_op<unsigned>(
          gate2.opType, gate2.param.value(), gate2.args);
    } else {
      test_circuit.add_op<unsigned>(gate1.opType, gate1.args);
      test_circuit.add_op<unsigned>(gate2.opType, gate2.args);
    }
    Circuit untransformed_circuit(test_circuit);
    WHEN("calling Transforms::remove_redundancies on circuit") {
      auto circuit_has_changed =
          Transforms::remove_redundancies().apply(test_circuit);
      if (test_case.gatesShouldCancel) {
        THEN("gates should be removed") {
          CHECK(circuit_has_changed);
          REQUIRE(test_circuit.circuit_equality(original_circuit));
        }
      } else {
        THEN("circuit should not change") {
          CHECK_FALSE(circuit_has_changed);
          REQUIRE(test_circuit.circuit_equality(untransformed_circuit));
        }
      }
    }
  }
}

SCENARIO("Transforms::remove_redundancies removes nested redundancies") {
  Circuit original_circuit(5);
  Circuit test_circuit(original_circuit);
  GIVEN("A circuit with nested redundancies") {
    test_circuit.add_op<unsigned>(OpType::H, {0});
    test_circuit.add_op<unsigned>(OpType::CnZ, {0, 1, 2, 3, 4});
    test_circuit.add_op<unsigned>(OpType::CU3, {0.0, 0.4, 0.2}, {2, 3});
    test_circuit.add_op<unsigned>(OpType::noop, {4});
    test_circuit.add_op<unsigned>(OpType::ISWAP, 0.5, {1, 4});
    test_circuit.add_op<unsigned>(OpType::ISWAP, -0.5, {1, 4});
    test_circuit.add_op<unsigned>(OpType::ISWAP, 0.5, {2, 3});
    test_circuit.add_op<unsigned>(OpType::ISWAP, -0.5, {3, 2});
    test_circuit.add_op<unsigned>(OpType::CU3, {0.0, -0.2, -0.4}, {3, 2});
    test_circuit.add_op<unsigned>(OpType::CnX, {0, 1, 2, 3, 4});
    test_circuit.add_op<unsigned>(OpType::CnX, {0, 3, 2, 1, 4});
    test_circuit.add_op<unsigned>(OpType::CnZ, {0, 1, 4, 2, 3});
    test_circuit.add_op<unsigned>(OpType::H, {0});
    WHEN("calling Transforms::remove_redundancies on circuit") {
      Transforms::remove_redundancies().apply(test_circuit);
      THEN("all gates should be removed") {
        REQUIRE(test_circuit.circuit_equality(original_circuit));
      }
    }
  }
}

SCENARIO(
    "Transforms::remove_redundancies reduces gate depth if some gates are "
    "redundant") {
  Circuit original_circuit(3);
  Circuit test_circuit(original_circuit);
  GIVEN("A circuit with nested redundancies") {
    test_circuit.add_op<unsigned>(OpType::CZ, {0, 1});
    test_circuit.add_op<unsigned>(OpType::H, {0});
    test_circuit.add_op<unsigned>(OpType::CZ, {1, 2});
    test_circuit.add_op<unsigned>(OpType::H, {2});
    test_circuit.add_op<unsigned>(OpType::H, {2});
    test_circuit.add_op<unsigned>(OpType::CZ, {2, 1});
    test_circuit.add_op<unsigned>(OpType::CU3, {0.0, 0.4, 0.2}, {0, 1});
    test_circuit.add_op<unsigned>(OpType::CU3, {0.0, -0.2, -0.4}, {1, 0});
    test_circuit.add_op<unsigned>(OpType::CCX, {0, 1, 2});
    test_circuit.add_op<unsigned>(OpType::CnX, {1, 0, 2});
    test_circuit.add_op<unsigned>(OpType::H, {0});
    test_circuit.add_op<unsigned>(OpType::CY, {0, 2});
    WHEN("calling Transforms::remove_redundancies on circuit") {
      Transforms::remove_redundancies().apply(test_circuit);
      THEN("all gates should be removed") {
        REQUIRE(test_circuit.depth() <= 6);
      }
    }
  }
}

}  // namespace test_RedundancyRemoval
}  // namespace tket
