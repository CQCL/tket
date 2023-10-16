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
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <tket/Circuit/Simulation/CircuitSimulator.hpp>
#include <tket/Gate/SymTable.hpp>

#include "../testutil.hpp"
#include "tket/Circuit/Circuit.hpp"
#include "tket/Circuit/PauliExpBoxes.hpp"

namespace tket {
namespace test_PauliExpBoxes {

SCENARIO("Pauli gadgets", "[boxes]") {
  GIVEN("Basis Circuit check") {
    PauliExpBox pbox(SymPauliTensor({Pauli::X}, 1.0));
    auto circ = pbox.to_circuit();
    circ->decompose_boxes_recursively();
    Circuit comp(1);
    comp.add_op<unsigned>(OpType::H, {0});
    comp.add_op<unsigned>(OpType::Rz, 1.0, {0});
    comp.add_op<unsigned>(OpType::H, {0});
    REQUIRE(*circ == comp);
  }
  GIVEN("Empty PauliExpBox compiles to empty circuit") {
    Circuit empty_circuit(0);
    PauliExpBox pbox;
    auto empty_pbox_circuit = pbox.to_circuit();
    REQUIRE(*empty_pbox_circuit == empty_circuit);
  }
  GIVEN("X") {
    // ---PauliExpBox([X], t)----Rx(-t)--- should be the identity
    double t = 1.687029013593215;
    Circuit c(1);
    DensePauliMap pauli_x = {Pauli::X};
    PauliExpBox pbox(SymPauliTensor(pauli_x, t));
    REQUIRE(pbox.get_paulis() == pauli_x);
    c.add_box(pbox, uvec{0});
    c.add_op<unsigned>(OpType::Rx, -t, {0});
    Eigen::Matrix2Xcd u = tket_sim::get_unitary(c);
    REQUIRE((u - Eigen::Matrix2cd::Identity()).cwiseAbs().sum() < ERR_EPS);
  }
  GIVEN("Y") {
    // ---PauliExpBox([Y], t)----Ry(-t)--- should be the identity
    double t = 1.6791969622440162;
    Circuit c(1);
    PauliExpBox pbox(SymPauliTensor({Pauli::Y}, t));
    c.add_box(pbox, uvec{0});
    c.add_op<unsigned>(OpType::Ry, -t, {0});
    Eigen::Matrix2Xcd u = tket_sim::get_unitary(c);
    REQUIRE((u - Eigen::Matrix2cd::Identity()).cwiseAbs().sum() < ERR_EPS);
  }
  GIVEN("Z") {
    // ---PauliExpBox([Z], t)----Rz(-t)--- should be the identity
    double t = 1.7811410013115163;
    Circuit c(1);
    PauliExpBox pbox(SymPauliTensor({Pauli::Z}, t));
    c.add_box(pbox, uvec{0});
    c.add_op<unsigned>(OpType::Rz, -t, {0});
    Eigen::Matrix2Xcd u = tket_sim::get_unitary(c);
    REQUIRE((u - Eigen::Matrix2cd::Identity()).cwiseAbs().sum() < ERR_EPS);
  }
  GIVEN("II") {
    double t = 0.10154905537993009;
    Eigen::Matrix4cd a;
    a << 1., 0., 0., 0., 0., 1., 0., 0., 0., 0., 1., 0., 0., 0., 0., 1.;
    ExpBox ebox(a, +0.5 * PI * t);
    Circuit c(2);
    c.add_box(ebox, {0, 1});
    PauliExpBox pbox(SymPauliTensor({Pauli::I, Pauli::I}, t));
    c.add_box(pbox, {0, 1});
    Eigen::MatrixXcd u = tket_sim::get_unitary(c);
    REQUIRE((u - Eigen::Matrix4cd::Identity()).cwiseAbs().sum() < ERR_EPS);
  }
  GIVEN("IX") {
    double t = -0.9124813027056411;
    Eigen::Matrix4cd a;
    a << 0., 1., 0., 0., 1., 0., 0., 0., 0., 0., 0., 1., 0., 0., 1., 0.;
    ExpBox ebox(a, +0.5 * PI * t);
    Circuit c(2);
    c.add_box(ebox, {0, 1});
    PauliExpBox pbox(SymPauliTensor({Pauli::I, Pauli::X}, t));
    c.add_box(pbox, {0, 1});
    Eigen::MatrixXcd u = tket_sim::get_unitary(c);
    REQUIRE((u - Eigen::Matrix4cd::Identity()).cwiseAbs().sum() < ERR_EPS);
  }
  GIVEN("IY") {
    double t = 0.4906808577976969;
    Eigen::Matrix4cd a;
    a << 0., -i_, 0., 0., i_, 0., 0., 0., 0., 0., 0., -i_, 0., 0., i_, 0.;
    ExpBox ebox(a, +0.5 * PI * t);
    Circuit c(2);
    c.add_box(ebox, {0, 1});
    PauliExpBox pbox(SymPauliTensor({Pauli::I, Pauli::Y}, t));
    c.add_box(pbox, {0, 1});
    Eigen::MatrixXcd u = tket_sim::get_unitary(c);
    REQUIRE((u - Eigen::Matrix4cd::Identity()).cwiseAbs().sum() < ERR_EPS);
  }
  GIVEN("IZ") {
    double t = -0.9536579982905538;
    Eigen::Matrix4cd a;
    a << 1., 0., 0., 0., 0., -1., 0., 0., 0., 0., 1., 0., 0., 0., 0., -1.;
    ExpBox ebox(a, +0.5 * PI * t);
    Circuit c(2);
    c.add_box(ebox, {0, 1});
    PauliExpBox pbox(SymPauliTensor({Pauli::I, Pauli::Z}, t));
    c.add_box(pbox, {0, 1});
    Eigen::MatrixXcd u = tket_sim::get_unitary(c);
    REQUIRE((u - Eigen::Matrix4cd::Identity()).cwiseAbs().sum() < ERR_EPS);
  }
  GIVEN("XI") {
    double t = 0.9735728239081902;
    Eigen::Matrix4cd a;
    a << 0., 0., 1., 0., 0., 0., 0., 1., 1., 0., 0., 0., 0., 1., 0., 0.;
    ExpBox ebox(a, +0.5 * PI * t);
    Circuit c(2);
    c.add_box(ebox, {0, 1});
    PauliExpBox pbox(SymPauliTensor({Pauli::X, Pauli::I}, t));
    c.add_box(pbox, {0, 1});
    Eigen::MatrixXcd u = tket_sim::get_unitary(c);
    REQUIRE((u - Eigen::Matrix4cd::Identity()).cwiseAbs().sum() < ERR_EPS);
  }
  GIVEN("XX") {
    double t = 0.27251750245844586;
    Eigen::Matrix4cd a;
    a << 0., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0., 0.;
    ExpBox ebox(a, +0.5 * PI * t);
    Circuit c(2);
    c.add_box(ebox, {0, 1});
    PauliExpBox pbox(SymPauliTensor({Pauli::X, Pauli::X}, t));
    c.add_box(pbox, {0, 1});
    Eigen::MatrixXcd u = tket_sim::get_unitary(c);
    REQUIRE((u - Eigen::Matrix4cd::Identity()).cwiseAbs().sum() < ERR_EPS);
  }
  GIVEN("XY") {
    double t = -0.7252139115522431;
    Eigen::Matrix4cd a;
    a << 0., 0., 0., -i_, 0., 0., i_, 0., 0., -i_, 0., 0., i_, 0., 0., 0.;
    ExpBox ebox(a, +0.5 * PI * t);
    Circuit c(2);
    c.add_box(ebox, {0, 1});
    PauliExpBox pbox(SymPauliTensor({Pauli::X, Pauli::Y}, t));
    c.add_box(pbox, {0, 1});
    Eigen::MatrixXcd u = tket_sim::get_unitary(c);
    REQUIRE((u - Eigen::Matrix4cd::Identity()).cwiseAbs().sum() < ERR_EPS);
  }
  GIVEN("XZ") {
    double t = 0.7474044702065266;
    Eigen::Matrix4cd a;
    a << 0., 0., 1., 0., 0., 0., 0., -1., 1., 0., 0., 0., 0., -1., 0., 0.;
    ExpBox ebox(a, +0.5 * PI * t);
    Circuit c(2);
    c.add_box(ebox, {0, 1});
    PauliExpBox pbox(SymPauliTensor({Pauli::X, Pauli::Z}, t));
    c.add_box(pbox, {0, 1});
    Eigen::MatrixXcd u = tket_sim::get_unitary(c);
    REQUIRE((u - Eigen::Matrix4cd::Identity()).cwiseAbs().sum() < ERR_EPS);
  }
  GIVEN("YI") {
    double t = 0.31314409051199577;
    Eigen::Matrix4cd a;
    a << 0., 0., -i_, 0., 0., 0., 0., -i_, i_, 0., 0., 0., 0., i_, 0., 0.;
    ExpBox ebox(a, +0.5 * PI * t);
    Circuit c(2);
    c.add_box(ebox, {0, 1});
    PauliExpBox pbox(SymPauliTensor({Pauli::Y, Pauli::I}, t));
    c.add_box(pbox, {0, 1});
    Eigen::MatrixXcd u = tket_sim::get_unitary(c);
    REQUIRE((u - Eigen::Matrix4cd::Identity()).cwiseAbs().sum() < ERR_EPS);
  }
  GIVEN("YX") {
    double t = -0.4855765841278301;
    Eigen::Matrix4cd a;
    a << 0., 0., 0., -i_, 0., 0., -i_, 0., 0., i_, 0., 0., i_, 0., 0., 0.;
    ExpBox ebox(a, +0.5 * PI * t);
    Circuit c(2);
    c.add_box(ebox, {0, 1});
    PauliExpBox pbox(SymPauliTensor({Pauli::Y, Pauli::X}, t));
    c.add_box(pbox, {0, 1});
    Eigen::MatrixXcd u = tket_sim::get_unitary(c);
    REQUIRE((u - Eigen::Matrix4cd::Identity()).cwiseAbs().sum() < ERR_EPS);
  }
  GIVEN("YY") {
    double t = 0.3103588880238326;
    Eigen::Matrix4cd a;
    a << 0., 0., 0., -1., 0., 0., 1., 0., 0., 1., 0., 0., -1., 0., 0., 0.;
    ExpBox ebox(a, +0.5 * PI * t);
    Circuit c(2);
    c.add_box(ebox, {0, 1});
    PauliExpBox pbox(SymPauliTensor({Pauli::Y, Pauli::Y}, t));
    c.add_box(pbox, {0, 1});
    Eigen::MatrixXcd u = tket_sim::get_unitary(c);
    REQUIRE((u - Eigen::Matrix4cd::Identity()).cwiseAbs().sum() < ERR_EPS);
  }
  GIVEN("YZ") {
    double t = -0.1130806991828821;
    Eigen::Matrix4cd a;
    a << 0., 0., -i_, 0., 0., 0., 0., i_, i_, 0., 0., 0., 0., -i_, 0., 0.;
    ExpBox ebox(a, +0.5 * PI * t);
    Circuit c(2);
    c.add_box(ebox, {0, 1});
    PauliExpBox pbox(SymPauliTensor({Pauli::Y, Pauli::Z}, t));
    c.add_box(pbox, {0, 1});
    Eigen::MatrixXcd u = tket_sim::get_unitary(c);
    REQUIRE((u - Eigen::Matrix4cd::Identity()).cwiseAbs().sum() < ERR_EPS);
  }
  GIVEN("ZI") {
    double t = -0.21235736398463878;
    Eigen::Matrix4cd a;
    a << 1., 0., 0., 0., 0., 1., 0., 0., 0., 0., -1., 0., 0., 0., 0., -1.;
    ExpBox ebox(a, +0.5 * PI * t);
    Circuit c(2);
    c.add_box(ebox, {0, 1});
    PauliExpBox pbox(SymPauliTensor({Pauli::Z, Pauli::I}, t));
    c.add_box(pbox, {0, 1});
    Eigen::MatrixXcd u = tket_sim::get_unitary(c);
    REQUIRE((u - Eigen::Matrix4cd::Identity()).cwiseAbs().sum() < ERR_EPS);
  }
  GIVEN("ZX") {
    double t = 0.5841730428035412;
    Eigen::Matrix4cd a;
    a << 0., 1., 0., 0., 1., 0., 0., 0., 0., 0., 0., -1., 0., 0., -1., 0.;
    ExpBox ebox(a, +0.5 * PI * t);
    Circuit c(2);
    c.add_box(ebox, {0, 1});
    PauliExpBox pbox(SymPauliTensor({Pauli::Z, Pauli::X}, t));
    c.add_box(pbox, {0, 1});
    Eigen::MatrixXcd u = tket_sim::get_unitary(c);
    REQUIRE((u - Eigen::Matrix4cd::Identity()).cwiseAbs().sum() < ERR_EPS);
  }
  GIVEN("ZY") {
    double t = 0.4300676558283072;
    Eigen::Matrix4cd a;
    a << 0., -i_, 0., 0., i_, 0., 0., 0., 0., 0., 0., i_, 0., 0., -i_, 0.;
    ExpBox ebox(a, +0.5 * PI * t);
    Circuit c(2);
    c.add_box(ebox, {0, 1});
    PauliExpBox pbox(SymPauliTensor({Pauli::Z, Pauli::Y}, t));
    c.add_box(pbox, {0, 1});
    Eigen::MatrixXcd u = tket_sim::get_unitary(c);
    REQUIRE((u - Eigen::Matrix4cd::Identity()).cwiseAbs().sum() < ERR_EPS);
  }
  GIVEN("ZZ") {
    double t = -0.18497547540553927;
    Eigen::Matrix4cd a;
    a << 1., 0., 0., 0., 0., -1., 0., 0., 0., 0., -1., 0., 0., 0., 0., 1.;
    ExpBox ebox(a, +0.5 * PI * t);
    Circuit c(2);
    c.add_box(ebox, {0, 1});
    PauliExpBox pbox(SymPauliTensor({Pauli::Z, Pauli::Z}, t));
    c.add_box(pbox, {0, 1});
    Eigen::MatrixXcd u = tket_sim::get_unitary(c);
    REQUIRE((u - Eigen::Matrix4cd::Identity()).cwiseAbs().sum() < ERR_EPS);
  }
  GIVEN("complex coefficient") {
    Expr ei{SymEngine::I};
    PauliExpBox pebox(SymPauliTensor({Pauli::Z}, ei));
    Expr p = pebox.get_phase();
    REQUIRE(p == ei);
  }
}
SCENARIO("Pauli gadget pairs", "[boxes]") {
  GIVEN("Basis Circuit check") {
    PauliExpPairBox pbox(SymPauliTensor({Pauli::X}, 1.0), SymPauliTensor({Pauli::I}, 0.0));
    auto circ = pbox.to_circuit();
    circ->decompose_boxes_recursively();
    Circuit comp(1);
    comp.add_op<unsigned>(OpType::H, {0});
    comp.add_op<unsigned>(OpType::Rz, 1.0, {0});
    comp.add_op<unsigned>(OpType::H, {0});
    REQUIRE(*circ == comp);
  }
  GIVEN("Empty PauliExpPairBox compiles to empty circuit") {
    Circuit empty_circuit(0);
    PauliExpPairBox pbox;
    auto empty_pbox_circuit = pbox.to_circuit();
    empty_pbox_circuit->decompose_boxes_recursively();
    REQUIRE(*empty_pbox_circuit == empty_circuit);
  }
  GIVEN("Construction with two pauli strings of different length throws") {
    DensePauliMap pauli_string0{Pauli::X, Pauli::Z};
    DensePauliMap pauli_string1{Pauli::X, Pauli::Z, Pauli::I};
    REQUIRE_THROWS_AS(
        PauliExpPairBox(
            SymPauliTensor(pauli_string0, 1.0),
            SymPauliTensor(pauli_string1, 1.0)),
        PauliExpBoxInvalidity);
  }
  GIVEN("is_clifford test cases") {
    SECTION("Empty Paulis") {
      REQUIRE(PauliExpPairBox(SymPauliTensor({}, 1.2), SymPauliTensor({}, 0.1))
                  .is_clifford());
    }
    SECTION("Various phases") {
      auto phase_case = GENERATE(
          // (phase0, phase1, expected is_clifford result)
          std::make_tuple(0.0, 0.0, true), std::make_tuple(0.5, 0.0, true),
          std::make_tuple(1.0, 0.0, true), std::make_tuple(1.5, 0.0, true),
          std::make_tuple(2.0, 0.0, true), std::make_tuple(0.5, 0.5, true),
          std::make_tuple(0.5, 1.0, true), std::make_tuple(0.5, 1.5, true),
          std::make_tuple(0.5, 2.0, true), std::make_tuple(0.0, 0.3, false),
          std::make_tuple(0.1, 0.3, false), std::make_tuple(1.1, 2.0, false));
      auto pbox = PauliExpPairBox(
          SymPauliTensor({Pauli::X, Pauli::Y}, get<0>(phase_case)),
          SymPauliTensor({Pauli::X, Pauli::Y}, get<1>(phase_case)));
      REQUIRE(pbox.is_clifford() == get<2>(phase_case));
    }
  }
  GIVEN("free_symbols") {
    auto a = SymTable::fresh_symbol("a");
    auto b = SymTable::fresh_symbol("b");
    auto ea = Expr(a);
    auto eb = Expr(b);
    DensePauliMap paulis0{Pauli::X};
    DensePauliMap paulis1{Pauli::Z};
    REQUIRE(PauliExpPairBox(
                SymPauliTensor(paulis0, 0.2), SymPauliTensor(paulis1, 0.4))
                .free_symbols()
                .empty());
    REQUIRE(
        PauliExpPairBox(
            SymPauliTensor(paulis0, ea), SymPauliTensor(paulis1, 0.4))
            .free_symbols() == SymSet{a});
    REQUIRE(
        PauliExpPairBox(
            SymPauliTensor(paulis0, 1.0), SymPauliTensor(paulis1, eb))
            .free_symbols() == SymSet{b});
    REQUIRE(
        PauliExpPairBox(
            SymPauliTensor(paulis0, ea), SymPauliTensor(paulis1, eb))
            .free_symbols() == SymSet{a, b});
  }
  GIVEN("dagger") {
    auto ea = Expr(SymTable::fresh_symbol("a"));
    DensePauliMap paulis0{Pauli::X};
    DensePauliMap paulis1{Pauli::Z};
    auto cx_config = CXConfigType::Star;
    auto box = PauliExpPairBox(
        SymPauliTensor(paulis0, ea), SymPauliTensor(paulis1, 0.4), cx_config);
    auto dagger_box =
        std::dynamic_pointer_cast<const PauliExpPairBox>(box.dagger());

    const auto [actual_paulis0, actual_paulis1] = dagger_box->get_paulis_pair();
    const auto [actual_phase0, actual_phase1] = dagger_box->get_phase_pair();
    REQUIRE(actual_paulis0 == paulis1);
    REQUIRE(actual_phase0 == -Expr(0.4));
    REQUIRE(actual_paulis1 == paulis0);
    REQUIRE(actual_phase1 == -ea);
    REQUIRE(dagger_box->get_cx_config() == cx_config);
  }
  GIVEN("transpose") {
    auto ea = Expr(SymTable::fresh_symbol("a"));
    DensePauliMap paulis0{Pauli::X, Pauli::Y, Pauli::Z, Pauli::I, Pauli::Y};
    DensePauliMap paulis1{Pauli::Y, Pauli::Y, Pauli::Z, Pauli::I, Pauli::Y};
    auto cx_config = CXConfigType::MultiQGate;

    WHEN("paulis1 contains odd number of Y") {
      auto box = PauliExpPairBox(
          SymPauliTensor(paulis0, ea), SymPauliTensor(paulis1, 0.4), cx_config);
      auto transpose_box =
          std::dynamic_pointer_cast<const PauliExpPairBox>(box.transpose());
      const auto [actual_paulis0, actual_paulis1] =
          transpose_box->get_paulis_pair();
      const auto [actual_phase0, actual_phase1] =
          transpose_box->get_phase_pair();
      REQUIRE(actual_paulis0 == paulis1);
      REQUIRE(actual_phase0 == -Expr(0.4));
      REQUIRE(actual_paulis1 == paulis0);
      REQUIRE(actual_phase1 == ea);
      REQUIRE(transpose_box->get_cx_config() == cx_config);
    }

    std::swap(paulis0, paulis1);
    WHEN("paulis0 contains odd number of Y") {
      auto box = PauliExpPairBox(
          SymPauliTensor(paulis0, ea), SymPauliTensor(paulis1, 0.4), cx_config);
      auto transpose_box =
          std::dynamic_pointer_cast<const PauliExpPairBox>(box.transpose());
      const auto [actual_paulis0, actual_paulis1] =
          transpose_box->get_paulis_pair();
      const auto [actual_phase0, actual_phase1] =
          transpose_box->get_phase_pair();
      REQUIRE(actual_paulis0 == paulis1);
      REQUIRE(actual_phase0 == Expr(0.4));
      REQUIRE(actual_paulis1 == paulis0);
      REQUIRE(actual_phase1 == -ea);
      REQUIRE(transpose_box->get_cx_config() == cx_config);
    }
  }
  GIVEN("symbol_substitution") {
    auto a = SymTable::fresh_symbol("a");
    auto b = SymTable::fresh_symbol("b");
    auto ea = Expr(a);
    auto eb = Expr(b);
    DensePauliMap paulis0{Pauli::X};
    DensePauliMap paulis1{Pauli::Z};

    auto box = PauliExpPairBox(
        SymPauliTensor(paulis0, ea), SymPauliTensor(paulis1, eb));

    WHEN("only first phase is substituted") {
      SymEngine::map_basic_basic sub_map{std::make_pair(a, Expr(0.8))};
      auto sub_box = std::dynamic_pointer_cast<const PauliExpPairBox>(
          box.symbol_substitution(sub_map));
      const auto [actual_phase0, actual_phase1] = sub_box->get_phase_pair();
      REQUIRE(actual_phase0 == Expr(0.8));
      REQUIRE(actual_phase1 == eb);
    }

    WHEN("only second phase is substituted") {
      SymEngine::map_basic_basic sub_map{std::make_pair(b, Expr(0.3))};
      auto sub_box = std::dynamic_pointer_cast<const PauliExpPairBox>(
          box.symbol_substitution(sub_map));
      const auto [actual_phase0, actual_phase1] = sub_box->get_phase_pair();
      REQUIRE(actual_phase0 == ea);
      REQUIRE(actual_phase1 == Expr(0.3));
    }

    WHEN("both phases are substituted") {
      SymEngine::map_basic_basic sub_map{
          std::make_pair(a, Expr(0.8)), std::make_pair(b, Expr(0.3))};
      auto sub_box = std::dynamic_pointer_cast<const PauliExpPairBox>(
          box.symbol_substitution(sub_map));
      const auto [actual_phase0, actual_phase1] = sub_box->get_phase_pair();
      REQUIRE(actual_phase0 == Expr(0.8));
      REQUIRE(actual_phase1 == Expr(0.3));
    }
  }
}

SCENARIO("Pauli gadget commuting sets", "[boxes]") {
  GIVEN("Basis Circuit check") {
    PauliExpCommutingSetBox pbox(
        {{{Pauli::X}, 1.0}, {{Pauli::I}, 0.0}, {{Pauli::I}, 0.0}});
    auto circ = pbox.to_circuit();
    circ->decompose_boxes_recursively();
    Circuit comp(1);
    comp.add_op<unsigned>(OpType::H, {0});
    comp.add_op<unsigned>(OpType::Rz, 1.0, {0});
    comp.add_op<unsigned>(OpType::H, {0});
    REQUIRE(*circ == comp);
  }
  GIVEN("Empty PauliExpPairBox compiles to empty circuit") {
    Circuit empty_circuit(0);
    PauliExpCommutingSetBox pbox;
    auto empty_pbox_circuit = pbox.to_circuit();
    empty_pbox_circuit->decompose_boxes_recursively();
    REQUIRE(*empty_pbox_circuit == empty_circuit);
  }
  GIVEN("Construction with no gadgets throws") {
    REQUIRE_THROWS_AS(
        PauliExpCommutingSetBox(std::vector<SymPauliTensor>{}),
        PauliExpBoxInvalidity);
  }
  GIVEN("Construction with pauli strings of different length throws") {
    DensePauliMap pauli_string0{Pauli::X, Pauli::Z};
    DensePauliMap pauli_string1{Pauli::X, Pauli::I};
    DensePauliMap pauli_string2{Pauli::X, Pauli::Z, Pauli::I};
    REQUIRE_THROWS_AS(
        PauliExpCommutingSetBox({
            SymPauliTensor(pauli_string0, 1.0),
            SymPauliTensor(pauli_string1, 1.0),
            SymPauliTensor(pauli_string2, 1.0),
        }),
        PauliExpBoxInvalidity);
  }
  GIVEN("Construction with non-commuting pauli strings throws") {
    DensePauliMap pauli_string0{Pauli::X, Pauli::Z};
    DensePauliMap pauli_string1{Pauli::Z, Pauli::I};
    DensePauliMap pauli_string2{Pauli::X, Pauli::Z};
    REQUIRE_THROWS_AS(
        PauliExpCommutingSetBox({
            {pauli_string0, 1.0},
            {pauli_string1, 1.0},
            {pauli_string2, 1.0},
        }),
        PauliExpBoxInvalidity);
  }
  GIVEN("is_clifford test cases") {
    SECTION("Empty Paulis") {
      REQUIRE(PauliExpCommutingSetBox({{{}, 1.2}, {{}, 0.1}, {{}, 1.1}})
                  .is_clifford());
    }
    SECTION("Various phases") {
      auto phase_case = GENERATE(
          // (phase0, phase1, expected is_clifford result)
          std::make_tuple(0.0, 0.0, 1.0, true),
          std::make_tuple(0.5, 0.0, 0.0, true),
          std::make_tuple(1.0, 0.0, 2.0, true),
          std::make_tuple(1.5, 0.0, 0.0, true),
          std::make_tuple(2.0, 0.0, 0.5, true),
          std::make_tuple(0.5, 0.5, 0.5, true),
          std::make_tuple(0.5, 1.0, 1.0, true),
          std::make_tuple(0.5, 1.5, 1.5, true),
          std::make_tuple(0.5, 2.0, 2.0, true),
          std::make_tuple(0.0, 0.3, 0.3, false),
          std::make_tuple(0.1, 0.3, 0.3, false),
          std::make_tuple(0.0, 0.0, 0.3, false),
          std::make_tuple(0.1, 0.3, 0.3, false),
          std::make_tuple(0.0, 2.0, 1.1, false),
          std::make_tuple(0.1, 0.3, 0.3, false),
          std::make_tuple(1.1, 2.0, 2.0, false));
      auto pbox = PauliExpCommutingSetBox({
          {{Pauli::I, Pauli::Y, Pauli::I}, get<0>(phase_case)},
          {{Pauli::X, Pauli::Y, Pauli::Z}, get<1>(phase_case)},
          {{Pauli::X, Pauli::Y, Pauli::Z}, get<2>(phase_case)},
      });
      REQUIRE(pbox.is_clifford() == get<3>(phase_case));
    }
  }
  GIVEN("free_symbols") {
    auto a = SymTable::fresh_symbol("a");
    auto b = SymTable::fresh_symbol("b");
    auto c = SymTable::fresh_symbol("c");
    auto ea = Expr(a);
    auto eb = Expr(b);
    auto ec = Expr(c);
    auto paulis0 = std::vector<Pauli>{Pauli::X};
    auto paulis1 = std::vector<Pauli>{Pauli::X};
    auto paulis2 = std::vector<Pauli>{Pauli::I};
    REQUIRE(PauliExpCommutingSetBox(
                {{paulis0, 0.2}, {paulis1, 0.4}, {paulis2, 0.3}})
                .free_symbols()
                .empty());
    REQUIRE(
        PauliExpCommutingSetBox({{paulis0, ea}, {paulis1, 0.4}, {paulis2, 0.3}})
            .free_symbols() == SymSet{a});
    REQUIRE(
        PauliExpCommutingSetBox({{paulis0, 0.2}, {paulis1, eb}, {paulis2, 0.3}})
            .free_symbols() == SymSet{b});
    REQUIRE(
        PauliExpCommutingSetBox({{paulis0, 0.2}, {paulis1, 0.4}, {paulis2, ec}})
            .free_symbols() == SymSet{c});
    REQUIRE(
        PauliExpCommutingSetBox({{paulis0, ea}, {paulis1, eb}, {paulis2, 0.3}})
            .free_symbols() == SymSet{a, b});
    REQUIRE(
        PauliExpCommutingSetBox({{paulis0, 0.2}, {paulis1, eb}, {paulis2, ec}})
            .free_symbols() == SymSet{b, c});
    REQUIRE(
        PauliExpCommutingSetBox({{paulis0, ea}, {paulis1, 0.4}, {paulis2, ec}})
            .free_symbols() == SymSet{a, c});
    REQUIRE(
        PauliExpCommutingSetBox({{paulis0, ea}, {paulis1, eb}, {paulis2, ec}})
            .free_symbols() == SymSet{a, b, c});
  }
  GIVEN("dagger") {
    auto ea = Expr(SymTable::fresh_symbol("a"));
    auto paulis0 = std::vector<Pauli>{Pauli::Z};
    auto paulis1 = std::vector<Pauli>{Pauli::I};
    auto paulis2 = std::vector<Pauli>{Pauli::Z};
    auto phase0 = ea;
    auto phase1 = Expr(0.4);
    auto phase2 = Expr(1.3);
    auto cx_config = CXConfigType::Tree;
    auto box = PauliExpCommutingSetBox(
        {{paulis0, phase0}, {paulis1, phase1}, {paulis2, phase2}}, cx_config);
    auto dagger_box =
        std::dynamic_pointer_cast<const PauliExpCommutingSetBox>(box.dagger());
    auto pauli_gadgets = dagger_box->get_pauli_gadgets();

    REQUIRE(pauli_gadgets.size() == 3);
    REQUIRE(pauli_gadgets[0].string == paulis0);
    REQUIRE(pauli_gadgets[0].coeff == -phase0);
    REQUIRE(pauli_gadgets[1].string == paulis1);
    REQUIRE(pauli_gadgets[1].coeff == -phase1);
    REQUIRE(pauli_gadgets[2].string == paulis2);
    REQUIRE(pauli_gadgets[2].coeff == -phase2);
    REQUIRE(dagger_box->get_cx_config() == cx_config);
  }
  GIVEN("transpose") {
    auto ea = Expr(SymTable::fresh_symbol("a"));
    auto paulis0 = std::vector<Pauli>{Pauli::Y, Pauli::Y, Pauli::Y, Pauli::Y};
    auto paulis1 = std::vector<Pauli>{Pauli::I, Pauli::Y, Pauli::Y, Pauli::Y};
    auto paulis2 = std::vector<Pauli>{Pauli::Y, Pauli::Y, Pauli::I, Pauli::I};
    auto paulis3 = std::vector<Pauli>{Pauli::Y, Pauli::I, Pauli::I, Pauli::I};
    auto phase0 = ea;
    auto phase1 = Expr(0.4);
    auto phase2 = Expr(1.3);
    auto phase3 = ea;
    auto cx_config = CXConfigType::Snake;
    auto box = PauliExpCommutingSetBox(
        {{paulis0, phase0},
         {paulis1, phase1},
         {paulis2, phase2},
         {paulis3, phase3}},
        cx_config);
    auto transpose_box =
        std::dynamic_pointer_cast<const PauliExpCommutingSetBox>(
            box.transpose());
    auto pauli_gadgets = transpose_box->get_pauli_gadgets();

    REQUIRE(pauli_gadgets.size() == 4);
    REQUIRE(pauli_gadgets[0].string == paulis0);
    REQUIRE(pauli_gadgets[0].coeff == phase0);
    REQUIRE(pauli_gadgets[1].string == paulis1);
    REQUIRE(pauli_gadgets[1].coeff == -phase1);
    REQUIRE(pauli_gadgets[2].string == paulis2);
    REQUIRE(pauli_gadgets[2].coeff == phase2);
    REQUIRE(pauli_gadgets[3].string == paulis3);
    REQUIRE(pauli_gadgets[3].coeff == -phase3);
    REQUIRE(transpose_box->get_cx_config() == cx_config);
  }
  GIVEN("symbol_substitution") {
    auto a = SymTable::fresh_symbol("a");
    auto b = SymTable::fresh_symbol("b");
    auto c = SymTable::fresh_symbol("c");
    auto ea = Expr(a);
    auto eb = Expr(b);
    auto ec = Expr(c);
    auto sub_a = Expr(0.8);
    auto sub_b = Expr(0.3);
    auto sub_c = Expr(2.3);
    auto paulis0 = std::vector<Pauli>{Pauli::X};
    auto paulis1 = std::vector<Pauli>{Pauli::X};
    auto paulis2 = std::vector<Pauli>{Pauli::X};

    auto box =
        PauliExpCommutingSetBox({{paulis0, ea}, {paulis1, eb}, {paulis2, ec}});

    SymEngine::map_basic_basic sub_map1{std::make_pair(a, sub_a)};
    auto sub_box1 = std::dynamic_pointer_cast<const PauliExpCommutingSetBox>(
        box.symbol_substitution(sub_map1));
    auto pauli_gadgets1 = sub_box1->get_pauli_gadgets();
    REQUIRE(pauli_gadgets1[0].coeff == sub_a);
    REQUIRE(pauli_gadgets1[1].coeff == eb);
    REQUIRE(pauli_gadgets1[2].coeff == ec);

    SymEngine::map_basic_basic sub_map2{std::make_pair(b, sub_b)};
    auto sub_box2 = std::dynamic_pointer_cast<const PauliExpCommutingSetBox>(
        box.symbol_substitution(sub_map2));
    auto pauli_gadgets2 = sub_box2->get_pauli_gadgets();
    REQUIRE(pauli_gadgets2[0].coeff == ea);
    REQUIRE(pauli_gadgets2[1].coeff == sub_b);
    REQUIRE(pauli_gadgets2[2].coeff == ec);

    SymEngine::map_basic_basic sub_map3{std::make_pair(c, sub_c)};
    auto sub_box3 = std::dynamic_pointer_cast<const PauliExpCommutingSetBox>(
        box.symbol_substitution(sub_map3));
    auto pauli_gadgets3 = sub_box3->get_pauli_gadgets();
    REQUIRE(pauli_gadgets3[0].coeff == ea);
    REQUIRE(pauli_gadgets3[1].coeff == eb);
    REQUIRE(pauli_gadgets3[2].coeff == sub_c);

    sub_map1.merge(sub_map2);
    sub_map1.merge(sub_map3);
    auto sub_box4 = std::dynamic_pointer_cast<const PauliExpCommutingSetBox>(
        box.symbol_substitution(sub_map1));
    auto pauli_gadgets4 = sub_box4->get_pauli_gadgets();
    REQUIRE(pauli_gadgets4[0].coeff == sub_a);
    REQUIRE(pauli_gadgets4[1].coeff == sub_b);
    REQUIRE(pauli_gadgets4[2].coeff == sub_c);
  }
}

}  // namespace test_PauliExpBoxes
}  // namespace tket
