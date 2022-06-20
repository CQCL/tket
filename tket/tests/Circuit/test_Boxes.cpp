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
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../testutil.hpp"
#include "Circuit/CircUtils.hpp"
#include "Circuit/Circuit.hpp"
#include "Converters/PhasePoly.hpp"
#include "Eigen/src/Core/Matrix.h"
#include "Gate/SymTable.hpp"
#include "Simulation/CircuitSimulator.hpp"

using Catch::Matchers::StartsWith;

namespace tket {
namespace test_Boxes {

SCENARIO("CircBox requires simple circuits", "[boxes]") {
  Circuit circ(2);
  circ.add_op<unsigned>(OpType::Y, {0});
  circ.add_op<unsigned>(OpType::CX, {0, 1});
  REQUIRE(circ.is_simple());
  Qubit qb0(0);
  Qubit qb1(1);
  Qubit a0("a", 0);
  Qubit a1("a", 1);
  unit_map_t qubit_map = {{qb0, a0}, {qb1, a1}};
  circ.rename_units(qubit_map);
  REQUIRE(!circ.is_simple());
  REQUIRE_THROWS_AS(CircBox(circ), SimpleOnly);
}

SCENARIO("Using Boxes", "[boxes]") {
  GIVEN("CircBox manipulation") {
    // Small box
    Circuit u(2);
    u.add_op<unsigned>(OpType::Ry, -0.75, {0});
    u.add_op<unsigned>(OpType::CX, {0, 1});
    const CircBox ubox(u);
    Circuit v(2);
    v.add_box(ubox, {0, 1});
    {
      const auto raw_u_unitary = tket_sim::get_unitary(u);
      const auto v_unitary = tket_sim::get_unitary(v);
      CHECK(raw_u_unitary.isApprox(v_unitary));
    }
    Circuit c0(3);
    c0.add_op<unsigned>(OpType::Rx, 0.5, {0});
    c0.add_op<unsigned>(OpType::Ry, 1.5, {1});
    c0.add_op<unsigned>(OpType::Rz, 0.75, {2});
    c0.add_box(ubox, {1, 0});
    c0.add_op<unsigned>(OpType::CX, {1, 2});
    REQUIRE(c0.n_gates() == 5);
    CircBox c0box(c0);
    // Put them in a bigger circuit
    Circuit d(4, 3);
    d.add_box(c0box, {1, 2, 0});
    d.add_op<unsigned>(OpType::CX, {0, 3});
    REQUIRE(d.n_gates() == 2);
    d.add_box(c0box, {3, 2, 1});
    REQUIRE(d.n_gates() == 3);
    d.add_box(c0box, {2, 3, 1});
    REQUIRE(d.n_gates() == 4);
    // Box up the bigger circuit
    CircBox dbox(d);
    Circuit e(4, 3);
    e.add_box(dbox, {/*qbs*/ 0, 1, 2, 3, /*cbs*/ 0, 1, 2});
    e.add_box(dbox, {/*qbs*/ 1, 2, 3, 0, /*cbs*/ 1, 2, 0});
    REQUIRE(e.n_gates() == 2);
    REQUIRE(!e.is_symbolic());
    // A circuit equivalent to c0 without boxes
    Circuit c0a(3);
    c0a.add_op<unsigned>(OpType::Rx, 0.5, {0});
    c0a.add_op<unsigned>(OpType::Ry, 1.5, {1});
    c0a.add_op<unsigned>(OpType::Rz, 0.75, {2});
    c0a.add_op<unsigned>(OpType::Ry, -0.75, {1});
    c0a.add_op<unsigned>(OpType::CX, {1, 0});
    c0a.add_op<unsigned>(OpType::CX, {1, 2});
    // Check c0 and c0a are equivalent
    Eigen::MatrixXcd uc0 = tket_sim::get_unitary(c0);
    Eigen::MatrixXcd uc0a = tket_sim::get_unitary(c0a);
    REQUIRE((uc0 - uc0a).cwiseAbs().sum() < ERR_EPS);
  }
  GIVEN("Unitary1qBox manipulation") {
    // random 1qb gate
    Circuit setup(1);
    setup.add_op<unsigned>(OpType::TK1, {0.2374, 1.0353, 0.5372}, {0});
    Eigen::Matrix2cd m = get_matrix_from_circ(setup);
    Unitary1qBox mbox(m);
    Circuit c(1);
    c.add_box(mbox, qubit_vector_t{Qubit("q", 0)});
    REQUIRE(c.n_gates() == 1);
    // extract its circuit
    std::shared_ptr<Circuit> excirc = mbox.to_circuit();
    // check we extract the same circuit from the deserialized box
    VertexSet vset = c.get_gates_of_type(OpType::Unitary1qBox);
    REQUIRE(vset.size() == 1);
    Vertex v = *vset.begin();
    Op_ptr op = c.get_Op_ptr_from_Vertex(v);
    std::shared_ptr<const Unitary1qBox> b =
        std::dynamic_pointer_cast<const Unitary1qBox>(op);
    std::shared_ptr<Circuit> excirc1 = b->to_circuit();
    REQUIRE(*excirc1 == *excirc);
    // compose with inverse of box
    c.append(c.dagger());
    Eigen::MatrixXcd c1m = tket_sim::get_unitary(c);
    // check it's the identity
    REQUIRE((c1m - Eigen::Matrix2cd::Identity()).cwiseAbs().sum() < ERR_EPS);
  }
  GIVEN("Unitary2qBox manipulation") {
    // permutation matrix
    Eigen::Matrix4cd m;
    m << 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0;
    Unitary2qBox mbox(m);
    Circuit c(2);
    c.add_box(mbox, {0, 1});
    REQUIRE(c.n_gates() == 1);
    // make a more complicated 2-qubit circuit
    Circuit d(2);
    d.add_op<unsigned>(OpType::Rx, 0.2, {0});
    d.add_op<unsigned>(OpType::Ry, 1.2, {1});
    d.add_op<unsigned>(OpType::CX, {0, 1});
    d.add_op<unsigned>(OpType::Rz, 0.4, {1});
    d.add_op<unsigned>(OpType::H, {0});
    d.add_op<unsigned>(OpType::CX, {1, 0});
    // get its unitary
    Eigen::Matrix4cd dm = get_matrix_from_2qb_circ(d);
    // make a box out of this
    Unitary2qBox dbox(dm);
    // make this into a new circuit
    Circuit d1(2);
    d1.add_box(dbox, {0, 1});
    // compose with inverse of d
    d1.append(d.dagger());
    Eigen::MatrixXcd d1m = tket_sim::get_unitary(d1);
    // check it's the identity
    REQUIRE((d1m - Eigen::Matrix4cd::Identity()).cwiseAbs().sum() < ERR_EPS);
  }
  GIVEN("little-endian representation") {
    Eigen::Matrix4cd m0;
    m0 << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0;
    Unitary2qBox m0box(m0);
    Circuit c0(2);
    c0.add_box(m0box, {0, 1});
    Circuit c1(2);
    c1.add_op<unsigned>(OpType::CX, {0, 1});
    Eigen::Matrix4cd m1 = get_matrix_from_2qb_circ(c1);
    REQUIRE((m0 - m1).cwiseAbs().sum() < ERR_EPS);
  }
  GIVEN("ExpBox manipulation") {
    // random hermitian matrix
    Eigen::Matrix4cd A;
    A << 0., 1., 2., 3., 1., 2., 3. * i_, 4., 2., -3. * i_, 3, 2. - 3. * i_, 3.,
        4., 2. + 3. * i_, 5.;
    ExpBox ebox(A, -0.5);
    Circuit c(2);
    c.add_box(ebox, {0, 1});
    Eigen::Matrix4cd U = (+0.5 * i_ * A).exp();  // should be the inverse
    Unitary2qBox ubox(U);
    c.add_box(ubox, {0, 1});  // should act as the identity
    Eigen::MatrixXcd uc = tket_sim::get_unitary(c);
    REQUIRE((uc - Eigen::Matrix4cd::Identity()).cwiseAbs().sum() < ERR_EPS);
  }
}

SCENARIO("Pauli gadgets", "[boxes]") {
  GIVEN("X") {
    // ---PauliExpBox([X], t)----Rx(-t)--- should be the identity
    double t = 1.687029013593215;
    Circuit c(1);
    PauliExpBox pbox({Pauli::X}, t);
    c.add_box(pbox, uvec{0});
    c.add_op<unsigned>(OpType::Rx, -t, {0});
    Eigen::Matrix2Xcd u = tket_sim::get_unitary(c);
    REQUIRE((u - Eigen::Matrix2cd::Identity()).cwiseAbs().sum() < ERR_EPS);
  }
  GIVEN("Y") {
    // ---PauliExpBox([Y], t)----Ry(-t)--- should be the identity
    double t = 1.6791969622440162;
    Circuit c(1);
    PauliExpBox pbox({Pauli::Y}, t);
    c.add_box(pbox, uvec{0});
    c.add_op<unsigned>(OpType::Ry, -t, {0});
    Eigen::Matrix2Xcd u = tket_sim::get_unitary(c);
    REQUIRE((u - Eigen::Matrix2cd::Identity()).cwiseAbs().sum() < ERR_EPS);
  }
  GIVEN("Z") {
    // ---PauliExpBox([Z], t)----Rz(-t)--- should be the identity
    double t = 1.7811410013115163;
    Circuit c(1);
    PauliExpBox pbox({Pauli::Z}, t);
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
    PauliExpBox pbox({Pauli::I, Pauli::I}, t);
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
    PauliExpBox pbox({Pauli::I, Pauli::X}, t);
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
    PauliExpBox pbox({Pauli::I, Pauli::Y}, t);
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
    PauliExpBox pbox({Pauli::I, Pauli::Z}, t);
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
    PauliExpBox pbox({Pauli::X, Pauli::I}, t);
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
    PauliExpBox pbox({Pauli::X, Pauli::X}, t);
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
    PauliExpBox pbox({Pauli::X, Pauli::Y}, t);
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
    PauliExpBox pbox({Pauli::X, Pauli::Z}, t);
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
    PauliExpBox pbox({Pauli::Y, Pauli::I}, t);
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
    PauliExpBox pbox({Pauli::Y, Pauli::X}, t);
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
    PauliExpBox pbox({Pauli::Y, Pauli::Y}, t);
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
    PauliExpBox pbox({Pauli::Y, Pauli::Z}, t);
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
    PauliExpBox pbox({Pauli::Z, Pauli::I}, t);
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
    PauliExpBox pbox({Pauli::Z, Pauli::X}, t);
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
    PauliExpBox pbox({Pauli::Z, Pauli::Y}, t);
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
    PauliExpBox pbox({Pauli::Z, Pauli::Z}, t);
    c.add_box(pbox, {0, 1});
    Eigen::MatrixXcd u = tket_sim::get_unitary(c);
    REQUIRE((u - Eigen::Matrix4cd::Identity()).cwiseAbs().sum() < ERR_EPS);
  }
  GIVEN("complex coefficient") {
    Expr ei{SymEngine::I};
    PauliExpBox pebox({Pauli::Z}, ei);
    Expr p = pebox.get_phase();
    REQUIRE(p == ei);
  }
}

SCENARIO("box daggers", "[boxes]") {
  GIVEN("a circuit made of various boxes") {
    // CircuitBox
    Circuit c0(2);
    c0.add_op<unsigned>(OpType::Ry, -0.75, {0});
    c0.add_op<unsigned>(OpType::CX, {0, 1});
    CircBox cbox(c0);
    // Unitary2qBox
    Eigen::Matrix4cd m;
    m << 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0;
    Unitary2qBox ubox(m);
    // ExpBox
    Eigen::Matrix4cd A;
    A << 0., 1., 2., 3., 1., 2., 3. * i_, 4., 2., -3. * i_, 3, 2. - 3. * i_, 3.,
        4., 2. + 3. * i_, 5.;
    ExpBox ebox(A, -0.5);
    // PauliExpBox
    PauliExpBox pbox({Pauli::X, Pauli::Y, Pauli::Z}, 0.8);

    // Put all these boxes into a circuit
    Circuit w(3);
    w.add_op<unsigned>(OpType::Rx, 0.5, {0});
    w.add_op<unsigned>(OpType::CX, {0, 1});
    w.add_box(cbox, {1, 2});
    w.add_box(ubox, {1, 0});
    w.add_box(ebox, {2, 1});
    w.add_box(pbox, {1, 2, 0});

    // Compute the dagger
    Circuit wdag = w.dagger();

    // Check dagger is correct
    w.append(wdag);
    Eigen::MatrixXcd u = tket_sim::get_unitary(w);
    REQUIRE((u - Eigen::MatrixXcd::Identity(8, 8)).cwiseAbs().sum() < ERR_EPS);
  }
}

SCENARIO("QControlBox", "[boxes]") {
  GIVEN("controlled X") {
    Op_ptr op = get_op_ptr(OpType::X);
    QControlBox qcbox(op);
    std::shared_ptr<Circuit> c = qcbox.to_circuit();
    Circuit expected(2);
    expected.add_op<unsigned>(OpType::CX, {0, 1});
    REQUIRE(*c == expected);
  }
  GIVEN("controlled CX") {
    Op_ptr op = get_op_ptr(OpType::CX);
    QControlBox qcbox(op);
    std::shared_ptr<Circuit> c = qcbox.to_circuit();
    Circuit expected(3);
    expected.add_op<unsigned>(OpType::CCX, {0, 1, 2});
    REQUIRE(*c == expected);
  }
  GIVEN("controlled CCX") {
    Op_ptr op = get_op_ptr(OpType::CCX);
    QControlBox qcbox(op);
    std::shared_ptr<Circuit> c = qcbox.to_circuit();
    Circuit expected(4);
    expected.add_op<unsigned>(OpType::CnX, {0, 1, 2, 3});
    REQUIRE(*c == expected);
  }
  GIVEN("controlled CnX") {
    Circuit c0(4);
    c0.add_op<unsigned>(OpType::CnX, {0, 1, 2, 3});
    Op_ptr op = c0.get_commands()[0].get_op_ptr();
    QControlBox qcbox(op);
    std::shared_ptr<Circuit> c = qcbox.to_circuit();
    Circuit expected(5);
    expected.add_op<unsigned>(OpType::CnX, {0, 1, 2, 3, 4});
    REQUIRE(*c == expected);
  }
  GIVEN("controlled Rz") {
    double a = 0.125;
    Circuit c0(1);
    c0.add_op<unsigned>(OpType::Rz, a, {0});
    Op_ptr op = c0.get_commands()[0].get_op_ptr();
    QControlBox qcbox(op);
    std::shared_ptr<Circuit> c = qcbox.to_circuit();
    const Eigen::MatrixXcd U = tket_sim::get_unitary(*c);
    Eigen::MatrixXcd V = Eigen::MatrixXcd::Identity(4, 4);
    V(2, 2) = exp(-0.5 * i_ * PI * a);
    V(3, 3) = exp(0.5 * i_ * PI * a);
    REQUIRE(U.isApprox(V));
  }
  GIVEN("controlled Rx") {
    double a = 0.125;
    Circuit c0(1);
    c0.add_op<unsigned>(OpType::Rx, a, {0});
    Op_ptr op = c0.get_commands()[0].get_op_ptr();
    QControlBox qcbox(op);
    std::shared_ptr<Circuit> c = qcbox.to_circuit();
    const Eigen::MatrixXcd U = tket_sim::get_unitary(*c);
    Eigen::MatrixXcd V = Eigen::MatrixXcd::Identity(4, 4);
    V(2, 2) = std::cos(0.5 * PI * a);
    V(2, 3) = i_ * std::sin(-0.5 * PI * a);
    V(3, 2) = i_ * std::sin(-0.5 * PI * a);
    V(3, 3) = std::cos(0.5 * PI * a);
    REQUIRE(U.isApprox(V));
  }
  GIVEN("controlled Ry") {
    double a = 0.125;
    Circuit c0(1);
    c0.add_op<unsigned>(OpType::Ry, a, {0});
    Op_ptr op = c0.get_commands()[0].get_op_ptr();
    QControlBox qcbox(op);
    std::shared_ptr<Circuit> c = qcbox.to_circuit();
    const Eigen::MatrixXcd U = tket_sim::get_unitary(*c);
    Eigen::MatrixXcd V = Eigen::MatrixXcd::Identity(4, 4);
    V(2, 2) = std::cos(0.5 * PI * a);
    V(2, 3) = std::sin(-0.5 * PI * a);
    V(3, 2) = std::sin(0.5 * PI * a);
    V(3, 3) = std::cos(0.5 * PI * a);
    REQUIRE(U.isApprox(V));
  }
  GIVEN("controlled S") {
    Op_ptr op = get_op_ptr(OpType::S);
    QControlBox qcbox(op);
    std::shared_ptr<Circuit> c = qcbox.to_circuit();
    const Eigen::MatrixXcd U = tket_sim::get_unitary(*c);
    Eigen::MatrixXcd V = Eigen::MatrixXcd::Identity(4, 4);
    V(2, 2) = 1;
    V(3, 3) = i_;
    REQUIRE(U.isApprox(V));
  }
  GIVEN("controlled V") {
    const double sq = 1 / std::sqrt(2.);
    Op_ptr op = get_op_ptr(OpType::V);
    QControlBox qcbox(op);
    std::shared_ptr<Circuit> c = qcbox.to_circuit();
    const Eigen::MatrixXcd U = tket_sim::get_unitary(*c);
    Eigen::MatrixXcd V = Eigen::MatrixXcd::Identity(4, 4);
    V(2, 2) = sq;
    V(2, 3) = sq * -i_;
    V(3, 2) = sq * -i_;
    V(3, 3) = sq;
    REQUIRE(U.isApprox(V));
  }
  GIVEN("controlled SX") {
    Op_ptr op = get_op_ptr(OpType::SX);
    QControlBox qcbox(op);
    std::shared_ptr<Circuit> c = qcbox.to_circuit();
    const Eigen::MatrixXcd U = tket_sim::get_unitary(*c);
    Eigen::MatrixXcd V = Eigen::MatrixXcd::Identity(4, 4);
    V(2, 2) = 0.5 * (1. + i_);
    V(2, 3) = 0.5 * (1. - i_);
    V(3, 2) = 0.5 * (1. - i_);
    V(3, 3) = 0.5 * (1. + i_);
    REQUIRE(U.isApprox(V));
  }
  GIVEN("controlled Sycamore") {
    Op_ptr op = get_op_ptr(OpType::Sycamore);
    QControlBox qcbox(op);
    std::shared_ptr<Circuit> c = qcbox.to_circuit();
    const Eigen::MatrixXcd U = tket_sim::get_unitary(*c);
    Eigen::MatrixXcd V = Eigen::MatrixXcd::Identity(8, 8);
    V(5, 5) = V(6, 6) = 0;
    V(5, 6) = V(6, 5) = -i_;
    V(7, 7) = exp(-i_ * PI / 6.);
    REQUIRE(U.isApprox(V));
  }
  GIVEN("2-controlled X") {
    Op_ptr op = get_op_ptr(OpType::X);
    QControlBox qcbox(op, 2);
    std::shared_ptr<Circuit> c = qcbox.to_circuit();
    Circuit expected(3);
    expected.add_op<unsigned>(OpType::CCX, {0, 1, 2});
    REQUIRE(*c == expected);
  }
  GIVEN("controlled CircBox") {
    Circuit c0(2);
    c0.add_op<unsigned>(OpType::H, {0});
    c0.add_op<unsigned>(OpType::CX, {0, 1});
    const Eigen::MatrixXcd U0 = tket_sim::get_unitary(c0);
    CircBox cbox(c0);
    Op_ptr op = std::make_shared<CircBox>(cbox);
    QControlBox qcbox(op);
    std::shared_ptr<Circuit> c = qcbox.to_circuit();
    const Eigen::MatrixXcd U = tket_sim::get_unitary(*c);
    Eigen::MatrixXcd V = Eigen::MatrixXcd::Identity(8, 8);
    for (unsigned i = 0; i < 4; i++) {
      for (unsigned j = 0; j < 4; j++) {
        V(4 + i, 4 + j) = U0(i, j);
      }
    }
    REQUIRE(U.isApprox(V));
  }
  GIVEN("controlled Unitary1qBox") {
    Circuit c0(1);
    c0.add_op<unsigned>(OpType::TK1, {0.6, 0.7, 0.8}, {0});
    c0.add_phase(0.9);
    Eigen::Matrix2cd m0 = get_matrix_from_circ(c0);
    Unitary1qBox mbox(m0);
    Op_ptr op = std::make_shared<Unitary1qBox>(mbox);
    QControlBox qcbox(op);
    std::shared_ptr<Circuit> c = qcbox.to_circuit();
    const Eigen::MatrixXcd U = tket_sim::get_unitary(*c);
    Eigen::MatrixXcd V = Eigen::MatrixXcd::Identity(4, 4);
    for (unsigned i = 0; i < 2; i++) {
      for (unsigned j = 0; j < 2; j++) {
        V(2 + i, 2 + j) = m0(i, j);
      }
    }
    REQUIRE(U.isApprox(V));
  }
  GIVEN("controlled Unitary2qBox") {
    Circuit c0(2);
    c0.add_op<unsigned>(OpType::Rx, 0.2, {0});
    c0.add_op<unsigned>(OpType::Ry, 1.2, {1});
    c0.add_op<unsigned>(OpType::CX, {0, 1});
    c0.add_op<unsigned>(OpType::Rz, 0.4, {1});
    c0.add_op<unsigned>(OpType::H, {0});
    c0.add_op<unsigned>(OpType::CX, {1, 0});
    Eigen::Matrix4cd m0 = get_matrix_from_2qb_circ(c0);
    Unitary2qBox ubox(m0);
    Op_ptr op = std::make_shared<Unitary2qBox>(ubox);
    QControlBox qcbox(op);
    std::shared_ptr<Circuit> c = qcbox.to_circuit();
    const Eigen::MatrixXcd U = tket_sim::get_unitary(*c);
    Eigen::MatrixXcd V = Eigen::MatrixXcd::Identity(8, 8);
    for (unsigned i = 0; i < 4; i++) {
      for (unsigned j = 0; j < 4; j++) {
        V(4 + i, 4 + j) = m0(i, j);
      }
    }
    REQUIRE(U.isApprox(V));
  }
  GIVEN("2-controlled Unitary2qBox") {
    // https://cqc.atlassian.net/browse/TKET-1651
    Eigen::Matrix4cd M{{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, -1}};
    Unitary2qBox ubox(M);
    Op_ptr op = std::make_shared<Unitary2qBox>(ubox);
    QControlBox qcbox(op, 2);
    std::shared_ptr<Circuit> c = qcbox.to_circuit();
    Eigen::MatrixXcd U = tket_sim::get_unitary(*c);
    Eigen::MatrixXcd V = Eigen::MatrixXcd::Identity(16, 16);
    V(15, 15) = -1;
    REQUIRE(U.isApprox(V));
  }
  GIVEN("controlled symbolic operation") {
    Sym s = SymEngine::symbol("a");
    Expr a = Expr(s);
    Op_ptr op = get_op_ptr(OpType::Rx, a);
    QControlBox qcbox(op);
    Circuit c(*qcbox.to_circuit());
    double v = 0.125;
    double x = std::cos(0.5 * PI * v), y = std::sin(0.5 * PI * v);
    symbol_map_t map = {{s, v}};
    c.symbol_substitution(map);
    const Eigen::MatrixXcd U = tket_sim::get_unitary(c);
    Eigen::Matrix4cd V = Eigen::Matrix4cd::Identity();
    V(2, 2) = V(3, 3) = x;
    V(2, 3) = V(3, 2) = -i_ * y;
    REQUIRE(U.isApprox(V));
  }
  GIVEN("nested QControlBox") {
    Op_ptr op = get_op_ptr(OpType::S);
    QControlBox qcbox(op);
    Circuit c(2);
    c.add_op<unsigned>(OpType::H, {0});
    c.add_box(qcbox, {0, 1});
    const Eigen::MatrixXcd U = tket_sim::get_unitary(c);
    CircBox cbox(c);
    Op_ptr op1 = std::make_shared<CircBox>(cbox);
    QControlBox qcbox1(op1);
    std::shared_ptr<Circuit> c1 = qcbox1.to_circuit();
    const Eigen::MatrixXcd U1 = tket_sim::get_unitary(*c1);
    Eigen::MatrixXcd V = Eigen::MatrixXcd::Identity(8, 8);
    for (unsigned i = 0; i < 4; i++) {
      for (unsigned j = 0; j < 4; j++) {
        V(4 + i, 4 + j) = U(i, j);
      }
    }
    REQUIRE(U1.isApprox(V));
  }
}

SCENARIO("Unitary3qBox", "[boxes]") {
  GIVEN("An 8x8 permutation matrix") {
    Eigen::MatrixXcd U = Eigen::MatrixXcd::Zero(8, 8);
    U(0, 3) = 1;
    U(1, 1) = 1;
    U(2, 7) = 1;
    U(3, 5) = 1;
    U(4, 0) = 1;
    U(5, 4) = 1;
    U(6, 2) = 1;
    U(7, 6) = 1;
    Unitary3qBox ubox(U);
    std::shared_ptr<const Circuit> c = ubox.to_circuit();
    REQUIRE(c->count_gates(OpType::CX) <= 24);
    Eigen::MatrixXcd U1 = tket_sim::get_unitary(*c);
    REQUIRE(U1.isApprox(U));
  }
}

SCENARIO("Checking equality", "[boxes]") {
  GIVEN("Some different types") {
    Circuit u(2);
    u.add_op<unsigned>(OpType::Rz, -0.75, {0});
    u.add_op<unsigned>(OpType::CX, {0, 1});
    const CircBox ubox(u);
    Eigen::Matrix4cd m = tket_sim::get_unitary(u);
    Unitary2qBox mbox(m);
    REQUIRE(ubox != mbox);

    Op_ptr op1 = get_op_ptr(OpType::X);
    Op_ptr op2 = get_op_ptr(OpType::Z);
    REQUIRE(op1 == op1);
    REQUIRE(op1 != op2);

    Eigen::Matrix4cd A;
    A << 0., 1., 2., 3., 1., 2., 3. * i_, 4., 2., -3. * i_, 3, 2. - 3. * i_, 3.,
        4., 2. + 3. * i_, 5.;
    ExpBox ebox(A, -0.5);
    REQUIRE(ebox != mbox);

    PhasePolyBox ppbox(u);
    REQUIRE(ppbox != mbox);
  }
  GIVEN("CircBoxes") {
    Circuit u(2);
    u.add_op<unsigned>(OpType::Ry, -0.75, {0});
    u.add_op<unsigned>(OpType::CX, {0, 1});
    const CircBox ubox(u);

    WHEN("both arguments are equal") { REQUIRE(ubox == ubox); }
    WHEN("both arguments are different") { REQUIRE(ubox != CircBox(u)); }
  }
  GIVEN("Unitary1qBox") {
    Circuit setup(1);
    setup.add_op<unsigned>(OpType::TK1, {0.2374, 1.0353, 0.5372}, {0});
    Eigen::Matrix2cd m = tket_sim::get_unitary(setup);
    Unitary1qBox mbox(m);

    WHEN("both arguments are equal") { REQUIRE(mbox == mbox); }
    WHEN("both arguments are different") {
      setup.add_op<unsigned>(OpType::TK1, {0.2374, 1.0353, 0.5372}, {0});
      Eigen::Matrix2cd m2 = tket_sim::get_unitary(setup);
      Unitary1qBox mbox2(m2);
      REQUIRE(mbox != mbox2);
    }
  }
  GIVEN("Unitary2qBox") {
    Circuit setup(2);
    setup.add_op<unsigned>(OpType::TK1, {0.2374, 1.0353, 0.5372}, {0});
    setup.add_op<unsigned>(OpType::CX, {0, 1});
    Eigen::Matrix4cd m = tket_sim::get_unitary(setup);
    Unitary2qBox mbox(m);

    WHEN("both arguments are equal") { REQUIRE(mbox == mbox); }
    WHEN("both arguments are different") {
      setup.add_op<unsigned>(OpType::CX, {1, 0});
      Eigen::Matrix4cd m2 = tket_sim::get_unitary(setup);
      Unitary2qBox mbox2(m2);
      REQUIRE(mbox != mbox2);
    }
  }
  GIVEN("Unitary3qBox") {
    Circuit setup(3);
    setup.add_op<unsigned>(OpType::TK1, {0.2374, 1.0353, 0.5372}, {0});
    setup.add_op<unsigned>(OpType::CX, {0, 1});
    setup.add_op<unsigned>(OpType::CX, {1, 2});
    Eigen::MatrixXcd m = tket_sim::get_unitary(setup);
    Unitary3qBox mbox(m);

    WHEN("both arguments are equal") { REQUIRE(mbox == mbox); }
    WHEN("both arguments are different") {
      setup.add_op<unsigned>(OpType::CX, {0, 2});
      Eigen::MatrixXcd m2 = tket_sim::get_unitary(setup);
      Unitary3qBox mbox2(m2);
      REQUIRE(mbox != mbox2);
    }
  }
  GIVEN("ExpBox") {
    // random hermitian matrix
    Eigen::Matrix4cd A;
    A << 0., 1., 2., 3., 1., 2., 3. * i_, 4., 2., -3. * i_, 3, 2. - 3. * i_, 3.,
        4., 2. + 3. * i_, 5.;
    ExpBox ebox(A, -0.5);
    WHEN("both arguments are equal") { REQUIRE(ebox == ebox); }
    WHEN("both arguments are different") {
      ExpBox ebox2(A, -0.2);
      REQUIRE(ebox != ebox2);
    }
  }
  GIVEN("Pauli gadgets") {
    double t = 1.687029013593215;
    PauliExpBox pbox({Pauli::X}, t);
    WHEN("both arguments are equal") { REQUIRE(pbox == pbox); }
    WHEN("both arguments are different") {
      PauliExpBox pbox2({Pauli::Y}, t);
      REQUIRE(pbox != pbox2);
    }
  }
  GIVEN("QControlBox") {
    Op_ptr op = get_op_ptr(OpType::X);
    QControlBox qcbox(op);
    WHEN("both arguments are equal") { REQUIRE(qcbox == qcbox); }
    WHEN("both arguments are different") {
      Op_ptr op2 = get_op_ptr(OpType::Y);
      QControlBox qcbox2(op2);
      REQUIRE(qcbox != qcbox2);
    }
  }
  GIVEN("PhasePolyBox") {
    Circuit u(2);
    u.add_op<unsigned>(OpType::Rz, -0.75, {0});
    u.add_op<unsigned>(OpType::CX, {0, 1});
    PhasePolyBox ppbox(u);
    WHEN("both arguments are equal") { REQUIRE(ppbox == ppbox); }
    WHEN("both arguments are different") {
      u.add_op<unsigned>(OpType::CX, {1, 0});
      PhasePolyBox ppbox2(u);
      REQUIRE(ppbox != ppbox2);
    }
  }
  GIVEN("CustomGate") {
    Circuit setup(1);
    Sym a = SymTable::fresh_symbol("a");
    Expr ea(a);

    // "random" 1qb gate.
    const double param1 = 1.23323;
    const double param2 = 0.42323;
    const double param3 = 0.34212;
    const std::string name1{"gate name1"};
    const std::string name2{"gate name2"};
    setup.add_op<unsigned>(OpType::TK1, {ea, param1, param2}, {0});

    composite_def_ptr_t def1 = CompositeGateDef::define_gate(name1, setup, {a});
    composite_def_ptr_t def2 = CompositeGateDef::define_gate(name2, setup, {a});
    const CustomGate g1(def1, {param3});
    const CustomGate g1_repeated(def1, {param3});
    const CustomGate g1_wrong(def1, {param1});
    const CustomGate g2(def2, {param3});

    // Check that all IDs are different.
    const std::set<boost::uuids::uuid> ids{
        g1.get_id(), g1_repeated.get_id(), g1_wrong.get_id(), g2.get_id()};
    CHECK(ids.size() == 4);
    CHECK(g1 == g1);
    CHECK(g1 == g1_repeated);
    CHECK(g1 != g2);
    CHECK(g1 != g1_wrong);
    CHECK(g1_repeated != g1_wrong);
    CHECK_THROWS_AS(CustomGate(nullptr, {param3}), std::runtime_error);
  }
}

SCENARIO("Checking box names", "[boxes]") {
  GIVEN("CustomGate without parameters") {
    Circuit setup(1);
    setup.add_op<unsigned>(OpType::TK1, {0.3333, 1.111, 0.5555}, {0});
    const std::string name("gate without params");
    composite_def_ptr_t def = CompositeGateDef::define_gate(name, setup, {});
    CustomGate g(def, {});
    CHECK(g.get_name() == name);
  }
  GIVEN("CustomGate with 1 parameter") {
    Circuit setup(1);
    Sym a = SymTable::fresh_symbol("a");
    Expr ea(a);
    setup.add_op<unsigned>(OpType::TK1, {ea, 0.3333, 1.111}, {0});
    const std::string prefix("gate with params");
    composite_def_ptr_t def = CompositeGateDef::define_gate(prefix, setup, {a});
    CustomGate g(def, {0.4444});

    // Of course, 0.4444 is NOT exactly represented by a double,
    // so it might print something like 0.4443999... or 0.4440000...1.
    // This test will still pass even if so.
    CHECK_THAT(g.get_name(), StartsWith(prefix + "(0.444"));
  }
  GIVEN("CustomGate with 3 parameters") {
    Circuit setup(1);
    Sym a = SymTable::fresh_symbol("a");
    Sym b = SymTable::fresh_symbol("b");
    Sym c = SymTable::fresh_symbol("c");
    Expr ea(a);
    Expr eb(b);
    Expr ec(c);
    setup.add_op<unsigned>(OpType::TK1, {ea, eb, ec}, {0});
    const std::string prefix("gate with 3 params");
    composite_def_ptr_t def =
        CompositeGateDef::define_gate(prefix, setup, {a, b, c});
    CustomGate g(def, {0.1111, 0.2222, 0.4444});
    const std::string name = g.get_name();
    CHECK(name == "gate with 3 params(0.1111,0.2222,0.4444)");
  }
}

}  // namespace test_Boxes
}  // namespace tket
