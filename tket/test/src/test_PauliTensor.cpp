// Copyright Quantinuum
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

#include "tket/Utils/PauliTensor.hpp"

namespace tket {
namespace test_PauliTensor {

SCENARIO("Testing equality of sparse PauliTensor variants") {
  Qubit q0 = Qubit("q", 0);
  Qubit q1 = Qubit("q", 1);
  Qubit q2 = Qubit("r", 0);
  Qubit q3 = Qubit("s");
  Qubit q4 = Qubit("t", 0, 1);
  Qubit q5 = Qubit("t", 0, 0);
  GIVEN("Two exactly identical Pauli strings") {
    QubitPauliMap map = {
        {q0, Pauli::I}, {q1, Pauli::X}, {q2, Pauli::Y}, {q3, Pauli::Z}};
    SpCxPauliTensor a(map, i_);
    SpCxPauliTensor b(map, i_);
    REQUIRE(a == b);
    THEN("We add some extra Is on each one") {
      a.set(q4, Pauli::I);
      b.set(q5, Pauli::I);
      REQUIRE(a == b);
    }
  }
  GIVEN("Two Pauli strings with different Paulis but same coefficient") {
    SpPauliString a({{q0, Pauli::X}});
    SpPauliString b({{q0, Pauli::Y}});
    REQUIRE(a != b);
    REQUIRE(a < b);
  }
  GIVEN("Two Pauli strings with disjoint Paulis but same coefficient") {
    SpPauliString a(q0, Pauli::X);
    SpPauliString b(q1, Pauli::X);
    REQUIRE(a != b);
    REQUIRE(b < a);
  }
  GIVEN("Two Pauli strings with same Paulis but different coefficient") {
    SpCxPauliTensor a(q0, Pauli::X, 1.);
    SpCxPauliTensor b(q0, Pauli::X, i_);
    REQUIRE(a != b);
    REQUIRE(b < a);
  }
  GIVEN("Two completely different Pauli strings") {
    QubitPauliMap qpm_a(
        {{q0, Pauli::I}, {q1, Pauli::X}, {q2, Pauli::Y}, {q3, Pauli::Z}});
    QubitPauliMap qpm_b(
        {{q0, Pauli::X}, {q1, Pauli::I}, {q2, Pauli::Z}, {q4, Pauli::Y}});
    SpCxPauliTensor a(qpm_a, 1.);
    SpCxPauliTensor b(qpm_b, i_);
    REQUIRE(a != b);
    REQUIRE(a < b);
  }
}

SCENARIO("Testing equality of dense PauliTensor variants") {
  GIVEN("Two exactly identical Pauli strings") {
    DensePauliMap map = {Pauli::I, Pauli::X, Pauli::Y, Pauli::Z};
    CxPauliTensor a(map, i_);
    CxPauliTensor b(map, i_);
    REQUIRE(a == b);
    THEN("We add some extra Is on each one") {
      a.set(4, Pauli::I);
      b.set(5, Pauli::I);
      REQUIRE(a == b);
    }
  }
  GIVEN("Two Pauli strings with different Paulis but same coefficient") {
    PauliString a({Pauli::X});
    PauliString b({Pauli::Y});
    REQUIRE(a != b);
    REQUIRE(a < b);
  }
  GIVEN("Two Pauli strings with disjoint Paulis but same coefficient") {
    PauliString a({Pauli::X});
    PauliString b({Pauli::I, Pauli::X});
    REQUIRE(a != b);
    REQUIRE(b < a);
  }
  GIVEN("Two Pauli strings with same Paulis but different coefficient") {
    CxPauliTensor a({Pauli::X}, 1.);
    CxPauliTensor b({Pauli::X}, i_);
    REQUIRE(a != b);
    REQUIRE(b < a);
  }
  GIVEN("Two completely different Pauli strings") {
    DensePauliMap qpm_a({Pauli::I, Pauli::X, Pauli::Y, Pauli::Z});
    DensePauliMap qpm_b({Pauli::X, Pauli::I, Pauli::Z, Pauli::Y});
    CxPauliTensor a(qpm_a, 1.);
    CxPauliTensor b(qpm_b, i_);
    REQUIRE(a != b);
    REQUIRE(a < b);
  }
}

SCENARIO("Testing casting between different PauliTensor variants") {
  GIVEN("Casting between different coefficient types") {
    SpPauliString a1{};
    SpPauliStabiliser b1{};
    SpPauliStabiliser bi({}, 1);
    SpPauliStabiliser bm1({}, 2);
    SpPauliStabiliser bmi({}, 3);
    SpCxPauliTensor c1{};
    SpCxPauliTensor ci({}, i_);
    SpCxPauliTensor cm1({}, -1.);
    SpCxPauliTensor cmi({}, -i_);
    SpCxPauliTensor cval({}, 0.48 - 2.3 * i_);
    SpSymPauliTensor d1{};
    SpSymPauliTensor di({}, Expr(SymEngine::I));
    SpSymPauliTensor dm1({}, -1.);
    SpSymPauliTensor dmi({}, -Expr(SymEngine::I));
    SpSymPauliTensor dval({}, 0.48 - 2.3 * i_);
    SpSymPauliTensor dsym({}, Expr("a"));

    CHECK((SpPauliString)a1 == a1);
    CHECK((SpPauliString)b1 == a1);
    CHECK((SpPauliString)c1 == a1);
    CHECK((SpPauliString)d1 == a1);
    CHECK((SpPauliStabiliser)a1 == b1);
    CHECK((SpPauliStabiliser)b1 == b1);
    CHECK((SpPauliStabiliser)c1 == b1);
    CHECK((SpPauliStabiliser)d1 == b1);
    CHECK((SpCxPauliTensor)a1 == c1);
    CHECK((SpCxPauliTensor)b1 == c1);
    CHECK((SpCxPauliTensor)c1 == c1);
    CHECK((SpCxPauliTensor)d1 == c1);
    CHECK((SpSymPauliTensor)a1 == d1);
    CHECK((SpSymPauliTensor)b1 == d1);
    CHECK((SpSymPauliTensor)c1 == d1);
    CHECK((SpSymPauliTensor)d1 == d1);

    CHECK((SpPauliString)bi == a1);
    CHECK((SpPauliString)ci == a1);
    CHECK((SpPauliString)di == a1);
    CHECK((SpPauliStabiliser)bi == bi);
    CHECK((SpPauliStabiliser)ci == bi);
    CHECK((SpPauliStabiliser)di == bi);
    CHECK((SpCxPauliTensor)bi == ci);
    CHECK((SpCxPauliTensor)ci == ci);
    CHECK((SpCxPauliTensor)di == ci);
    CHECK((SpSymPauliTensor)bi == di);
    CHECK((SpSymPauliTensor)ci == di);
    CHECK((SpSymPauliTensor)di == di);

    CHECK((SpPauliString)bm1 == a1);
    CHECK((SpPauliString)cm1 == a1);
    CHECK((SpPauliString)dm1 == a1);
    CHECK((SpPauliStabiliser)bm1 == bm1);
    CHECK((SpPauliStabiliser)cm1 == bm1);
    CHECK((SpPauliStabiliser)dm1 == bm1);
    CHECK((SpCxPauliTensor)bm1 == cm1);
    CHECK((SpCxPauliTensor)cm1 == cm1);
    CHECK((SpCxPauliTensor)dm1 == cm1);
    CHECK((SpSymPauliTensor)bm1 == dm1);
    CHECK((SpSymPauliTensor)cm1 == dm1);
    CHECK((SpSymPauliTensor)dm1 == dm1);

    CHECK((SpPauliString)bmi == a1);
    CHECK((SpPauliString)cmi == a1);
    CHECK((SpPauliString)dmi == a1);
    CHECK((SpPauliStabiliser)bmi == bmi);
    CHECK((SpPauliStabiliser)cmi == bmi);
    CHECK((SpPauliStabiliser)dmi == bmi);
    CHECK((SpCxPauliTensor)bmi == cmi);
    CHECK((SpCxPauliTensor)cmi == cmi);
    CHECK((SpCxPauliTensor)dmi == cmi);
    CHECK((SpSymPauliTensor)bmi == dmi);
    CHECK((SpSymPauliTensor)cmi == dmi);
    CHECK((SpSymPauliTensor)dmi == dmi);

    CHECK((SpPauliString)cval == a1);
    CHECK((SpPauliString)dval == a1);
    REQUIRE_THROWS((SpPauliStabiliser)cval);
    REQUIRE_THROWS((SpPauliStabiliser)dval);
    CHECK((SpCxPauliTensor)cval == cval);
    CHECK((SpCxPauliTensor)dval == cval);
    CHECK((SpSymPauliTensor)cval == dval);
    CHECK((SpSymPauliTensor)dval == dval);

    CHECK((SpPauliString)dsym == a1);
    REQUIRE_THROWS((SpPauliStabiliser)dsym);
    REQUIRE_THROWS((SpCxPauliTensor)dsym);
    CHECK((SpSymPauliTensor)dsym == dsym);
  }
  GIVEN("Casting between different Pauli containers") {
    PauliString ps({Pauli::I, Pauli::X, Pauli::Y});
    SpPauliString sps({Qubit(1), Qubit(2)}, {Pauli::X, Pauli::Y});
    SpPauliString non_default(Qubit("a", 0), Pauli::Z);

    CHECK((SpPauliString)ps == sps);
    CHECK((SpPauliString)sps == sps);
    CHECK((PauliString)ps == ps);
    CHECK((PauliString)sps == ps);
    CHECK((SpPauliString)non_default == non_default);
    REQUIRE_THROWS((PauliString)non_default);
  }
  GIVEN("Casting coefficient, keeping dense container") {
    DensePauliMap dpm{Pauli::X, Pauli::I, Pauli::Z};
    CHECK((PauliStabiliser)PauliString(dpm) == PauliStabiliser(dpm));
    CHECK(
        (SymPauliTensor)CxPauliTensor(dpm, 0.87 + 1.2 * i_) ==
        SymPauliTensor(dpm, Expr(0.87 + 1.2 * i_)));
  }
}

SCENARIO("Qubit partitions") {
  GIVEN("Sparse PauliTensors") {
    std::vector<Qubit> qs{Qubit(0),           Qubit("a", 0), Qubit("a", 1),
                          Qubit("b", {0, 0}), Qubit("c", 4), Qubit("p", 12),
                          Qubit("anc", 0),    Qubit(2)};
    std::list<Qubit> qsl{qs.begin(), qs.end()};
    SpPauliString xxyyzzii(
        qsl, {Pauli::X, Pauli::X, Pauli::Y, Pauli::Y, Pauli::Z, Pauli::Z,
              Pauli::I, Pauli::I});
    SpPauliString ixyxyziz(
        qsl, {Pauli::I, Pauli::X, Pauli::Y, Pauli::X, Pauli::Y, Pauli::Z,
              Pauli::I, Pauli::Z});
    xxyyzzii.compress();

    // Common qubits should ignore Pauli::I matches
    CHECK(
        xxyyzzii.common_qubits(ixyxyziz) ==
        std::set<Qubit>{qs[1], qs[2], qs[5]});
    CHECK(
        xxyyzzii.conflicting_qubits(ixyxyziz) == std::set<Qubit>{qs[3], qs[4]});
    CHECK(xxyyzzii.own_qubits(ixyxyziz) == std::set<Qubit>{qs[0]});
    CHECK(ixyxyziz.own_qubits(xxyyzzii) == std::set<Qubit>{qs[7]});
  }
  GIVEN("Dense PauliTensors") {
    PauliString xxyyzzii(
        {Pauli::X, Pauli::X, Pauli::Y, Pauli::Y, Pauli::Z, Pauli::Z, Pauli::I});
    PauliString ixyxyziz(
        {Pauli::I, Pauli::X, Pauli::Y, Pauli::X, Pauli::Y, Pauli::Z, Pauli::I,
         Pauli::Z});

    // Common indices should ignore Pauli::I matches
    CHECK(xxyyzzii.common_indices(ixyxyziz) == std::set<unsigned>{1, 2, 5});
    CHECK(xxyyzzii.conflicting_indices(ixyxyziz) == std::set<unsigned>{3, 4});
    CHECK(xxyyzzii.own_indices(ixyxyziz) == std::set<unsigned>{0});
    CHECK(ixyxyziz.own_indices(xxyyzzii) == std::set<unsigned>{7});
  }
}

SCENARIO("String formatting of PauliTensor") {
  CHECK(SpPauliString().to_str() == "()");
  CHECK(SpPauliStabiliser({}, 0).to_str() == "()");
  CHECK(SpPauliStabiliser({}, 1).to_str() == "i*()");
  CHECK(SpPauliStabiliser({}, 2).to_str() == "-()");
  CHECK(SpPauliStabiliser({}, 3).to_str() == "-i*()");
  CHECK(SpCxPauliTensor({}, 1.).to_str() == "()");
  CHECK(SpCxPauliTensor({}, -1.).to_str() == "-()");
  CHECK(SpCxPauliTensor({}, 4.2 + 0.87 * i_).to_str() == "(4.2,0.87)*()");
  CHECK(SpSymPauliTensor({}, 1.).to_str() == "()");
  CHECK(SpSymPauliTensor({}, -1.).to_str() == "-()");
  CHECK(
      SpSymPauliTensor({}, 4.2 + 0.87 * Expr(SymEngine::I)).to_str() ==
      "(4.2 + 0.87*I)*()");
  CHECK(SpSymPauliTensor({}, Expr("2*a")).to_str() == "(2*a)*()");

  CHECK(
      SpPauliString({{Qubit("a", 2), Pauli::X},
                     {Qubit("a", 0), Pauli::Z},
                     {Qubit("b", 0), Pauli::I},
                     {Qubit("b", 1), Pauli::Y}})
          .to_str() == "(Za[0], Xa[2], Ib[0], Yb[1])");
  CHECK(
      PauliString({Pauli::I, Pauli::Z, Pauli::X, Pauli::Y, Pauli::I})
          .to_str() == "IZXYI");
  CHECK(PauliStabiliser({Pauli::X, Pauli::Y}, 2).to_str() == "-XY");
  CHECK(
      CxPauliTensor({Pauli::Z, Pauli::Z, Pauli::I}, 3.1 - 0.1 * i_).to_str() ==
      "(3.1,-0.1)*ZZI");
  CHECK(
      SymPauliTensor(DensePauliMap(5, Pauli::Y), Expr("k")).to_str() ==
      "(k)*YYYYY");
}

SCENARIO("Testing multiplication of sparse PauliTensor") {
  Qubit q0 = Qubit("q", 0);
  Qubit q1 = Qubit("q", 1);
  Qubit q2 = Qubit("r", 0);
  Qubit q3 = Qubit("s");
  Qubit q4 = Qubit("t", 0, 1);
  GIVEN("Two Pauli strings with disjoint non-trivial components") {
    SpPauliString a(q0, Pauli::X);
    SpPauliString b(q1, Pauli::Y);
    SpPauliString c({{q0, Pauli::X}, {q1, Pauli::Y}});
    REQUIRE((a * b) == c);
  }
  GIVEN("Multiplying by a trivial Pauli string") {
    SpCxPauliTensor a(q0, Pauli::X, 2.);
    SpCxPauliTensor b(q0, Pauli::X, 3. * i_);
    REQUIRE((a * SpCxPauliTensor({}, 1.5 * i_)) == b);
  }
  GIVEN("Two exactly identical Pauli strings") {
    QubitPauliMap map = {
        {q0, Pauli::I}, {q1, Pauli::X}, {q2, Pauli::Y}, {q3, Pauli::Z}};
    SpPauliStabiliser a(map, 3);
    SpPauliStabiliser b({}, 2);
    REQUIRE((a * a).get(q0) == Pauli::I);
    REQUIRE((a * a) == b);
  }
  GIVEN("Each individual Pauli combination") {
    SpPauliStabiliser I(q0, Pauli::I);
    SpPauliStabiliser X(q0, Pauli::X);
    SpPauliStabiliser Y(q0, Pauli::Y);
    SpPauliStabiliser Z(q0, Pauli::Z);
    SpPauliStabiliser i({}, 1);
    SpPauliStabiliser mi({}, 3);
    REQUIRE((I * I) == I);
    REQUIRE((I * X) == X);
    REQUIRE((I * Y) == Y);
    REQUIRE((I * Z) == Z);
    REQUIRE((X * I) == X);
    REQUIRE((X * X) == I);
    REQUIRE((X * Y) == (i * Z));
    REQUIRE((X * Z) == (mi * Y));
    REQUIRE((Y * I) == Y);
    REQUIRE((Y * X) == (mi * Z));
    REQUIRE((Y * Y) == I);
    REQUIRE((Y * Z) == (i * X));
    REQUIRE((Z * I) == Z);
    REQUIRE((Z * X) == (i * Y));
    REQUIRE((Z * Y) == (mi * X));
    REQUIRE((Z * Z) == I);
  }
  GIVEN("2*IXYZ(I) * -1.5i*XIZ(I)Y") {
    QubitPauliMap tensor_a(
        {{q0, Pauli::I}, {q1, Pauli::X}, {q2, Pauli::Y}, {q3, Pauli::Z}});
    QubitPauliMap tensor_b(
        {{q0, Pauli::X}, {q1, Pauli::I}, {q2, Pauli::Z}, {q4, Pauli::Y}});
    SpCxPauliTensor a(tensor_a, 2.);
    SpCxPauliTensor b(tensor_b, -1.5 * i_);
    QubitPauliMap tensor_c(
        {{q0, Pauli::X},
         {q1, Pauli::X},
         {q2, Pauli::X},
         {q3, Pauli::Z},
         {q4, Pauli::Y}});
    SpCxPauliTensor c(tensor_c, 3.);
    REQUIRE((a * b) == c);
  }
}

SCENARIO("Testing multiplication of dense PauliTensor") {
  GIVEN("Two Pauli strings with disjoint non-trivial components") {
    PauliString a({Pauli::X});
    PauliString b({Pauli::I, Pauli::Y});
    PauliString c({Pauli::X, Pauli::Y});
    REQUIRE((a * b) == c);
  }
  GIVEN("Multiplying by a trivial Pauli string") {
    CxPauliTensor a({Pauli::X}, 2.);
    CxPauliTensor b({Pauli::X}, 3. * i_);
    REQUIRE((a * CxPauliTensor({}, 1.5 * i_)) == b);
  }
  GIVEN("Two exactly identical Pauli strings") {
    DensePauliMap map = {Pauli::I, Pauli::X, Pauli::Y, Pauli::Z};
    PauliStabiliser a(map, 3);
    PauliStabiliser b({}, 2);
    REQUIRE((a * a) == b);
  }
  GIVEN("Each individual Pauli combination") {
    PauliStabiliser I({Pauli::I});
    PauliStabiliser X({Pauli::X});
    PauliStabiliser Y({Pauli::Y});
    PauliStabiliser Z({Pauli::Z});
    PauliStabiliser i({}, 1);
    PauliStabiliser mi({}, 3);
    REQUIRE((I * I) == I);
    REQUIRE((I * X) == X);
    REQUIRE((I * Y) == Y);
    REQUIRE((I * Z) == Z);
    REQUIRE((X * I) == X);
    REQUIRE((X * X) == I);
    REQUIRE((X * Y) == (i * Z));
    REQUIRE((X * Z) == (mi * Y));
    REQUIRE((Y * I) == Y);
    REQUIRE((Y * X) == (mi * Z));
    REQUIRE((Y * Y) == I);
    REQUIRE((Y * Z) == (i * X));
    REQUIRE((Z * I) == Z);
    REQUIRE((Z * X) == (i * Y));
    REQUIRE((Z * Y) == (mi * X));
    REQUIRE((Z * Z) == I);
  }
  GIVEN("2*IXYZ(I) * -1.5i*XIZ(I)Y") {
    DensePauliMap tensor_a({Pauli::I, Pauli::X, Pauli::Y, Pauli::Z});
    DensePauliMap tensor_b({Pauli::X, Pauli::I, Pauli::Z, Pauli::I, Pauli::Y});
    CxPauliTensor a(tensor_a, 2.);
    CxPauliTensor b(tensor_b, -1.5 * i_);
    DensePauliMap tensor_c({Pauli::X, Pauli::X, Pauli::X, Pauli::Z, Pauli::Y});
    CxPauliTensor c(tensor_c, 3.);
    REQUIRE((a * b) == c);
  }
}

SCENARIO("Test hashing for sparse PauliTensor") {
  GIVEN("Trivial strings") {
    SpPauliString qps1;
    SpPauliString qps2;
    REQUIRE(qps1.hash_value() == qps2.hash_value());
    WHEN("Add I Pauli") {
      qps1.set(Qubit(0), Pauli::I);
      REQUIRE(qps1.hash_value() == qps2.hash_value());
    }
  }
  GIVEN("Nontrivial strings") {
    QubitPauliMap qpm{
        {Qubit(0), Pauli::Z},
        {Qubit(1), Pauli::Y},
        {Qubit(2), Pauli::X},
        {Qubit(3), Pauli::I}};
    SpPauliString qps1(qpm);
    SpPauliString qps2(qpm);
    qps1.set(Qubit(4), Pauli::X);
    qps2.set(Qubit(4), Pauli::X);
    qps2.set(Qubit(5), Pauli::I);
    REQUIRE(qps1.hash_value() == qps2.hash_value());
  }
  GIVEN("Trivial tensor") {
    SpCxPauliTensor qpt1;
    SpCxPauliTensor qpt2;
    REQUIRE(qpt1.hash_value() == qpt2.hash_value());
    WHEN("Add I Pauli") {
      qpt1.set(Qubit(0), Pauli::I);
      REQUIRE(qpt1.hash_value() == qpt2.hash_value());
    }
  }
  GIVEN("Nontrivial tensors") {
    QubitPauliMap qpm{
        {Qubit(0), Pauli::Z},
        {Qubit(1), Pauli::Y},
        {Qubit(2), Pauli::X},
        {Qubit(3), Pauli::I}};
    SpSymPauliTensor qpt1(qpm, .5 * i_);
    SpSymPauliTensor qpt2(qpm, .5 * i_);
    qpt1.set(Qubit(4), Pauli::X);
    qpt2.set(Qubit(4), Pauli::X);
    qpt2.set(Qubit(5), Pauli::I);
    qpt2.set(Qubit(6), Pauli::I);
    REQUIRE(qpt1.hash_value() == qpt2.hash_value());
  }
}

SCENARIO("Test hashing for dense PauliTensor") {
  GIVEN("Trivial strings") {
    PauliString qps1;
    PauliString qps2;
    REQUIRE(qps1.hash_value() == qps2.hash_value());
    WHEN("Add I Pauli") {
      qps1.set(0, Pauli::I);
      REQUIRE(qps1.hash_value() == qps2.hash_value());
    }
  }
  GIVEN("Nontrivial strings") {
    DensePauliMap qpm{Pauli::Z, Pauli::Y, Pauli::X, Pauli::I};
    PauliString qps1(qpm);
    PauliString qps2(qpm);
    qps1.set(4, Pauli::X);
    qps2.set(4, Pauli::X);
    qps2.set(5, Pauli::I);
    REQUIRE(qps1.hash_value() == qps2.hash_value());
  }
  GIVEN("Stabilisers") {
    DensePauliMap pm{Pauli::Z, Pauli::Y, Pauli::X, Pauli::I};
    PauliStabiliser ps1(pm, 2);
    PauliStabiliser ps2(pm, 6);
    REQUIRE(ps1.hash_value() == ps2.hash_value());
  }
  GIVEN("Trivial tensor") {
    CxPauliTensor qpt1;
    CxPauliTensor qpt2;
    REQUIRE(qpt1.hash_value() == qpt2.hash_value());
    WHEN("Add I Pauli") {
      qpt1.set(0, Pauli::I);
      REQUIRE(qpt1.hash_value() == qpt2.hash_value());
    }
  }
  GIVEN("Nontrivial tensors") {
    DensePauliMap qpm{Pauli::Z, Pauli::Y, Pauli::X, Pauli::I};
    SymPauliTensor qpt1(qpm, .5 * i_);
    SymPauliTensor qpt2(qpm, .5 * i_);
    qpt1.set(4, Pauli::X);
    qpt2.set(4, Pauli::X);
    qpt2.set(5, Pauli::I);
    qpt2.set(6, Pauli::I);
    REQUIRE(qpt1.hash_value() == qpt2.hash_value());
  }
}

SCENARIO("json serialisation of PauliTensor") {
  PauliString xyz({Pauli::X, Pauli::Y, Pauli::Z});
  nlohmann::json j = xyz;
  CHECK(j.get<PauliString>() == xyz);
  SpPauliString za(Qubit("a", 0), Pauli::Z);
  j = za;
  CHECK(j.get<SpPauliString>() == za);
  PauliStabiliser zz({Pauli::Z, Pauli::Z}, 3);
  j = zz;
  CHECK(j.get<PauliStabiliser>() == zz);
  SpPauliStabiliser ziz({Pauli::Z, Pauli::I, Pauli::Z}, 2);
  j = ziz;
  CHECK(j.get<SpPauliStabiliser>() == ziz);
  CxPauliTensor yiy({Pauli::Y, Pauli::I, Pauli::Y}, 0.2 * i_);
  j = yiy;
  CHECK(j.get<CxPauliTensor>() == yiy);
  SpCxPauliTensor xb(Qubit("b", {1, 0}), Pauli::X, -2.3);
  j = xb;
  CHECK(j.get<SpCxPauliTensor>() == xb);
  SymPauliTensor izyx({Pauli::I, Pauli::Z, Pauli::Y, Pauli::X}, Expr("g"));
  j = izyx;
  CHECK(j.get<SymPauliTensor>() == izyx);
  SpSymPauliTensor xaxb(
      {Qubit("a", 0), Qubit("b")}, {Pauli::X, Pauli::X}, -1.98);
  j = xaxb;
  CHECK(j.get<SpSymPauliTensor>() == xaxb);
}

SCENARIO("Test matrix evaluation") {
  GIVEN("Default ordering") {
    SpPauliString ixs({Qubit(0), Qubit(1)}, {Pauli::I, Pauli::X});
    PauliString ixd({Pauli::I, Pauli::X});
    CmplxSpMat ix(4, 4);
    ix.insert(0, 1) = 1.;
    ix.insert(1, 0) = 1.;
    ix.insert(2, 3) = 1.;
    ix.insert(3, 2) = 1.;
    // Eigen sparse matrices don't have an equality check, isApprox is the
    // nearest
    CHECK(ixs.to_sparse_matrix().isApprox(ix));
    CHECK(ixd.to_sparse_matrix().isApprox(ix));
    SpPauliString ixq({{Qubit("b", 0), Pauli::X}, {Qubit("a", 1), Pauli::I}});
    CHECK(ixq.to_sparse_matrix().isApprox(ix));
  }
  GIVEN("Padding to n qubits") {
    SpPauliString ixs(DensePauliMap{Pauli::I, Pauli::X});
    PauliString ixd({Pauli::I, Pauli::X});
    CmplxSpMat ixi(8, 8);
    ixi.insert(0, 2) = 1.;
    ixi.insert(1, 3) = 1.;
    ixi.insert(2, 0) = 1.;
    ixi.insert(3, 1) = 1.;
    ixi.insert(4, 6) = 1.;
    ixi.insert(5, 7) = 1.;
    ixi.insert(6, 4) = 1.;
    ixi.insert(7, 5) = 1.;
    CHECK(ixs.to_sparse_matrix(3).isApprox(ixi));
    CHECK(ixd.to_sparse_matrix(3).isApprox(ixi));
  }
  GIVEN("Custom qubit ordering") {
    SpPauliString ixs(DensePauliMap{Pauli::I, Pauli::X});
    PauliString ixd({Pauli::I, Pauli::X});
    CmplxSpMat xi(4, 4);
    xi.insert(0, 2) = 1.;
    xi.insert(1, 3) = 1.;
    xi.insert(2, 0) = 1.;
    xi.insert(3, 1) = 1.;
    CmplxSpMat ixi(8, 8);
    ixi.insert(0, 2) = 1.;
    ixi.insert(1, 3) = 1.;
    ixi.insert(2, 0) = 1.;
    ixi.insert(3, 1) = 1.;
    ixi.insert(4, 6) = 1.;
    ixi.insert(5, 7) = 1.;
    ixi.insert(6, 4) = 1.;
    ixi.insert(7, 5) = 1.;
    CHECK(ixs.to_sparse_matrix({Qubit(1), Qubit(0)}).isApprox(xi));
    CHECK(ixs.to_sparse_matrix({Qubit(2), Qubit(1), Qubit(0)}).isApprox(ixi));
    CHECK(ixd.to_sparse_matrix({Qubit(1), Qubit(0)}).isApprox(xi));
    CHECK(ixd.to_sparse_matrix({Qubit(2), Qubit(1), Qubit(0)}).isApprox(ixi));
  }
  GIVEN("Different strings") {
    PauliString xyzd({Pauli::X, Pauli::Y, Pauli::Z});
    CmplxSpMat xyz(8, 8);
    xyz.insert(0, 6) = -i_;
    xyz.insert(1, 7) = i_;
    xyz.insert(2, 4) = i_;
    xyz.insert(3, 5) = -i_;
    xyz.insert(4, 2) = -i_;
    xyz.insert(5, 3) = i_;
    xyz.insert(6, 0) = i_;
    xyz.insert(7, 1) = -i_;
    CHECK(xyzd.to_sparse_matrix().isApprox(xyz));
  }
  GIVEN("Different coefficients") {
    CmplxSpMat ix(4, 4);
    ix.insert(0, 1) = 1.;
    ix.insert(1, 0) = 1.;
    ix.insert(2, 3) = 1.;
    ix.insert(3, 2) = 1.;
    DensePauliMap ixd{Pauli::I, Pauli::X};
    CHECK(PauliString(ixd).to_sparse_matrix().isApprox(ix));
    CHECK(PauliStabiliser(ixd, 0).to_sparse_matrix().isApprox(ix));
    CHECK(PauliStabiliser(ixd, 1).to_sparse_matrix().isApprox(i_ * ix));
    CHECK(PauliStabiliser(ixd, 2).to_sparse_matrix().isApprox(-ix));
    CHECK(PauliStabiliser(ixd, 3).to_sparse_matrix().isApprox(-i_ * ix));
    CHECK(CxPauliTensor(ixd, 4.2 + 0.1 * i_)
              .to_sparse_matrix()
              .isApprox(Complex(4.2 + 0.1 * i_) * ix));
    CHECK(SymPauliTensor(ixd, 4.2 + 0.1 * i_)
              .to_sparse_matrix()
              .isApprox(Complex(4.2 + 0.1 * i_) * ix));
  }
}

}  // namespace test_PauliTensor
}  // namespace tket
