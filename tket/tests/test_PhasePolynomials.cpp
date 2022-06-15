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

#include "Circuit/Boxes.hpp"
#include "Circuit/CircUtils.hpp"
#include "CircuitsForTesting.hpp"
#include "Converters/PhasePoly.hpp"
#include "Predicates/CompilerPass.hpp"
#include "Transformations/Decomposition.hpp"
#include "Transformations/Rebase.hpp"
#include "Transformations/Transform.hpp"
#include "Utils/HelperFunctions.hpp"
#include "Utils/MatrixAnalysis.hpp"
#include "testutil.hpp"

namespace tket {
namespace test_PhasePolynomials {

/**
 * generate circuit for tests with large structures of Rz+CX gates
 * @param n size of the circuit
 * @return circuit
 */
Circuit generate_test_circuit(unsigned n) {
  Circuit circ(n);
  if (n >= 2) {
    for (unsigned i = 0; i < n; ++i) {
      circ.add_op<unsigned>(OpType::H, {i});
    }

    circ.add_op<unsigned>(OpType::Rz, 0.01, {0});
    for (unsigned i = 1; i < n; ++i) {
      circ.add_op<unsigned>(OpType::CX, {i, 0});
      circ.add_op<unsigned>(OpType::Rz, 0.1 * i, {0});
    }

    for (unsigned i = 0; i < 2; ++i) {
      circ.add_op<unsigned>(OpType::H, {i});
    }
    circ.add_barrier({0, 1});

    circ.add_op<unsigned>(OpType::Rz, 0.01, {0});
    for (unsigned i = 1; i < n; ++i) {
      circ.add_op<unsigned>(OpType::CX, {i, 0});
      circ.add_op<unsigned>(OpType::Rz, 0.1 * i, {0});
    }

    for (unsigned i = 0; i < n; ++i) {
      circ.add_op<unsigned>(OpType::H, {i});
    }
    if (n > 3) {
      circ.add_op<unsigned>(OpType::CX, {0, 1});
      circ.add_op<unsigned>(OpType::Rz, 0.1 * 2, {0});
      circ.add_op<unsigned>(OpType::Rz, 0.1 * 3, {1});
    }
  }
  return circ;
}

SCENARIO("Test basic phase polynomial creation") {
  GIVEN("A 2qb CX+Rz phase gadget") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.3, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    PhasePolyBox ppbox(circ);

    const PhasePolynomial& phasepoly = ppbox.get_phase_polynomial();
    REQUIRE(phasepoly.size() == 1);
    std::vector<bool> correct{1, 1};
    REQUIRE(phasepoly.begin()->first == correct);
    REQUIRE(test_equiv_val(phasepoly.begin()->second, 0.3));
    const MatrixXb& basis_map = ppbox.get_linear_transformation();
    REQUIRE(basis_map == MatrixXb::Identity(2, 2));
  }
  GIVEN("A larger CX+Rz phase polynomial") {
    unsigned n = 5;
    Circuit circ(n);
    circ.add_op<unsigned>(OpType::Rz, 0.1, {0});
    for (unsigned i = 1; i < n; ++i) {
      circ.add_op<unsigned>(OpType::CX, {i, 0});
      circ.add_op<unsigned>(OpType::Rz, 0.1, {0});
    }
    PhasePolyBox ppbox(circ);
    const PhasePolynomial& phasepoly = ppbox.get_phase_polynomial();
    REQUIRE(phasepoly.size() == n);
    unsigned i = 1;
    for (auto pp_it = phasepoly.begin(); pp_it != phasepoly.end(); ++pp_it) {
      std::vector<bool> correct(n, false);
      for (unsigned j = 0; j < i; ++j) correct[j] = true;
      REQUIRE(pp_it->first == correct);
      REQUIRE(test_equiv_val(pp_it->second, 0.1));
      ++i;
    }
    const MatrixXb& basis_map = ppbox.get_linear_transformation();
    MatrixXb correct_basis_map = MatrixXb::Identity(n, n);
    for (unsigned i = 0; i < n; ++i) {
      correct_basis_map.row(0)[i] = 1;
    }
    REQUIRE(basis_map == correct_basis_map);
  }
}

SCENARIO("Test diagonal Phase Polynomial circuit generation") {
  GIVEN("A simple Phase Polynomial") {
    const auto& prepend = CircuitsForTesting::get().prepend_2qb_circuit;
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.89, {1});
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    Circuit after1 = prepend >> circ;
    PhasePolyBox ppbox(circ);
    std::shared_ptr<Circuit> circptr = ppbox.to_circuit();
    Circuit circ2 = *circptr;
    Circuit after2 = prepend >> circ2;
    REQUIRE(test_statevector_comparison(after1, after2));
  }
  GIVEN("A more complicated Phase Polynomial") {
    auto prepend = CircuitsForTesting::get_prepend_circuit(5);
    prepend.add_op<unsigned>(OpType::Rx, 0.41, {2});
    prepend.add_op<unsigned>(OpType::Rz, 1.1, {2});
    prepend.add_op<unsigned>(OpType::Rx, 0.3, {3});
    prepend.add_op<unsigned>(OpType::Rz, 1.9, {3});

    Circuit circ(5);
    circ.add_op<unsigned>(OpType::Rz, 1.346, {1});
    add_2qb_gates(circ, OpType::CX, {{0, 1}, {1, 2}, {2, 3}});
    circ.add_op<unsigned>(OpType::Rz, 0.76, {3});
    add_2qb_gates(circ, OpType::CX, {{2, 3}, {1, 2}, {0, 1}});
    circ.add_op<unsigned>(OpType::Rz, 0.76, {0});
    add_2qb_gates(circ, OpType::CX, {{3, 2}, {0, 2}, {2, 4}});
    circ.add_op<unsigned>(OpType::Rz, 1.346, {4});
    add_2qb_gates(circ, OpType::CX, {{2, 4}, {0, 2}, {3, 2}});
    Circuit after1 = prepend >> circ;
    PhasePolyBox ppbox(circ);
    std::shared_ptr<Circuit> circptr = ppbox.to_circuit();
    Circuit circ2 = *circptr;
    Circuit after2 = prepend >> circ2;
    REQUIRE(test_statevector_comparison(after1, after2));
  }
}

SCENARIO("Test affine Phase Polynomial circuit generation") {
  GIVEN("A simple Phase Polynomial") {
    const auto& prepend = CircuitsForTesting::get().prepend_2qb_circuit;
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    circ.add_op<unsigned>(OpType::Rz, 0.32, {1});

    Circuit after1 = prepend >> circ;
    PhasePolyBox ppbox(circ);
    std::shared_ptr<Circuit> circptr = ppbox.to_circuit();
    Circuit circ2 = *circptr;
    Circuit after2 = prepend >> circ2;
    REQUIRE(test_statevector_comparison(after1, after2));
  }
  GIVEN("A 3qb Phase Polynomial") {
    auto prepend = CircuitsForTesting::get_prepend_circuit(3);
    prepend.add_op<unsigned>(OpType::Rx, 0.41, {2});
    prepend.add_op<unsigned>(OpType::Rz, 1.1, {2});

    Circuit circ(3);
    circ.add_op<unsigned>(OpType::Rz, 1.346, {1});
    add_2qb_gates(circ, OpType::CX, {{0, 1}, {1, 2}});
    circ.add_op<unsigned>(OpType::Rz, 0.76, {0});
    circ.add_op<unsigned>(OpType::CX, {0, 2});
    Circuit after1 = prepend >> circ;
    PhasePolyBox ppbox(circ);
    std::shared_ptr<Circuit> circptr = ppbox.to_circuit();
    Circuit circ2 = *circptr;
    Circuit after2 = prepend >> circ2;
    REQUIRE(test_statevector_comparison(after1, after2));
  }
  GIVEN("A 5qb Phase Polynomial") {
    auto prepend = CircuitsForTesting::get_prepend_circuit(5);
    prepend.add_op<unsigned>(OpType::Rx, 0.41, {2});
    prepend.add_op<unsigned>(OpType::Rz, 1.1, {2});
    prepend.add_op<unsigned>(OpType::Rx, 0.3, {3});
    prepend.add_op<unsigned>(OpType::Rz, 1.9, {3});
    prepend.add_op<unsigned>(OpType::Rx, 0.52, {4});
    prepend.add_op<unsigned>(OpType::Rz, 1.5, {4});

    Circuit circ(5);
    circ.add_op<unsigned>(OpType::Rz, 1.346, {1});
    add_2qb_gates(circ, OpType::CX, {{0, 1}, {1, 2}, {2, 3}});
    circ.add_op<unsigned>(OpType::Rz, 0.76, {3});
    circ.add_op<unsigned>(OpType::Rz, 0.76, {0});
    add_2qb_gates(circ, OpType::CX, {{3, 2}, {0, 2}, {2, 4}});
    circ.add_op<unsigned>(OpType::Rz, 1.346, {4});
    Circuit after1 = prepend >> circ;
    PhasePolyBox ppbox(circ);
    std::shared_ptr<Circuit> circptr = ppbox.to_circuit();
    Circuit circ2 = *circptr;
    Circuit after2 = prepend >> circ2;
    REQUIRE(test_statevector_comparison(after1, after2));
  }
}

SCENARIO("Test assertion in PhasePolyBox creation") {
  GIVEN("check invalid qubit_indices i") {
    unsigned n_qubits = 2;

    boost::bimap<Qubit, unsigned> qubit_indices;
    for (unsigned i = 0; i < 3; ++i) {
      qubit_indices.insert({Qubit(i), i});
    }

    PhasePolynomial phase_polynomial = {{{true, false}, 1.0}};

    MatrixXb linear_transformation(2, 2);
    linear_transformation << 0, 1,  //
        1, 0;

    REQUIRE_THROWS_AS(
        PhasePolyBox(
            n_qubits, qubit_indices, phase_polynomial, linear_transformation),
        std::invalid_argument);
  }
  GIVEN("check invalid qubit_indices ii") {
    unsigned n_qubits = 2;

    boost::bimap<Qubit, unsigned> qubit_indices;
    for (unsigned i = 0; i < 2; ++i) {
      qubit_indices.insert({Qubit(i), (i + 1)});
    }

    PhasePolynomial phase_polynomial = {{{true, false}, 1.0}};

    MatrixXb linear_transformation(2, 2);
    linear_transformation << 0, 1,  //
        1, 0;

    REQUIRE_THROWS_AS(
        PhasePolyBox(
            n_qubits, qubit_indices, phase_polynomial, linear_transformation),
        std::invalid_argument);
  }
  GIVEN("check invalid phase_polynomial i") {
    unsigned n_qubits = 2;

    boost::bimap<Qubit, unsigned> qubit_indices;
    for (unsigned i = 0; i < 2; ++i) {
      qubit_indices.insert({Qubit(i), i});
    }

    PhasePolynomial phase_polynomial = {{{true}, 1.0}};

    MatrixXb linear_transformation(2, 2);
    linear_transformation << 0, 1,  //
        1, 0;

    REQUIRE_THROWS_AS(
        PhasePolyBox(
            n_qubits, qubit_indices, phase_polynomial, linear_transformation),
        std::invalid_argument);
  }
  GIVEN("check invalid phase_polynomial ii") {
    unsigned n_qubits = 2;

    boost::bimap<Qubit, unsigned> qubit_indices;
    for (unsigned i = 0; i < 2; ++i) {
      qubit_indices.insert({Qubit(i), i});
    }

    PhasePolynomial phase_polynomial = {{{false, false}, 1.0}};

    MatrixXb linear_transformation(2, 2);
    linear_transformation << 0, 1,  //
        1, 0;

    REQUIRE_THROWS_AS(
        PhasePolyBox(
            n_qubits, qubit_indices, phase_polynomial, linear_transformation),
        std::invalid_argument);
  }
  GIVEN("check invalid linear_transformation i") {
    unsigned n_qubits = 2;

    boost::bimap<Qubit, unsigned> qubit_indices;
    for (unsigned i = 0; i < 2; ++i) {
      qubit_indices.insert({Qubit(i), i});
    }

    PhasePolynomial phase_polynomial = {{{false, false}, 1.0}};

    MatrixXb linear_transformation(1, 1);
    linear_transformation << 1;

    REQUIRE_THROWS_AS(
        PhasePolyBox(
            n_qubits, qubit_indices, phase_polynomial, linear_transformation),
        std::invalid_argument);
  }
}

SCENARIO("Test conversion of circuit to circuit with phase poly boxes") {
  GIVEN("convert_to_phase_poly simple") {
    unsigned n = 3;
    Circuit circ(n);

    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::H, {2});

    circ.add_op<unsigned>(OpType::Rz, 0.01, {0});
    for (unsigned i = 1; i < n; ++i) {
      circ.add_op<unsigned>(OpType::CX, {i, 0});
      circ.add_op<unsigned>(OpType::Rz, 0.1 * i, {0});
    }

    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::H, {2});

    CircToPhasePolyConversion conv = CircToPhasePolyConversion(circ);
    conv.convert();
    Circuit result = conv.get_circuit();

    REQUIRE(test_statevector_comparison(circ, result));
  }
  GIVEN("convert_to_phase_poly complex") {
    unsigned n = 11;
    Circuit circ(n);

    for (unsigned i = 0; i < n; ++i) {
      circ.add_op<unsigned>(OpType::H, {i});
    }

    circ.add_op<unsigned>(OpType::Rz, 0.01, {0});
    for (unsigned i = 1; i < n; ++i) {
      circ.add_op<unsigned>(OpType::CX, {i, 0});
      circ.add_op<unsigned>(OpType::Rz, 0.1 * i, {0});
    }

    for (unsigned i = 0; i < n; ++i) {
      circ.add_op<unsigned>(OpType::H, {i});
    }

    CircToPhasePolyConversion conv = CircToPhasePolyConversion(circ);
    conv.convert();
    Circuit result = conv.get_circuit();

    REQUIRE(test_statevector_comparison(circ, result));
  }
  GIVEN("convert_to_phase_poly complex II") {
    unsigned n = 11;
    Circuit circ(n);

    for (unsigned i = 0; i < n; ++i) {
      circ.add_op<unsigned>(OpType::H, {i});
    }

    circ.add_op<unsigned>(OpType::Rz, 0.01, {0});
    for (unsigned i = 1; i < n; ++i) {
      circ.add_op<unsigned>(OpType::CX, {i, 0});
      circ.add_op<unsigned>(OpType::Rz, 0.1 * i, {0});
    }

    for (unsigned i = 0; i < n; ++i) {
      circ.add_op<unsigned>(OpType::H, {i});
    }

    circ.add_op<unsigned>(OpType::Rz, 0.01, {0});
    for (unsigned i = 1; i < n / 2; ++i) {
      circ.add_op<unsigned>(OpType::CX, {i, 0});
      circ.add_op<unsigned>(OpType::Rz, 0.1 * i, {0});
    }

    for (unsigned i = 0; i < n; ++i) {
      circ.add_op<unsigned>(OpType::H, {i});
    }

    CircToPhasePolyConversion conv = CircToPhasePolyConversion(circ);
    conv.convert();
    Circuit result = conv.get_circuit();

    REQUIRE(test_statevector_comparison(circ, result));
  }
  GIVEN("convert_to_phase_poly simple II ") {
    unsigned n = 2;
    Circuit circ(n);

    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});

    circ.add_op<unsigned>(OpType::Rz, 0.01, {0});
    for (unsigned i = 1; i < n; ++i) {
      circ.add_op<unsigned>(OpType::CX, {i, 0});
      circ.add_op<unsigned>(OpType::Rz, 0.1 * i, {0});
    }

    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});

    CircToPhasePolyConversion conv = CircToPhasePolyConversion(circ);
    conv.convert();
    Circuit result = conv.get_circuit();

    std::vector<OpType> qubit_types = std::vector<OpType>(5, OpType::H);

    qubit_types[2] = OpType::PhasePolyBox;

    int count = 0;
    for (const Command& com : result) {
      OpType ot = com.get_op_ptr()->get_type();
      REQUIRE(ot == qubit_types[count]);
      ++count;
    }

    REQUIRE(test_statevector_comparison(circ, result));
  }
  GIVEN("convert_to_phase_poly simple III ") {
    unsigned n = 2;
    Circuit circ(n);

    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});

    circ.add_op<unsigned>(OpType::Rz, 0.01, {0});
    for (unsigned i = 1; i < n; ++i) {
      circ.add_op<unsigned>(OpType::CX, {i, 0});
      circ.add_op<unsigned>(OpType::Rz, 0.1 * i, {0});
    }

    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});

    circ.add_op<unsigned>(OpType::Rz, 0.01, {0});
    for (unsigned i = 1; i < n; ++i) {
      circ.add_op<unsigned>(OpType::CX, {i, 0});
      circ.add_op<unsigned>(OpType::Rz, 0.1 * i, {0});
    }

    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});

    CircToPhasePolyConversion conv = CircToPhasePolyConversion(circ);
    conv.convert();
    Circuit result = conv.get_circuit();

    std::vector<OpType> qubit_types = std::vector<OpType>(8, OpType::H);

    qubit_types[2] = OpType::PhasePolyBox;
    qubit_types[5] = OpType::PhasePolyBox;

    int count = 0;

    for (const Command& com : result) {
      OpType ot = com.get_op_ptr()->get_type();
      REQUIRE(ot == qubit_types[count]);
      ++count;
    }

    REQUIRE(test_statevector_comparison(circ, result));
  }

  GIVEN("convert_to_phase_poly simple IV ") {
    unsigned n = 2;
    Circuit circ(n);

    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});

    circ.add_op<unsigned>(OpType::Rz, 0.01, {0});
    for (unsigned i = 1; i < n; ++i) {
      circ.add_op<unsigned>(OpType::CX, {i, 0});
      circ.add_op<unsigned>(OpType::Rz, 0.1 * i, {0});
    }

    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_barrier({0, 1});

    circ.add_op<unsigned>(OpType::Rz, 0.01, {0});
    for (unsigned i = 1; i < n; ++i) {
      circ.add_op<unsigned>(OpType::CX, {i, 0});
      circ.add_op<unsigned>(OpType::Rz, 0.1 * i, {0});
    }

    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});

    CircToPhasePolyConversion conv = CircToPhasePolyConversion(circ, 1);
    conv.convert();
    Circuit result = conv.get_circuit();

    std::vector<OpType> qubit_types = std::vector<OpType>(9, OpType::H);

    qubit_types[2] = OpType::PhasePolyBox;
    qubit_types[5] = OpType::Barrier;
    qubit_types[6] = OpType::PhasePolyBox;

    int count = 0;

    for (const Command& com : result) {
      OpType ot = com.get_op_ptr()->get_type();
      REQUIRE(ot == qubit_types[count]);
      ++count;
    }
  }

  GIVEN("convert_to_phase_poly simple V ") {
    unsigned n = 2;
    Circuit circ(n, n);

    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});

    circ.add_op<unsigned>(OpType::Rz, 0.01, {0});
    for (unsigned i = 1; i < n; ++i) {
      circ.add_op<unsigned>(OpType::CX, {i, 0});
      circ.add_op<unsigned>(OpType::Rz, 0.1 * i, {0});
    }

    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::Measure, {0, 0});

    circ.add_op<unsigned>(OpType::Rz, 0.01, {0});
    for (unsigned i = 1; i < n; ++i) {
      circ.add_op<unsigned>(OpType::CX, {i, 0});
      circ.add_op<unsigned>(OpType::Rz, 0.1 * i, {0});
    }

    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});

    CircToPhasePolyConversion conv = CircToPhasePolyConversion(circ);
    conv.convert();
    Circuit result = conv.get_circuit();

    std::vector<OpType> qubit_types = std::vector<OpType>(9, OpType::H);

    qubit_types[2] = OpType::PhasePolyBox;
    qubit_types[5] = OpType::Measure;
    qubit_types[6] = OpType::PhasePolyBox;

    int count = 0;

    for (const Command& com : result) {
      OpType ot = com.get_op_ptr()->get_type();
      REQUIRE(ot == qubit_types[count]);
      ++count;
    }
  }
  GIVEN("convert_to_phase_poly simple VI ") {
    unsigned n = 2;
    Circuit circ(n, n);

    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});

    circ.add_op<unsigned>(OpType::Rz, 0.01, {0});
    for (unsigned i = 1; i < n; ++i) {
      circ.add_op<unsigned>(OpType::CX, {i, 0});
      circ.add_op<unsigned>(OpType::Rz, 0.1 * i, {0});
    }

    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::Measure, {0, 0});
    circ.add_op<unsigned>(OpType::Measure, {1, 1});

    circ.add_op<unsigned>(OpType::Rz, 0.01, {0});
    for (unsigned i = 1; i < n; ++i) {
      circ.add_op<unsigned>(OpType::CX, {i, 0});
      circ.add_op<unsigned>(OpType::Rz, 0.1 * i, {0});
    }

    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});

    CircToPhasePolyConversion conv = CircToPhasePolyConversion(circ);
    conv.convert();
    Circuit result = conv.get_circuit();

    std::vector<OpType> qubit_types = std::vector<OpType>(10, OpType::H);

    std::map<Qubit, unsigned> qubit_indices;
    std::map<Bit, unsigned> bit_indices;

    unsigned i = 0;
    for (const Qubit& qb : circ.all_qubits()) {
      qubit_indices.insert({qb, i});
      ++i;
    }

    i = 0;
    for (const Bit& b : circ.all_bits()) {
      bit_indices.insert({b, i});
      ++i;
    }

    qubit_types[2] = OpType::PhasePolyBox;
    qubit_types[5] = OpType::Measure;
    qubit_types[6] = OpType::Measure;
    qubit_types[7] = OpType::PhasePolyBox;

    int count = 0;
    int countmeasure = 0;

    std::vector<unsigned> measure_qu_bits = std::vector<unsigned>(2, 0);
    measure_qu_bits[1] = 1;

    for (const Command& com : result) {
      OpType ot = com.get_op_ptr()->get_type();
      unit_vector_t qbs = com.get_args();
      unsigned qb = qubit_indices.at(Qubit(qbs[0]));
      unsigned b = 11;
      if (ot == OpType::Measure) {
        b = bit_indices.at(Bit(qbs[1]));
        REQUIRE(qb == measure_qu_bits[countmeasure]);
        REQUIRE(b == measure_qu_bits[countmeasure]);
        ++countmeasure;
      }

      REQUIRE(ot == qubit_types[count]);
      ++count;
    }
  }
  GIVEN("convert_to_phase_poly simple VII ") {
    unsigned n = 2;
    Circuit circ(n, n);

    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});

    circ.add_op<unsigned>(OpType::Rz, 0.01, {0});
    for (unsigned i = 1; i < n; ++i) {
      circ.add_op<unsigned>(OpType::CX, {i, 0});
      circ.add_op<unsigned>(OpType::Rz, 0.1 * i, {0});
    }

    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::Measure, {0, 0});
    circ.add_op<unsigned>(OpType::Measure, {1, 1});
    circ.add_op<unsigned>(OpType::Reset, {0});
    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});
    circ.add_op<unsigned>(OpType::Collapse, {0});

    circ.add_op<unsigned>(OpType::Rz, 0.01, {0});
    for (unsigned i = 1; i < n; ++i) {
      circ.add_op<unsigned>(OpType::CX, {i, 0});
      circ.add_op<unsigned>(OpType::Rz, 0.1 * i, {0});
    }

    circ.add_op<unsigned>(OpType::H, {0});
    circ.add_op<unsigned>(OpType::H, {1});

    CircToPhasePolyConversion conv = CircToPhasePolyConversion(circ);
    conv.convert();
    Circuit result = conv.get_circuit();

    std::vector<OpType> qubit_types = std::vector<OpType>(14, OpType::H);

    std::map<Qubit, unsigned> qubit_indices;
    std::map<Bit, unsigned> bit_indices;

    unsigned i = 0;
    for (const Qubit& qb : circ.all_qubits()) {
      qubit_indices.insert({qb, i});
      ++i;
    }

    i = 0;
    for (const Bit& b : circ.all_bits()) {
      bit_indices.insert({b, i});
      ++i;
    }

    qubit_types[2] = OpType::PhasePolyBox;
    qubit_types[5] = OpType::Measure;
    qubit_types[6] = OpType::Measure;
    qubit_types[7] = OpType::Reset;
    qubit_types[10] = OpType::Collapse;
    qubit_types[11] = OpType::PhasePolyBox;

    int count = 0;
    int countmeasure = 0;

    std::vector<unsigned> measure_qu_bits = std::vector<unsigned>(2, 0);
    measure_qu_bits[1] = 1;

    for (const Command& com : result) {
      OpType ot = com.get_op_ptr()->get_type();
      unit_vector_t qbs = com.get_args();
      unsigned qb = qubit_indices.at(Qubit(qbs[0]));
      unsigned b = 11;
      if (ot == OpType::Measure) {
        b = bit_indices.at(Bit(qbs[1]));
        REQUIRE(qb == measure_qu_bits[countmeasure]);
        REQUIRE(b == measure_qu_bits[countmeasure]);
        ++countmeasure;
      }

      REQUIRE(ot == qubit_types[count]);
      ++count;
    }
  }
  GIVEN("convert_to_phase_poly minimal box size") {
    Circuit circ = generate_test_circuit(2);

    CircToPhasePolyConversion conv = CircToPhasePolyConversion(circ, 2);
    conv.convert();
    Circuit result = conv.get_circuit();

    std::vector<OpType> qubit_types = std::vector<OpType>(13, OpType::H);

    qubit_types[2] = OpType::Rz;
    qubit_types[3] = OpType::CX;
    qubit_types[4] = OpType::Rz;
    qubit_types[7] = OpType::Barrier;
    qubit_types[8] = OpType::Rz;
    qubit_types[9] = OpType::CX;
    qubit_types[10] = OpType::Rz;

    int count = 0;

    for (const Command& com : result) {
      OpType ot = com.get_op_ptr()->get_type();
      REQUIRE(ot == qubit_types[count]);
      ++count;
    }
  }
  GIVEN("convert_to_phase_poly minimal box size II") {
    Circuit circ = generate_test_circuit(4);

    CircToPhasePolyConversion conv = CircToPhasePolyConversion(circ, 2);
    conv.convert();
    Circuit result = conv.get_circuit();

    std::vector<OpType> qubit_types = std::vector<OpType>(16, OpType::H);

    qubit_types[4] = OpType::PhasePolyBox;
    qubit_types[7] = OpType::Barrier;
    qubit_types[8] = OpType::PhasePolyBox;
    qubit_types[13] = OpType::CX;
    qubit_types[14] = OpType::Rz;
    qubit_types[15] = OpType::Rz;

    int count = 0;

    for (const Command& com : result) {
      OpType ot = com.get_op_ptr()->get_type();
      REQUIRE(ot == qubit_types[count]);
      ++count;
    }
  }
  GIVEN("compilerpass - phase poly box I") {
    unsigned n = 2;
    Circuit circ(n);

    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::X, {1});

    circ.add_op<unsigned>(OpType::Rz, 0.01, {0});
    for (unsigned i = 1; i < n; ++i) {
      circ.add_op<unsigned>(OpType::CX, {i, 0});
      circ.add_op<unsigned>(OpType::Rz, 0.1 * i, {0});
    }

    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::X, {1});
    circ.add_op<unsigned>(OpType::Y, {0});
    circ.add_op<unsigned>(OpType::Y, {1});

    circ.add_op<unsigned>(OpType::Rz, 0.01, {0});
    for (unsigned i = 1; i < n; ++i) {
      circ.add_op<unsigned>(OpType::CX, {i, 0});
      circ.add_op<unsigned>(OpType::Rz, 0.1 * i, {0});
    }

    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::X, {1});

    Circuit result = circ;

    Transform t =
        (Transforms::rebase_UFR() >> Transforms::compose_phase_poly_boxes());

    REQUIRE(t.apply(circ));
    REQUIRE(circ.count_gates(OpType::X) == 0);
    REQUIRE(circ.count_gates(OpType::Y) == 0);
    REQUIRE(circ.count_gates(OpType::CX) == 0);
    REQUIRE(circ.count_gates(OpType::Rz) == 0);
    REQUIRE(circ.count_gates(OpType::H) == 16);
    REQUIRE(circ.count_gates(OpType::PhasePolyBox) == 7);

    REQUIRE(test_statevector_comparison(circ, result));
  }
  GIVEN("compilerpass - phase poly box II") {
    unsigned n = 11;
    Circuit circ(n);

    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::X, {1});
    circ.add_op<unsigned>(OpType::X, {2});
    circ.add_op<unsigned>(OpType::X, {3});
    circ.add_op<unsigned>(OpType::X, {4});
    circ.add_op<unsigned>(OpType::X, {5});
    circ.add_op<unsigned>(OpType::X, {6});

    circ.add_op<unsigned>(OpType::Rz, 0.01, {0});
    for (unsigned i = 1; i < n; ++i) {
      circ.add_op<unsigned>(OpType::CX, {i, 0});
      circ.add_op<unsigned>(OpType::Rz, 0.1 * i, {0});
    }

    for (unsigned i = 1; i < n; ++i) {
      circ.add_op<unsigned>(OpType::X, {i});
      circ.add_op<unsigned>(OpType::Y, {i});
      circ.add_op<unsigned>(OpType::X, {i});
      circ.add_op<unsigned>(OpType::Y, {i});
    }

    circ.add_op<unsigned>(OpType::Rz, 0.01, {0});
    for (unsigned i = 1; i < n; ++i) {
      circ.add_op<unsigned>(OpType::CX, {i, 0});
      circ.add_op<unsigned>(OpType::Rz, 0.1 * i, {0});
    }

    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::X, {1});

    Circuit result = circ;

    Transform t =
        (Transforms::rebase_UFR() >> Transforms::compose_phase_poly_boxes());

    REQUIRE(t.apply(circ));
    REQUIRE(circ.count_gates(OpType::X) == 0);
    REQUIRE(circ.count_gates(OpType::Y) == 0);
    REQUIRE(circ.count_gates(OpType::CX) == 0);
    REQUIRE(circ.count_gates(OpType::Rz) == 0);
    REQUIRE(circ.count_gates(OpType::H) == 98);
    REQUIRE(circ.count_gates(OpType::PhasePolyBox) == 19);

    REQUIRE(test_statevector_comparison(circ, result));
  }
  GIVEN("compilerpass - phase poly box II") {
    unsigned n = 11;
    Circuit circ(n);

    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::X, {1});
    circ.add_op<unsigned>(OpType::X, {2});
    circ.add_op<unsigned>(OpType::X, {3});
    circ.add_op<unsigned>(OpType::X, {4});
    circ.add_op<unsigned>(OpType::X, {5});
    circ.add_op<unsigned>(OpType::X, {6});

    circ.add_op<unsigned>(OpType::Rz, 0.01, {0});
    for (unsigned i = 1; i < n; ++i) {
      circ.add_op<unsigned>(OpType::CX, {i, 0});
      circ.add_op<unsigned>(OpType::Rz, 0.1 * i, {0});
    }

    for (unsigned i = 1; i < n; ++i) {
      circ.add_op<unsigned>(OpType::X, {i});
      circ.add_op<unsigned>(OpType::Y, {i});
      circ.add_op<unsigned>(OpType::X, {i});
      circ.add_op<unsigned>(OpType::Y, {i});
    }

    circ.add_op<unsigned>(OpType::Rz, 0.01, {0});
    for (unsigned i = 1; i < n; ++i) {
      circ.add_op<unsigned>(OpType::CX, {i, 0});
      circ.add_op<unsigned>(OpType::Rz, 0.1 * i, {0});
    }

    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_op<unsigned>(OpType::X, {1});

    Circuit result = circ;

    Transform t =
        (Transforms::rebase_UFR() >> Transforms::compose_phase_poly_boxes());

    REQUIRE(t.apply(circ));
    REQUIRE(circ.count_gates(OpType::X) == 0);
    REQUIRE(circ.count_gates(OpType::Y) == 0);
    REQUIRE(circ.count_gates(OpType::CX) == 0);
    REQUIRE(circ.count_gates(OpType::Rz) == 0);
    REQUIRE(circ.count_gates(OpType::H) == 98);
    REQUIRE(circ.count_gates(OpType::PhasePolyBox) == 19);

    REQUIRE(test_statevector_comparison(circ, result));
  }
}
SCENARIO("Phase polynomial synthesis without architecture") {
  GIVEN("single SWAP circuit") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::SWAP, {0, 1});
    circ.replace_SWAPs();
    PhasePolyBox ppbox(circ);
    std::shared_ptr<Circuit> circptr = ppbox.to_circuit();
    Circuit circ2 = *circptr;
    REQUIRE(test_unitary_comparison(circ, circ2));
  }
  GIVEN("more SWAP circuit") {
    Circuit circ(5);
    circ.add_op<unsigned>(OpType::SWAP, {0, 1});
    circ.add_op<unsigned>(OpType::SWAP, {0, 2});
    circ.add_op<unsigned>(OpType::SWAP, {0, 3});
    circ.replace_SWAPs();
    PhasePolyBox ppbox(circ);
    std::shared_ptr<Circuit> circptr = ppbox.to_circuit();
    Circuit circ2 = *circptr;
    REQUIRE(test_unitary_comparison(circ, circ2));
  }
  GIVEN("more SWAP circuit II") {
    Circuit circ(5);
    circ.add_op<unsigned>(OpType::SWAP, {0, 1});
    circ.add_op<unsigned>(OpType::SWAP, {1, 2});
    circ.add_op<unsigned>(OpType::SWAP, {2, 3});
    circ.add_op<unsigned>(OpType::SWAP, {3, 4});
    circ.replace_SWAPs();
    PhasePolyBox ppbox(circ);
    std::shared_ptr<Circuit> circptr = ppbox.to_circuit();
    Circuit circ2 = *circptr;
    REQUIRE(test_unitary_comparison(circ, circ2));
  }
}

}  // namespace test_PhasePolynomials
}  // namespace tket
