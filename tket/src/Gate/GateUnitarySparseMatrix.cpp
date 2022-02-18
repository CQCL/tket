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

#include "GateUnitarySparseMatrix.hpp"

#include "Gate/Gate.hpp"
#include "GateUnitaryMatrix.hpp"
#include "GateUnitaryMatrixError.hpp"
#include "GateUnitaryMatrixImplementations.hpp"
#include "GateUnitaryMatrixUtils.hpp"
#include "Utils/Assert.hpp"

namespace tket {
namespace internal {

// Given a controlled type, which can be written as
// a number of control qubits applied to a more primitive type,
// return the primitive type. (e.g. CnX returns X).
// Returns "noop" if it's not a controlled type we know how to deal with.
// We are excluding things like, e.g., CX, which are only 4x4 and too small
// to be worth dealing with specially.
static OpType get_primitive_type(OpType type_without_controls) {
  switch (type_without_controls) {
    case OpType::CnX:
      // fall through
    case OpType::CCX:
      return OpType::X;

    case OpType::CnRy:
      return OpType::Ry;

    default:
      return OpType::noop;
  }
}

// We have a type acting on 1 qubit.
// Convert this to a type acting on n qubits (given by the Gate object)
// by adding controls.
static std::vector<TripletCd>
convert_1qb_type_to_controlled_type_and_get_triplets(
    const Gate& gate, OpType one_qubit_type, double abs_epsilon) {
  const Gate new_gate(one_qubit_type, gate.get_params(), 1);
  const auto small_unitary = GateUnitaryMatrix::get_unitary(new_gate);
  auto triplets = get_triplets(small_unitary, abs_epsilon);

  // e.g., if CnX or CnRy for n=3, then U is 2x2, but we are embedding into
  // the bottom right corner of an 8x8 identity matrix
  unsigned full_matr_size = get_matrix_size(gate.n_qubits());
  unsigned translation = full_matr_size - small_unitary.rows();
  if (translation > 0) {
    for (auto& triplet : triplets) {
      triplet = TripletCd(
          triplet.row() + translation, triplet.col() + translation,
          triplet.value());
    }
    triplets.reserve(triplets.size() + translation);
    for (unsigned ii = 0; ii < translation; ++ii) {
      triplets.emplace_back(ii, ii, 1.0);
    }
  }
  return triplets;
}

namespace {
struct FixedTripletsWithNoParameters {
  std::vector<TripletCd> bridge_triplets;
  std::vector<TripletCd> cswap_triplets;

  FixedTripletsWithNoParameters();

  // It just so happens that these gates all take
  // the same number of qubits and parameters.
  // Therefore, test the gate for this, as well as
  // returning the const data for use.
  static const FixedTripletsWithNoParameters& get(const Gate& gate);
};

FixedTripletsWithNoParameters::FixedTripletsWithNoParameters() {
  bridge_triplets.reserve(8);
  cswap_triplets.reserve(8);
  const auto& bridge_cols =
      GateUnitaryMatrixImplementations::get_bridge_columns();
  const auto& cswap_cols =
      GateUnitaryMatrixImplementations::get_cswap_columns();
  for (unsigned ii = 0; ii < 8; ++ii) {
    bridge_triplets.emplace_back(ii, bridge_cols[ii], 1.0);
    cswap_triplets.emplace_back(ii, cswap_cols[ii], 1.0);
  }
}

const FixedTripletsWithNoParameters& FixedTripletsWithNoParameters::get(
    const Gate& gate) {
  static const FixedTripletsWithNoParameters data;

  GateUnitaryMatrixUtils::check_and_throw_upon_wrong_number_of_parameters(
      gate.get_type(), gate.n_qubits(),
      GateUnitaryMatrixUtils::get_checked_parameters(gate), 0);
  TKET_ASSERT(gate.n_qubits() == 3);
  return data;
}
}  // namespace

static std::vector<TripletCd> get_phase_gadget_triplets(
    unsigned number_of_qubits, double param) {
  // All diagonal entries have abs value 1,
  // so no point in abs_epsilon.
  const auto entries =
      GateUnitaryMatrixImplementations::PhaseGadget_diagonal_entries(
          number_of_qubits, param);
  std::vector<TripletCd> triplets(entries.rows());
  for (unsigned ii = 0; ii < entries.rows(); ++ii) {
    triplets[ii] = TripletCd(ii, ii, entries(ii));
  }
  return triplets;
}

static std::vector<TripletCd> get_triplets_for_noncontrolled_gate(
    const Gate& gate) {
  switch (gate.get_type()) {
    case OpType::CSWAP:
      return FixedTripletsWithNoParameters::get(gate).cswap_triplets;
    case OpType::BRIDGE:
      return FixedTripletsWithNoParameters::get(gate).bridge_triplets;
    case OpType::PhaseGadget: {
      const auto params = GateUnitaryMatrixUtils::get_checked_parameters(gate);
      GateUnitaryMatrixUtils::check_and_throw_upon_wrong_number_of_parameters(
          gate.get_type(), gate.n_qubits(), params, 1);
      return get_phase_gadget_triplets(gate.n_qubits(), params[0]);
    }
    default:
      return {};
  }
}

std::vector<TripletCd> GateUnitarySparseMatrix::get_unitary_triplets(
    const Gate& gate, double abs_epsilon) {
  const auto type = gate.get_type();
  const auto primitive_type = get_primitive_type(type);
  if (primitive_type != OpType::noop) {
    try {
      return convert_1qb_type_to_controlled_type_and_get_triplets(
          gate, primitive_type, abs_epsilon);
    } catch (const GateUnitaryMatrixError& e) {
      // GCOVR_EXCL_START
      TKET_ASSERT(
          AssertMessage()
          << "Converting " << gate.get_name()
          << " to sparse unitary, via adding controls to gate type "
          << OpDesc(primitive_type).name() << ": " << e.what());
      // GCOVR_EXCL_STOP
    }
  }
  return get_triplets_for_noncontrolled_gate(gate);
}

}  // namespace internal
}  // namespace tket
