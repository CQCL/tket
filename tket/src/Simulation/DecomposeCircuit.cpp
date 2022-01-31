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

#include "DecomposeCircuit.hpp"

#include <sstream>

#include "Circuit/Boxes.hpp"
#include "Circuit/Circuit.hpp"
#include "Gate/Gate.hpp"
#include "Gate/GateUnitaryMatrix.hpp"
#include "GateNodesBuffer.hpp"
#include "PauliExpBoxUnitaryCalculator.hpp"
#include "Utils/Assert.hpp"
#include "Utils/Exceptions.hpp"

namespace tket {
namespace tket_sim {
namespace internal {

typedef std::map<Qubit, unsigned> QMap;

// No checks because only used internally
static QMap get_qmap_no_checks(
    const Circuit& circ,
    const std::vector<unsigned>& parent_circuit_qubit_indices) {
  std::map<Qubit, unsigned> qmap;
  unsigned index = 0;
  for (const Qubit& q : circ.all_qubits()) {
    TKET_ASSERT(qmap.count(q) == 0);
    qmap[q] = parent_circuit_qubit_indices[index];
    ++index;
  }
  TKET_ASSERT(qmap.size() <= parent_circuit_qubit_indices.size());
  return qmap;
}

// Already known to be a box, with a nonempty Op ptr.
// If possible (if the op is able to calculate its own unitary matrix),
// fill the node triplets with the raw unitary matrix represented by this box.
// Return true if it found the triplets.
static bool fill_triplets_directly_from_box(
    GateNode& node, const std::shared_ptr<const Box>& box_ptr, OpType type,
    double abs_epsilon) {
  switch (type) {
    case OpType::Unitary1qBox: {
      auto u1q_ptr = dynamic_cast<const Unitary1qBox*>(box_ptr.get());
      TKET_ASSERT(u1q_ptr);
      node.triplets = tket::get_triplets(u1q_ptr->get_matrix(), abs_epsilon);
      return true;
    }
    case OpType::Unitary2qBox: {
      auto u2q_ptr = dynamic_cast<const Unitary2qBox*>(box_ptr.get());
      TKET_ASSERT(u2q_ptr);
      node.triplets = tket::get_triplets(u2q_ptr->get_matrix(), abs_epsilon);
      return true;
    }
    case OpType::Unitary3qBox: {
      auto u3q_ptr = dynamic_cast<const Unitary3qBox*>(box_ptr.get());
      TKET_ASSERT(u3q_ptr);
      node.triplets = tket::get_triplets(u3q_ptr->get_matrix(), abs_epsilon);
      return true;
    }
    case OpType::PauliExpBox: {
      auto pauli_box_ptr = dynamic_cast<const PauliExpBox*>(box_ptr.get());
      TKET_ASSERT(pauli_box_ptr);
      node.triplets = get_triplets(*pauli_box_ptr);
      return true;
    }
    case OpType::ExpBox: {
      auto exp_box_ptr = dynamic_cast<const ExpBox*>(box_ptr.get());
      TKET_ASSERT(exp_box_ptr);
      const auto numerical_data = exp_box_ptr->get_matrix_and_phase();
      Eigen::Matrix4cd matrix =
          numerical_data.first * (i_ * numerical_data.second);
      matrix = matrix.exp();
      node.triplets = tket::get_triplets(matrix);
      return true;
    }
    default:
      return false;
  }
}

static void add_global_phase(const Circuit& circ, GateNodesBuffer& buffer) {
  const auto global_phase = eval_expr(circ.get_phase());
  if (!global_phase) {
    throw NotImplemented("Circuit has symbolic global phase");
  }
  const double phase_value = global_phase.value();
  buffer.add_global_phase(phase_value);
}

static void throw_with_op_error_message(
    const std::string& op_name, const QMap& qmap, const Circuit& circ,
    const std::string& extra_message) {
  std::stringstream ss;
  ss << "Subcircuit\n"
     << circ << "\nwith " << qmap.size() << " qubits, has op " << op_name
     << ". " << extra_message;
  throw NotImplemented(ss.str());
}

static void fill_qubit_indices(
    const unit_vector_t& args, const QMap& qmap, GateNode& node) {
  TKET_ASSERT(args.size() <= qmap.size());
  node.qubit_indices.resize(args.size());
  for (unsigned ii = 0; ii < args.size(); ++ii) {
    node.qubit_indices[ii] = qmap.at(Qubit(args[ii]));
  }
}

static void decompose_circuit_recursive(
    const Circuit& circ, GateNodesBuffer& buffer,
    const std::vector<unsigned>& parent_circuit_qubit_indices,
    double abs_epsilon) {
  const auto qmap = get_qmap_no_checks(circ, parent_circuit_qubit_indices);
  unit_vector_t args;
  GateNode node;

  for (const Command& command : circ) {
    const Op_ptr current_op = command.get_op_ptr();
    TKET_ASSERT(current_op.get());
    {
      const auto current_type = current_op->get_type();
      if (current_type == OpType::noop || current_type == OpType::Barrier ||
          current_type == OpType::Measure) {
        // Really, ignoring OpType::Measure is the wrong thing to do,
        // since get_unitary is completely meaningless for such circuits.
        // But this was the old behaviour in pytket tests with other
        // simulators, e.g. QuEST, purely for testing purposes.
        // Other OpTypes which are inherently meaningless
        // with get_unitary will throw automatically,
        // simply because they are unknown types.
        continue;
      }
    }
    const OpDesc desc = current_op->get_desc();
    args = command.get_args();
    if (args.empty()) {
      throw_with_op_error_message(desc.name(), qmap, circ, "No args!");
    }
    fill_qubit_indices(args, qmap, node);
    if (desc.is_gate()) {
      const Gate* gate = dynamic_cast<const Gate*>(current_op.get());
      TKET_ASSERT(gate);
      node.triplets =
          GateUnitaryMatrix::get_unitary_triplets(*gate, abs_epsilon);
      buffer.push(node);
      continue;
    }
    if (!desc.is_box()) {
      throw_with_op_error_message(
          desc.name(), qmap, circ, "This is not a gate or box type.");
    }
    std::shared_ptr<const Box> box_ptr =
        std::dynamic_pointer_cast<const Box>(current_op);
    TKET_ASSERT(box_ptr.get());

    node.triplets.clear();
    const bool found_unitary = fill_triplets_directly_from_box(
        node, box_ptr, current_op->get_type(), abs_epsilon);

    if (!node.triplets.empty()) {
      TKET_ASSERT(found_unitary);
      buffer.push(node);
      continue;
    }
    TKET_ASSERT(!found_unitary);
    // Break this box down, recursively.
    std::shared_ptr<Circuit> box_circ = box_ptr->to_circuit();
    if (!box_circ.get()) {
      throw_with_op_error_message(
          desc.name(), qmap, circ,
          "This is a box, which couldn't be "
          "broken down into a circuit");
    }
    decompose_circuit_recursive(
        *box_circ, buffer, node.qubit_indices, abs_epsilon);
  }
  add_global_phase(circ, buffer);
}

void decompose_circuit(
    const Circuit& circ, GateNodesBuffer& buffer, double abs_epsilon) {
  // The qubits are just [0,1,2,...].
  std::vector<unsigned> iota(circ.n_qubits());
  std::iota(iota.begin(), iota.end(), 0);

  decompose_circuit_recursive(circ, buffer, iota, abs_epsilon);
  buffer.flush();
}

}  // namespace internal
}  // namespace tket_sim
}  // namespace tket
