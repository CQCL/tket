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

#include "tket/ArchAwareSynth/SteinerForest.hpp"

#include <algorithm>
#include <vector>

#include "tket/ArchAwareSynth/SteinerTree.hpp"
#include "tket/Architecture/Architecture.hpp"

namespace tket {
namespace aas {

ParityList parity_column_to_list(const std::vector<bool> &parity_column) {
  ParityList parity_list;
  for (unsigned i = 0; i != parity_column.size(); ++i) {
    if (parity_column[i]) {
      parity_list.push_back(i);
    }
  }
  return parity_list;
}

SteinerForest::SteinerForest(
    const PathHandler &paths, const PhasePolyBox &phasepolybox) {
  global_cost = 0;
  unsigned id_counter = 0;
  const PhasePolynomial &phasepoly = phasepolybox.get_phase_polynomial();
  const MatrixXb &linear_fn_matrix_ppb =
      phasepolybox.get_linear_transformation();
  unsigned n_qubits = phasepolybox.get_n_qubits();
  synth_circuit = Circuit(n_qubits);
  linear_function = DiagMatrix(linear_fn_matrix_ppb);

  for (const auto &parity : phasepoly) {
    ParityList parity_list = parity_column_to_list(parity.first);
    SteinerTree tree(paths, parity_list, *parity_list.begin());
    global_cost += tree.tree_cost;
    std::pair<SteinerTree, Expr> parity_on_graph{tree, parity.second};

    CostedTrees::iterator iter = current_trees.find(tree.tree_cost);
    if (iter == current_trees.end())
      current_trees.insert({tree.tree_cost, {parity_on_graph}});
    else
      iter->second.push_back(parity_on_graph);

    ++id_counter;
  }
  tree_count = id_counter;

  CostedTrees new_trees;
  for (const auto &cost_trees : current_trees) {
    for (auto tree_expr : cost_trees.second) {
      /* If we come across a tree which has only 1 node, check that its
      properties are correct
      and remove it from the forest, leaving an Rz in the circuit. */
      if (tree_expr.first.fully_reduced()) {
        unsigned qubit_id = UINT_MAX;
        for (SteinerNodeType q : tree_expr.first.node_types) {
          if (q == SteinerNodeType::Leaf) {
            TKET_ASSERT(qubit_id == UINT_MAX);  // Two nodes found to add the Rz
            qubit_id = tree_expr.first.root;
          }
        }

        TKET_ASSERT(qubit_id != UINT_MAX);  // No node found to add the Rz

        std::vector<unsigned> qubit{qubit_id};

        synth_circuit.add_op(OpType::Rz, tree_expr.second, qubit);
        --tree_count;
        /* Otherwise, the tree stays in. */
      } else {
        CostedTrees::iterator iter = new_trees.find(tree_expr.first.tree_cost);
        if (iter == new_trees.end())
          new_trees.insert({tree_expr.first.tree_cost, {tree_expr}});
        else
          iter->second.push_back(tree_expr);
        unsigned cost_change = tree_expr.first.last_operation_cost;
        global_cost += cost_change;
      }
    }
  }
  current_trees = new_trees;
}

void SteinerForest::add_row_globally(unsigned i, unsigned j) {
  /* CNOT with control j and target i. Which way round the indices are is a
   * wee bit fiddly. */
  std::vector<unsigned> qbs = {j, i};
  synth_circuit.add_op(OpType::CX, qbs);
  linear_function.col_add(
      i, j);  // prepend a CNOT to the linear reversible function

  CostedTrees new_trees;

  for (const auto &cost_trees : current_trees) {
    for (auto tree_expr : cost_trees.second) {
      tree_expr.first.add_row(i, j);
      /* If the tree has been fully reduced, remove it from the
      forest */
      if (tree_expr.first.fully_reduced()) {
        std::vector<unsigned> qubit{i};
        synth_circuit.add_op(OpType::Rz, tree_expr.second, qubit);
        --tree_count;
      }
      /* Otherwise, update the forest */
      else {
        CostedTrees::iterator iter = new_trees.find(tree_expr.first.tree_cost);
        if (iter == new_trees.end())
          new_trees.insert({tree_expr.first.tree_cost, {tree_expr}});
        else
          iter->second.push_back(tree_expr);
        unsigned cost_change = tree_expr.first.last_operation_cost;
        global_cost += cost_change;
      }
    }
  }
  current_trees = new_trees;
}

void SteinerForest::add_operation_list(const OperationList &oper_list) {
  for (const Operation &operation : oper_list) {
    add_row_globally(operation.first, operation.second);
  }
}

OperationList SteinerForest::operations_available_under_the_index(
    const PathHandler &path, unsigned index) const {
  OperationList operations;
  for (unsigned i = 0; i != index; ++i) {
    CostedTrees::const_iterator iter = current_trees.find(i);
    if (iter == current_trees.end()) continue;
    for (const auto &tree_pair : iter->second) {
      operations.splice(
          operations.begin(), tree_pair.first.operations_available(path));
    }
  }
  return operations;
}

OperationList SteinerForest::operations_available_at_index(
    const PathHandler &path, unsigned index) const {
  OperationList operations;
  CostedTrees::const_iterator iter = current_trees.find(index);
  if (iter != current_trees.end()) {
    for (const auto &tree_pair : iter->second) {
      operations.splice(
          operations.begin(), tree_pair.first.operations_available(path));
    }
  }
  return operations;
}

OperationList SteinerForest::operations_available_at_min_costs(
    const PathHandler &path) const {
  OperationList operations;

  // search in all trees with the minimal costs for operations
  for (const auto &tree_expr : current_trees.begin()->second) {
    OperationList op_to_add = tree_expr.first.operations_available(path);

    operations.insert(operations.end(), op_to_add.begin(), op_to_add.end());
  }

  return operations;
}

CostedOperations best_operations_lookahead(
    const PathHandler &path, const SteinerForest &forest, unsigned lookahead) {
  if (lookahead == 0) {
    throw std::logic_error("Must look ahead at least one step");
  }
  CostedTrees::const_reverse_iterator r_iter = forest.current_trees.rbegin();
  if (r_iter == forest.current_trees.rend()) {
    throw std::logic_error("Forest is empty");
  }

  OperationList operations_available =
      forest.operations_available_at_min_costs(path);

  TKET_ASSERT(!operations_available.empty());  // Cannot find any operations

  // set up base case
  OperationList ops{operations_available.front()};
  CostedOperations costed_operations =
      recursive_operation_search(path, forest, lookahead - 1, ops);
  operations_available.pop_front();

  for (OperationList::iterator op_iter = operations_available.begin();
       op_iter != operations_available.end(); ++op_iter) {
    ops = {*op_iter};
    CostedOperations candidate_operations =
        recursive_operation_search(path, forest, lookahead - 1, ops);

    if ((candidate_operations.first < costed_operations.first) ||
        ((candidate_operations.first == costed_operations.first) &&
         (candidate_operations.second.size() <
          costed_operations.second.size()))) {
      costed_operations = std::move(candidate_operations);
    }
  }

  return costed_operations;
}

CostedOperations recursive_operation_search(
    const PathHandler &path, SteinerForest forest, unsigned lookahead,
    OperationList row_operations) {
  CostedOperations costed_operations;
  CostedOperations candidate_operations;

  forest.add_row_globally(
      row_operations.back().first, row_operations.back().second);

  if ((lookahead == 0) || (forest.current_trees.empty())) {
    return {forest.global_cost, row_operations};
  }
  CostedTrees::const_reverse_iterator r_iter = forest.current_trees.rbegin();
  unsigned index = r_iter->first;
  OperationList operations_available =
      forest.operations_available_under_the_index(path, index);
  if (operations_available.empty()) {
    return {forest.global_cost, row_operations};
  } else {
    row_operations.push_back(operations_available.front());
    costed_operations =
        recursive_operation_search(path, forest, lookahead - 1, row_operations);
    row_operations.pop_back();
    operations_available.pop_front();

    for (OperationList::iterator op_iter = operations_available.begin();
         op_iter != operations_available.end(); ++op_iter) {
      row_operations.push_back(*op_iter);
      candidate_operations = recursive_operation_search(
          path, forest, lookahead - 1, row_operations);

      row_operations.pop_back();

      if ((candidate_operations.first < costed_operations.first) ||
          ((candidate_operations.first == costed_operations.first) &&
           (candidate_operations.second.size() <
            costed_operations.second.size()))) {
        costed_operations = std::move(candidate_operations);
      }
    }
    return costed_operations;
  }
}

Circuit phase_poly_synthesis_int(
    const Architecture &arch, const PhasePolyBox &phasepolybox,
    unsigned lookahead, CNotSynthType cnottype) {
  if (lookahead == 0)
    throw std::logic_error(
        "[AAS] the lookahead of the phase polynomial synthesis has to be "
        "greater than 0");

  CostedOperations best_operations;

  PathHandler path(arch);

  PathHandler acyclic_path = path.construct_acyclic_handler();

  SteinerForest forest(acyclic_path, phasepolybox);

  while (!forest.current_trees.empty()) {
    best_operations =
        best_operations_lookahead(acyclic_path, forest, lookahead);
    forest.add_operation_list(best_operations.second);
  }
  Circuit cnot_circ(path.get_size());
  switch (cnottype) {
    case CNotSynthType::HamPath: {
      cnot_circ =
          aas_CNOT_synth(forest.linear_function, path, CNotSynthType::HamPath);

      // check if matrix is id, if not the result will be wrong
      TKET_ASSERT(forest.linear_function.is_id());
      break;
    }
    case CNotSynthType::Rec: {
      Circuit resultofstep =
          aas_CNOT_synth(forest.linear_function, path, CNotSynthType::Rec);
      cnot_circ = cnot_circ >> resultofstep;

      // check if matrix is id, if not the result will be wrong
      TKET_ASSERT(forest.linear_function.is_id());
      break;
    }
    case CNotSynthType::SWAP: {
      cnot_circ = aas_CNOT_synth_SWAP(forest.linear_function, path);
      // the check if the forest.linear_function is id is executed in the
      // aas_CNOT_synth_SWAP function
      break;
    }
    default: {
      TKET_ASSERT(!"[AAS]: unknown type of cnot synth");
    }
  }

  return forest.synth_circuit >> cnot_circ.dagger();
}

static PhasePolyBox make_placed_ppb(
    const Architecture &arch, const PhasePolyBox &phasepolybox) {
  Circuit circuit_ppb_place =
      phasepolybox.generate_circuit_with_original_placement();

  const std::string register_name = "surplus";

  unsigned qb_counter = circuit_ppb_place.n_qubits();
  while (arch.n_nodes() > circuit_ppb_place.n_qubits()) {
    Qubit qb = Qubit(register_name, qb_counter);
    circuit_ppb_place.add_qubit(qb);
    ++qb_counter;
  }

  TKET_ASSERT(circuit_ppb_place.n_qubits() == arch.n_nodes());

  qubit_vector_t q_vec_place = circuit_ppb_place.all_qubits();
  std::map<Qubit, Node> qubit_to_nodes_place;
  unsigned counter_place = 0;

  for (Node no_place : arch.nodes()) {
    if (counter_place < circuit_ppb_place.n_qubits()) {
      qubit_to_nodes_place.insert({q_vec_place[counter_place], no_place});
      ++counter_place;
    }
  }

  circuit_ppb_place.rename_units(qubit_to_nodes_place);

  return PhasePolyBox(circuit_ppb_place);
}

class PhasePolySynthesizer {
 public:
  PhasePolySynthesizer(
      const Architecture &arch, const PhasePolyBox &phasepolybox,
      unsigned lookahead, CNotSynthType cnottype)
      : arch(arch),
        placed_ppb(make_placed_ppb(arch, phasepolybox)),
        lookahead(lookahead),
        cnottype(cnottype) {}

  Circuit get_result() {
    return (cnottype == CNotSynthType::HamPath) ? get_result_using_hampath()
                                                : get_result_standard();
  }

 private:
  Circuit get_result_from_hampath(const std::vector<Node> &hampath) {
    // maps from qubits/node to int
    std::map<UnitID, UnitID> forward_contiguous_uids_q;
    std::map<UnitID, UnitID> backward_contiguous_uids_n;
    // extra map with node type needed for the creation of the architecture
    std::map<UnitID, Node> unitid_to_int_nodes;
    for (unsigned i = 0; i < hampath.size(); i++) {
      UnitID qu_no = UnitID(hampath[i]);
      Qubit q = Qubit(i);
      Node n = Node(i);
      unitid_to_int_nodes.insert({qu_no, n});
      forward_contiguous_uids_q.insert({q, qu_no});
      backward_contiguous_uids_n.insert({qu_no, n});
    }
    // define new architecture
    std::vector<Architecture::Connection> new_con;
    for (auto pair : arch.get_all_edges_vec()) {
      new_con.push_back(
          {unitid_to_int_nodes[UnitID(pair.first)],
           unitid_to_int_nodes[UnitID(pair.second)]});
    }
    Architecture con_arch(new_con);
    return make_result_from_con_arch(
        con_arch, forward_contiguous_uids_q, backward_contiguous_uids_n);
  }

  Circuit get_result_using_hampath() {
    std::vector<Node> hampath = find_hampath(arch);  // using default timeout
    Circuit c0 = get_result_from_hampath(hampath);
    // Sometimes the reversed path gives a better circuit. Try both!
    std::reverse(hampath.begin(), hampath.end());
    Circuit c1 = get_result_from_hampath(hampath);
    return (c0.depth() < c1.depth()) ? c0 : c1;
  }

  Circuit get_result_standard() {
    // calculate iteration order
    IterationOrder iter_order(arch);
    std::vector<Node> node_order = iter_order.get_iterationorder();
    // maps from qubits/node to int
    std::map<UnitID, UnitID> forward_contiguous_uids_q;
    std::map<UnitID, UnitID> backward_contiguous_uids_n;
    // extra map with node type needed for the creation of the architecture
    std::map<UnitID, Node> unitid_to_int_nodes;
    for (unsigned i = 0; i < node_order.size(); i++) {
      // convert node to superclass type of qubit/node
      UnitID qu_no = UnitID(node_order[i]);
      Qubit q = Qubit(i);
      Node n = Node(i);
      forward_contiguous_uids_q.insert({q, qu_no});
      backward_contiguous_uids_n.insert({qu_no, n});
      unitid_to_int_nodes.insert({qu_no, n});
    }
    // define new architecture: include only the tree edges
    std::vector<Architecture::Connection> new_con;
    for (auto pair : iter_order.get_edgelist()) {
      new_con.push_back(
          {unitid_to_int_nodes[UnitID(pair.first)],
           unitid_to_int_nodes[UnitID(pair.second)]});
    }
    Architecture con_arch(new_con);
    return make_result_from_con_arch(
        con_arch, forward_contiguous_uids_q, backward_contiguous_uids_n);
  }

  Circuit make_result_from_con_arch(
      const Architecture &con_arch,
      const std::map<UnitID, UnitID> &forward_contiguous_uids_q,
      const std::map<UnitID, UnitID> &backward_contiguous_uids_n) {
    // define new phase poly box
    Circuit circuit_ppb = placed_ppb.generate_circuit_with_original_placement();
    // the aas code is implemented under the assumption that all qubits in the
    // circuit are named from 0 to n. The same assumption was made for the nodes
    // of the architecture. To make sure that this condition is fulfilled the
    // qubits in the circuit and the architecture are renamed. The new names are
    // reverted at the end of the aas procedure. The qubits and the nodes have
    // the same name in the input.
    circuit_ppb.rename_units(backward_contiguous_uids_n);
    PhasePolyBox new_ppb(circuit_ppb);
    Circuit result =
        phase_poly_synthesis_int(con_arch, new_ppb, lookahead, cnottype);
    result.rename_units(forward_contiguous_uids_q);
    return result;
  }

  Architecture arch;
  PhasePolyBox placed_ppb;
  unsigned lookahead;
  CNotSynthType cnottype;
};

Circuit phase_poly_synthesis(
    const Architecture &arch, const PhasePolyBox &phasepolybox,
    unsigned lookahead, CNotSynthType cnottype) {
  PhasePolySynthesizer pps(arch, phasepolybox, lookahead, cnottype);
  return pps.get_result();
}

}  // namespace aas
}  // namespace tket
