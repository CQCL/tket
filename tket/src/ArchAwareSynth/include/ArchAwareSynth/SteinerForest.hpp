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

#pragma once
#include "Circuit/Circuit.hpp"
#include "Converters/PhasePoly.hpp"
#include "Graphs/DirectedGraph.hpp"
#include "Graphs/Utils.hpp"
#include "SteinerTree.hpp"
#include "Utils/GraphHeaders.hpp"
#include "Utils/UnitID.hpp"

namespace tket {
namespace aas {

typedef std::map<unsigned, std::list<std::pair<SteinerTree, Expr>>> CostedTrees;
typedef std::vector<CostedTrees> TrialCostedTrees;
typedef std::list<unsigned> ParityList;
typedef std::pair<unsigned, OperationList> CostedOperations;

class NoHamiltonPath : public std::logic_error {
 public:
  explicit NoHamiltonPath(const std::string &message)
      : std::logic_error(message) {}
};

/**
 * Represents a series of sequential operations on an architecture,
 * each of which are represented by Steiner trees. The prototypical
 * example is a phase polynomial.
 */
class SteinerForest {
 public:
  /**
   * construct a steinerforst from a architecture and a phasepolybox
   * @param arch architecture for construction
   * @param phasepolybox box for construction
   */
  explicit SteinerForest(
      const Architecture &arch, const PhasePolyBox &phasepolybox)
      : SteinerForest(PathHandler(arch), phasepolybox) {}

  /**
   * construct a steinerforst from a PathHandler and a phasepolybox
   * @param paths PathHandler for construction
   * @param phasepolybox box for construction
   */
  explicit SteinerForest(
      const PathHandler &paths, const PhasePolyBox &phasepolybox);

  CostedTrees current_trees;
  TrialCostedTrees test_trees;
  Circuit synth_circuit;
  DiagMatrix linear_function;

  unsigned global_cost;
  unsigned tree_count;

  /**
   * add operation to all trees of the forest
   * @param i control index for operation
   * @param j target index for operation
   */
  void add_row_globally(unsigned i, unsigned j);

  /**
   * add a list of operation to the forest
   * @param oper_list list of operation to be added
   */
  void add_operation_list(const OperationList &oper_list);

  /**
   * finds an exhaustive list of operations which may be performed for trees
   * under a specified cost index
   * @param path pathhandler used for the search
   * @param index maximum cost of the trees included in the calculation
   */
  OperationList operations_available_under_the_index(
      const PathHandler &path, unsigned index) const;

  /**
   * finds an exhaustive list of operations which may be performed for trees
   * at a specified cost index. If there are more trees with the same costs, all
   * are included in the search
   * @param path pathhandler used for the search
   * @param index cost of the trees included in the calculation
   */
  OperationList operations_available_at_index(
      const PathHandler &path, unsigned index) const;

  /**
   * finds an exhaustive list of operations which may be performed for trees
   * with the minimal cost index. If there are more trees with the same costs,
   * all are included in the search
   * @param path pathhandler used for the search
   */
  OperationList operations_available_at_min_costs(
      const PathHandler &path) const;
};

/**
 * searches for the best operation in the given forest
 * @param path pathhandler used for the calculation
 * @param forest steinerforest used for the calculation
 * @param lookahead maximum steps of recursion used for the iteration
 */
CostedOperations best_operations_lookahead(
    const PathHandler &path, const SteinerForest &forest, unsigned lookahead);

/**
 * searches for the best operation in the given forest with operation which are
 * applied to the forest before the search is started
 * @param path pathhandler used for the calculation
 * @param forest steinerforest used for the calculation
 * @param lookahead maximum steps of recursion used for the iteration
 * @param row_operations operations which are executed before the search starts
 */
CostedOperations recursive_operation_search(
    const PathHandler &path, SteinerForest forest, unsigned lookahead,
    OperationList row_operations);

/**
 * function for architecture aware synthesis without the rename of the qubits
 * and nodes of the architecture. this function is asserting, that all qubits
 * are placed to nodes and they both are named with increasing integers. This
 * function will return a routed version of a given phase poly box. The
 * algorithm used for the cnot synthesis can be given by a parameter. The option
 * are a recursive algorithm, a swap based algorithm and a iterative algorithm.
 * The iterative algorithm needs a Hamilton path in the architecture. If there
 * is no Hamilton path in the architecture the function will throw a
 * logic_error.
 * @param arch architecture including all allowed edges for the routing
 * @param phasepolybox giving the phase poly box
 * @param lookahead giving the maximum iteration depth for the cnot+rz synthesis
 * @param cnottype type of cnot synthesis, allowing CNotSynthType::Rec,
 * CNotSynthType::HamPath or CNotSynthType::SWAP
 * @return routed circuit
 */
Circuit phase_poly_synthesis_int(
    const Architecture &arch, const PhasePolyBox &phasepolybox,
    unsigned lookahead = 1, CNotSynthType cnottype = CNotSynthType::Rec);

/**
 * main function for architecture aware synthesis in tket at the moment.
 * This function will return a routed version of a given phase poly box. The
 * algorithm used for the cnot synthesis can be given by a parameter. The option
 * are a recursive algorithm, a swap based algorithm and a iterative algorithm.
 * The iterative algorithm needs a Hamilton path in the architecture. If there
 * is no Hamilton path in the architecture the function will throw a
 * logic_error.
 * @param arch architecture including all allowed edges for the routing
 * @param phasepolybox giving the phase poly box
 * @param lookahead giving the maximum iteration depth for the cnot+rz synthesis
 * @param cnottype type of cnot synthesis, allowing CNotSynthType::Rec,
 * CNotSynthType::HamPath or CNotSynthType::SWAP
 * @return routed circuit
 */
Circuit phase_poly_synthesis(
    const Architecture &arch, const PhasePolyBox &phasepolybox,
    unsigned lookahead, CNotSynthType cnottype = CNotSynthType::Rec);

}  // namespace aas
}  // namespace tket
