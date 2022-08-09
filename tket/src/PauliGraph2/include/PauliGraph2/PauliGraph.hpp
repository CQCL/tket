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

#include "Clifford/SymplecticTableau.hpp"
#include "Utils/Expression.hpp"
#include "Utils/GraphHeaders.hpp"
#include "Utils/PauliStrings.hpp"
#include "Utils/SequencedContainers.hpp"

namespace tket {

// Forwards declarations for friends and internal uses
namespace pg {
class PauliGraph;
class PGOp;
typedef std::shared_ptr<const PGOp> PGOp_ptr;
}
class Circuit;
class Op;
typedef std::shared_ptr<const Op> Op_ptr;

pg::PauliGraph circuit_to_pauli_graph2(const Circuit &circ);
Circuit pauli_graph2_to_circuit_individual(
    const pg::PauliGraph &pg, CXConfigType cx_config);
Circuit pauli_graph2_to_circuit_pairwise(
    const pg::PauliGraph &pg, CXConfigType cx_config);
Circuit pauli_graph2_to_circuit_sets(
    const pg::PauliGraph &pg, CXConfigType cx_config);
std::vector<pg::PGOp_ptr> op_to_pgops(const Op_ptr& op, const unit_vector_t& args, pg::PauliGraph& pg, bool allow_tableau);

namespace pg {

class PGError : public std::logic_error {
 public:
  explicit PGError(const std::string& message) : std::logic_error(message) {}
};

enum class PGOpType {
  // Conventional Pauli Gadget, a rotation formed by exponentiating a Pauli tensor
  Rotation,

  // Clifford-angled Pauli Gadget
  CliffordRot,

  // A measurement in a multi-qubit Pauli basis
  Measure,

  // Decoherence in a multi-qubit Pauli basis (measurement ignoring the outcome)
  Decoherence,

  // Reset of a qubit, conjugated by a Clifford circuit
  Reset,

  // Some other PGOp conditioned on classical data
  Conditional,

  // An opaque boxed circuit component; treated as a local barrier
  // Defined in Converters module to have access to Circuit components
  Box
};

/**
 * Abstract class for a Pauli Graph Op.
 * Each PGOpType has a single possible subclass that can realise it, allowing 
 * us to statically cast to a subclass once that is determined.
 * 
 * Currently, each subclass of PGOp has a unique interpretation, with each
 * associated to a PGOpType for easy dynamic inspection.
 */
class PGOp {
 public:
  PGOpType get_type() const;

  virtual SymSet free_symbols() const = 0;

  virtual PGOp_ptr symbol_substitution(const SymEngine::map_basic_basic& sub_map) const = 0;

  virtual std::string get_name(bool latex = false) const = 0;

  bool operator==(const PGOp& other) const;

  /**
   * Checks equality between two instances of the same class.
   * The PGOp object passed as parameter must always be of the same type as this.
   *
   * For base class PGOp, it is sufficient that they have same type
   */
  virtual bool is_equal(const PGOp& other) const;

  bool commutes_with(const PGOp& other) const;

  virtual std::vector<QubitPauliTensor> active_paulis() const = 0;

  virtual bit_vector_t read_bits() const;

  virtual bit_vector_t write_bits() const;

  virtual ~PGOp();

 protected:
  PGOp(PGOpType type);
  const PGOpType type_;
};

class PGRotation : public PGOp {
 public:
  const QubitPauliTensor& get_tensor() const;

  const Expr& get_angle() const;

  PGRotation(const QubitPauliTensor& tensor, const Expr& angle);

  // Overrides from PGOp
  virtual SymSet free_symbols() const override;
  virtual PGOp_ptr symbol_substitution(const SymEngine::map_basic_basic& sub_map) const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual std::vector<QubitPauliTensor> active_paulis() const override;

 protected:
  QubitPauliTensor tensor_;
  Expr angle_;
};

class PGCliffordRot : public PGOp {
 public:
  const QubitPauliTensor& get_tensor() const;
  unsigned get_angle() const;

  PGCliffordRot(const QubitPauliTensor& tensor, unsigned angle);

  // Overrides from PGOp
  virtual SymSet free_symbols() const override;
  virtual PGOp_ptr symbol_substitution(const SymEngine::map_basic_basic& sub_map) const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual std::vector<QubitPauliTensor> active_paulis() const override;

 protected:
  QubitPauliTensor tensor_;
  unsigned angle_;
};

class PGMeasure : public PGOp {
 public:
  const QubitPauliTensor& get_tensor() const;
  const Bit& get_target() const;

  PGMeasure(const QubitPauliTensor& tensor, const Bit& target);

  // Overrides from PGOp
  virtual SymSet free_symbols() const override;
  virtual PGOp_ptr symbol_substitution(const SymEngine::map_basic_basic& sub_map) const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual std::vector<QubitPauliTensor> active_paulis() const override;
  virtual bit_vector_t write_bits() const override;

 protected:
  QubitPauliTensor tensor_;
  Bit target_;
};

class PGDecoherence : public PGOp {
 public:
  const QubitPauliTensor& get_tensor() const;

  PGDecoherence(const QubitPauliTensor& tensor);

  // Overrides from PGOp
  virtual SymSet free_symbols() const override;
  virtual PGOp_ptr symbol_substitution(const SymEngine::map_basic_basic& sub_map) const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual std::vector<QubitPauliTensor> active_paulis() const override;

 protected:
  QubitPauliTensor tensor_;
};

class PGReset : public PGOp {
 public:
  const QubitPauliTensor& get_stab() const;
  const QubitPauliTensor& get_destab() const;

  PGReset(const QubitPauliTensor& stab, const QubitPauliTensor& destab);

  // Overrides from PGOp
  virtual SymSet free_symbols() const override;
  virtual PGOp_ptr symbol_substitution(const SymEngine::map_basic_basic& sub_map) const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual std::vector<QubitPauliTensor> active_paulis() const override;

 protected:
  QubitPauliTensor stab_;
  QubitPauliTensor destab_;
};

class PGConditional : public PGOp {
 public:
  PGOp_ptr get_inner_op() const;
  const bit_vector_t& get_args() const;
  unsigned get_value() const;

  PGConditional(PGOp_ptr inner, const bit_vector_t& args, unsigned value);

  // Overrides from PGOp
  virtual SymSet free_symbols() const override;
  virtual PGOp_ptr symbol_substitution(const SymEngine::map_basic_basic& sub_map) const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual std::vector<QubitPauliTensor> active_paulis() const override;
  virtual bit_vector_t read_bits() const override;
  virtual bit_vector_t write_bits() const override;

 protected:
  PGOp_ptr inner_;
  bit_vector_t args_;
  unsigned value_;
};

struct PGVertProperties {
  PGOp_ptr op_;
};

struct PGEdgeProperties {};

typedef boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS, PGVertProperties, PGEdgeProperties> PGDAG;
typedef boost::graph_traits<PGDAG>::vertex_descriptor PGVert;
typedef boost::graph_traits<PGDAG>::edge_descriptor PGEdge;

typedef sequence_set_t<PGVert> PGVertSet;
typedef sequence_set_t<PGEdge> PGEdgeSet;

/**
 * Distinguish between Z and X rows of a tableau for a given qubit
 */
enum TableauRowType {
  ZRow,
  XRow,
};
typedef boost::bimap<std::pair<Qubit, TableauRowType>, unsigned> tableau_row_index_t;
typedef boost::bimap<Qubit, unsigned> tableau_col_index_t;

/**
 * Dependency graph of a circuit designed to abstract away Cliffords by focussing on Pauli gadgets.
 * Consists of a dependency graph of operations with Clifford tableaux at both the start and end for compatibility with either normal form.
 * The initial tableau is an isometry tableau to handle both input and zero qubits.
 * The final tableau is a partial tableau, only specifying the Pauli strings mapped to the live outputs, ignoring those that are discarded.
 */
class PauliGraph {
 public:
  /** Construct an empty dependency graph for the identity over n qubits. */
  explicit PauliGraph(unsigned n);

  /** Construct an empty dependency graph for the identity over given qubits. */
  explicit PauliGraph(const qubit_vector_t &qbs, const bit_vector_t &bits = {});

  /** Visualisation of the dependency graph */
  void to_graphviz(std::ostream &out) const;

  unsigned n_vertices() const;

  friend PauliGraph tket::circuit_to_pauli_graph2(const tket::Circuit &circ);
  friend tket::Circuit tket::pauli_graph2_to_circuit_individual(
      const PauliGraph &pg, CXConfigType cx_config);
  friend tket::Circuit tket::pauli_graph2_to_circuit_pairwise(
      const PauliGraph &pg, CXConfigType cx_config);
  friend tket::Circuit tket::pauli_graph2_to_circuit_sets(
      const PauliGraph &pg, CXConfigType cx_config);
  friend std::vector<PGOp_ptr> tket::op_to_pgops(const tket::Op_ptr& op, const unit_vector_t& args, PauliGraph& pg, bool allow_tableau);

 private:
  PGDAG graph_;
  SymplecticTableau initial_;
  tableau_row_index_t initial_rows_;
  tableau_col_index_t initial_cols_;
  SymplecticTableau final_;
  tableau_row_index_t final_rows_;
  tableau_col_index_t final_cols_;

  bit_vector_t bits_;

  PGVertSet start_line_;
  PGVertSet end_line_;

  PGVertSet get_successors(const PGVert& vert) const;
  PGVertSet get_predecessors(const PGVert& vert) const;
  PGEdgeSet get_in_edges(const PGVert& vert) const;
  PGEdgeSet get_out_edges(const PGVert& vert) const;
  PGVert source(const PGEdge& edge) const;
  PGVert target(const PGEdge& edge) const;

  /**
   * Appends a vertex the end of the dependency graph.
   * Assumes this is the result AFTER pushing it through the final
   * tableau.
   */
  void add_vertex_at_end(PGOp_ptr op);

  /**
   * Iterates through the vertices of a PauliGraph in a topological ordering.
   * When there are multiple commuting vertices that could be emitted, the
   * order of emission is deterministic but arbitrary.
   */
  class TopSortIterator {
   public:
    TopSortIterator();
    explicit TopSortIterator(const PauliGraph &pg);

    const PGVert &operator*() const;
    const PGVert *operator->() const;
    bool operator==(const TopSortIterator &other) const;
    bool operator!=(const TopSortIterator &other) const;

    TopSortIterator operator++(int);
    TopSortIterator &operator++();

   private:
    const PauliGraph *pg_;
    PGVert current_vert_;
    PGVertSet search_set_;
    std::unordered_set<PGVert> visited_;
  };

  TopSortIterator begin() const;
  TopSortIterator end() const;

};

}  // namespace pg
}  // namespace tket
