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

#include "tket/Clifford/ChoiMixTableau.hpp"
#include "tket/Utils/GraphHeaders.hpp"
#include "tket/Utils/PauliStrings.hpp"

namespace tket {

// Forwards declarations for friends and internal uses
namespace pg {
class PauliGraph;
class PGOp;
// Not const as we wish for these to be updated in-place
typedef std::shared_ptr<PGOp> PGOp_ptr;
}  // namespace pg
class Circuit;
class Op;
typedef std::shared_ptr<const Op> Op_ptr;

pg::PauliGraph circuit_to_pauli_graph3(const Circuit& circ);
Circuit pauli_graph3_to_circuit_individual(
    const pg::PauliGraph& pg, CXConfigType cx_config);

namespace pg {

class PGError : public std::logic_error {
 public:
  explicit PGError(const std::string& message) : std::logic_error(message) {}
};

enum class PGOpType {
  // Conventional Pauli Gadget, a rotation formed by exponentiating a Pauli
  // tensor
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
  Box,

  // A time-symmetric view of an imposed stabilizer/projection into a subspace
  // At synthesis, we can choose whether to impose this stabilizer by inclusion
  // as a row of the initial tableau or explicitly via a mid- or end-of-circuit
  // postselection
  Stabilizer,

  // The initial tableau
  // The active QubitPauliTensors are from the output segment of the tableau,
  // i.e. the segment that connects to the interior of the Pauli Graph
  InputTableau,

  // The final tableau
  // The active QubitPauliTensors are from the input segment of the tableau,
  // i.e. the segment that connects to the interior of the Pauli Graph
  OutputTableau,
};

/**
 * Abstract class for a Pauli Graph Op.
 * Each PGOpType has a single possible subclass that can realise it, allowing
 * us to statically cast to a subclass once that is determined.
 *
 * Currently, each subclass of PGOp has a unique interpretation, with each
 * associated to a PGOpType for easy dynamic inspection.
 *
 * This falls in line more so with Command than Op as each instance of a PGOp
 * relates to a specific cluster of Paulis within a given PauliGraph.
 */
class PGOp {
 public:
  PGOpType get_type() const;

  virtual SymSet free_symbols() const = 0;

  virtual PGOp_ptr symbol_substitution(
      const SymEngine::map_basic_basic& sub_map) const = 0;

  virtual std::string get_name(bool latex = false) const = 0;

  bool operator==(const PGOp& other) const;

  /**
   * Checks equality between two instances of the same class.
   * The PGOp object passed as parameter must always be of the same type as
   * this.
   *
   * For base class PGOp, it is sufficient that they have same type
   */
  virtual bool is_equal(const PGOp& other) const = 0;

  bool commutes_with(const PGOp& other) const;

  virtual unsigned n_paulis() const;

  virtual std::vector<QubitPauliTensor> active_paulis() const = 0;

  virtual QubitPauliTensor& port(unsigned p) = 0;

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
  virtual PGOp_ptr symbol_substitution(
      const SymEngine::map_basic_basic& sub_map) const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual std::vector<QubitPauliTensor> active_paulis() const override;
  virtual QubitPauliTensor& port(unsigned p) override;

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
  virtual PGOp_ptr symbol_substitution(
      const SymEngine::map_basic_basic& sub_map) const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual std::vector<QubitPauliTensor> active_paulis() const override;
  virtual QubitPauliTensor& port(unsigned p) override;

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
  virtual PGOp_ptr symbol_substitution(
      const SymEngine::map_basic_basic& sub_map) const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual std::vector<QubitPauliTensor> active_paulis() const override;
  virtual QubitPauliTensor& port(unsigned p) override;
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
  virtual PGOp_ptr symbol_substitution(
      const SymEngine::map_basic_basic& sub_map) const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual std::vector<QubitPauliTensor> active_paulis() const override;
  virtual QubitPauliTensor& port(unsigned p) override;

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
  virtual PGOp_ptr symbol_substitution(
      const SymEngine::map_basic_basic& sub_map) const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual unsigned n_paulis() const override;
  virtual std::vector<QubitPauliTensor> active_paulis() const override;
  virtual QubitPauliTensor& port(unsigned p) override;

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
  virtual PGOp_ptr symbol_substitution(
      const SymEngine::map_basic_basic& sub_map) const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual unsigned n_paulis() const override;
  virtual std::vector<QubitPauliTensor> active_paulis() const override;
  virtual QubitPauliTensor& port(unsigned p) override;
  virtual bit_vector_t read_bits() const override;
  virtual bit_vector_t write_bits() const override;

 protected:
  PGOp_ptr inner_;
  bit_vector_t args_;
  unsigned value_;
};

class PGStabilizer : public PGOp {
 public:
  const QubitPauliTensor& get_stab() const;

  PGStabilizer(const QubitPauliTensor& stab);

  // Overrides from PGOp
  virtual SymSet free_symbols() const override;
  virtual PGOp_ptr symbol_substitution(
      const SymEngine::map_basic_basic& sub_map) const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual std::vector<QubitPauliTensor> active_paulis() const override;
  virtual QubitPauliTensor& port(unsigned p) override;

 protected:
  QubitPauliTensor stab_;
};

class PGInputTableau : public PGOp {
 public:
  // Returns the row tensor as from the tableau; first component is for the
  // input segment, second for the output component (the active paulis)
  const ChoiMixTableau::row_tensor_t& get_full_row(unsigned p) const;

  PGInputTableau(const ChoiMixTableau& tableau);

  // Overrides from PGOp
  virtual SymSet free_symbols() const override;
  virtual PGOp_ptr symbol_substitution(
      const SymEngine::map_basic_basic& sub_map) const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual unsigned n_paulis() const override;
  virtual std::vector<QubitPauliTensor> active_paulis() const override;
  virtual QubitPauliTensor& port(unsigned p) override;

 protected:
  // Store the rows as QubitPauliTensors rather than an actual tableau object
  // for easier modification of individual rows in the same way as for rewriting
  // on other PGOps
  std::vector<ChoiMixTableau::row_tensor_t> rows_;
};

class PGOutputTableau : public PGOp {
 public:
  // Returns the row tensor as from the tableau; first component is for the
  // input segment (the active paulis), second for the output component
  const ChoiMixTableau::row_tensor_t& get_full_row(unsigned p) const;

  PGOutputTableau(const ChoiMixTableau& tableau);

  // Overrides from PGOp
  virtual SymSet free_symbols() const override;
  virtual PGOp_ptr symbol_substitution(
      const SymEngine::map_basic_basic& sub_map) const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual unsigned n_paulis() const override;
  virtual std::vector<QubitPauliTensor> active_paulis() const override;
  virtual QubitPauliTensor& port(unsigned p) override;

 protected:
  // Store the rows as QubitPauliTensors rather than an actual tableau object
  // for easier modification of individual rows in the same way as for rewriting
  // on other PGOps
  std::vector<ChoiMixTableau::row_tensor_t> rows_;
};

/**
 * PauliGraph
 *
 * This data structure provides a balance between the simple rewriting of an
 * instruction graph (with arcs between operations sharing the same physical
 * resource, e.g. Circuit) and the abstraction of a dependency DAG (abstracts
 * away all commutations).
 *
 * We attribute each instruction to a small number of Pauli strings, with the
 * guarantee that if each string from A commutes with each string of B then A
 * and B commute (this is a safe under-approximation of commutativity - there
 * may be commutations this doesn't identify). Rewriting requires us to update
 * the Pauli strings and the relation of anticommutations between the strings.
 *
 * We separately use a true dependency DAG for the classical dependencies (i.e.
 * there is a single edge between two operations if reordering them would cause
 * a RAW, WAR, or WAW hazard).
 *
 * We intend to support the following rewrites during optimisation:
 * - Reordering commuting operations
 * - Pauli reorder rules (just updating phases of strings)
 * - Clifford reorder rules (updating Pauli strings by multiplication)
 * - Merging compatible vertices (rotations, measurements, discards, etc.)
 * - "Product Rotation Lemma" actions (multiplies a Pauli string by a
 * stabilizer; see Simmons 2021)
 * - Deletion of identity vertices
 * - Deletions of vertices at start and end
 * - Absorbing Cliffords into the start and end tableaux
 * - Changing vertex types (e.g. continuously-parameterised rotation to discrete
 * Clifford rotation, reset expansion)
 *
 * Each operation corresponds to exactly one node in the classical graph but may
 * use multiple Pauli strings, so we attach operation details to the vertices of
 * the classical graph. The heterogeneity of contents for different kinds of
 * operations encourages an object-oriented structure for node contents, similar
 * to Ops in Circuits. Unlike Ops, the large variability in Pauli strings means
 * we won't benefit significantly from reusing immutable objects, so we instead
 * store separate objects for each vertex and allow them to be mutable to update
 * in-place where possible.
 *
 * Few rewrites will update the classical data so maintaining the classical
 * dependency for fast lookup is best (as opposed to maintaining a candidate
 * temporal ordering of the operations and determining classical dependencies on
 * the fly). Dependencies are typically sparse, so a directed adjacency list is
 * suitable.
 *
 * Some additional lookup maps maintain the most recent reads and writes to each
 * classical Bit to aid vertex insertion. These will be largely unimportant when
 * it comes to rewriting though.
 *
 * We store the anticommutation between the Pauli strings of different
 * operations to save recalculating them a lot on the fly. We specifically store
 * a directed form of the anticommutation relation that also factors in the
 * ordering of the operations, i.e. (P, Q) means both P and Q anticommute and
 * P's operation occurs after Q's. This can be a relatively dense relation and
 * updates due to multiplying strings involve taking XOR or symmetric difference
 * between the ancestors/descendants, so we store it as a Binary matrix for easy
 * updating via row/column updates. Row i indicates the anticommuting ancestors
 * (earlier in the circuit) of Pauli i, and column i indicates the anticommuting
 * descendants (later in the circuit).
 *
 * During rewrites, once we have decided on a vertex to rewrite around, we will
 * need to both find the rows/columns in the anticommutation matrix
 * corresponding to a particular vertex. Often the entries in the matrix will
 * then inform which other vertices need to be rewritten, e.g. when moving a
 * Clifford instruction to the start of the circuit, the positive indices in its
 * row give the ancestors that need to be updated, so we also need a reverse
 * lookup from the table indices. It is easiest to maintain this mapping as a
 * multi-indexed container, allowing other data to also be attached if needed in
 * the future.
 *
 * Each Pauli string within the PauliGraph can be uniquely identified either by
 * its index in the anticommutation matrix, or by a combination of the vertex
 * and index of the PauliString within the PGOp, referred to as its port. The
 * number of ports and their ordering/interpretation is fixed based on the
 * PGOpType/subclass of PGOp.
 *
 * During rewrites which eliminate vertices, we leave unused rows/columns in the
 * anticommutation matrix rather than attempt to reduce it at every opportunity.
 * A cleanup method can be written if we wish to run this occasionally during
 * long rewrite procedures.
 *
 * Whilst previous iterations of PauliGraph contained an explicit Clifford
 * tableau at the start or end of the circuit, we choose to represent these
 * within the graph itself, since including them in the anticommutation matrix
 * allows for easy identification of opportunities for eliminating instructions
 * around discards or stabilizers, or applying PRL actions. In the case where we
 * need to relate Pauli strings to inputs or outputs, we follow the style of
 * ChoiMixedTableau in describing pairs of related Pauli strings over the inputs
 * and interior or over the interior and outputs. However, we only care about
 * the interior Pauli strings in the anticommutation matrix.
 *
 * When a vertex may contain multiple ports, such as InputTableau and
 * OutputTableau, we view the actions on the ports as happening simultaneously,
 * so the anticommutation matrix will read false in the corresponding entries
 * even if the Pauli strings anticommute.
 */

typedef boost::adjacency_list<
    boost::setS, boost::listS, boost::bidirectionalS, PGOp_ptr>
    PGClassicalGraph;
typedef PGClassicalGraph::vertex_descriptor PGVert;

struct PGPauli {
  unsigned index;
  PGVert vert;
  unsigned port;
};

struct TagID {};
struct TagOp {};

typedef boost::multi_index::multi_index_container<
    PGPauli,
    boost::multi_index::indexed_by<
        boost::multi_index::ordered_unique<
            boost::multi_index::tag<TagID>,
            boost::multi_index::member<PGPauli, unsigned, &PGPauli::index>>,
        boost::multi_index::hashed_non_unique<
            boost::multi_index::tag<TagOp>,
            boost::multi_index::member<PGPauli, PGVert, &PGPauli::vert>>>>
    PGIndex;

class PauliGraph {
 public:
  explicit PauliGraph();
  explicit PauliGraph(
      const std::set<Qubit>& qubits, const std::set<Bit>& bits = {});
  void to_graphviz(std::ostream& out) const;

  PGVert add_vertex_at_end(PGOp_ptr op);

  // Verification of validity of the data structure
  // Expensive so intended for use in debugging and tests, but not live code
  void verify() const;

  friend PauliGraph tket::circuit_to_pauli_graph3(const tket::Circuit& circ);
  friend tket::Circuit tket::pauli_graph3_to_circuit_individual(
      const PauliGraph& pg, CXConfigType cx_config);
  friend std::vector<PGOp_ptr> tket::op_to_pgops(
      const tket::Op_ptr& op, const unit_vector_t& args, PauliGraph& pg,
      bool allow_tableau);

 private:
  MatrixXb pauli_ac_;
  PGIndex pauli_index_;
  PGClassicalGraph c_graph_;
  std::set<Qubit> qubits_;
  std::set<Bit> bits_;

  std::map<Bit, PGVert> last_writes_;
  std::map<Bit, std::unordered_set<PGVert>> last_reads_;

  std::optional<PGVert> initial_tableau_;
  std::optional<PGVert> final_tableau_;

  // Replaces the QubitPauliString of row target_r with coeff * source * target
  // and updates pauli_ac_ accordingly
  void multiply_strings(
      unsigned source_r, unsigned target_r, Complex coeff = 1.);
};

}  // namespace pg
}  // namespace tket
