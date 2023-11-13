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
#include "tket/Clifford/UnitaryTableau.hpp"
#include "tket/Utils/GraphHeaders.hpp"
#include "tket/Utils/PauliTensor.hpp"

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

  // An embedding of a StabiliserAssertionBox
  // Describes an ancilla qubit state, a target measurement bit, and a Pauli
  // string across the rest.
  // The semantics is that the ancilla qubit is reset, then the Pauli string
  // measured along it and recorded in the target bit.
  StabAssertion,

  // The initial tableau
  // The active SpPauliStabilisers are from the output segment of the tableau,
  // i.e. the segment that connects to the interior of the Pauli Graph
  InputTableau,

  // The final tableau
  // The active SpPauliStabilisers are from the input segment of the tableau,
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
  /**
   * Returns the type of PGOp, allowing us to determine the subclass of an
   * instance at runtime.
   */
  PGOpType get_type() const;

  /**
   * Returns the set of symbols used in any symbolic parameters of the PGOp.
   */
  virtual SymSet free_symbols() const = 0;

  /**
   * Performs symbolic substitution in any symbolic parameters of the PGOp.
   *
   * If the PGOp subclass uses symbolic parameters, this returns the result of
   * the substitution as a new PGOp. Otherwise, this returns an empty pointer.
   */
  virtual PGOp_ptr symbol_substitution(
      const SymEngine::map_basic_basic& sub_map) const = 0;

  /**
   * A human-readable summary of the PGOp.
   */
  virtual std::string get_name(bool latex = false) const = 0;

  /**
   * Equality check between any two PGOps.
   *
   * First compares their PGOpType to determine whether or not they are
   * instances of the same subclass, and if so uses the relevant is_equal()
   * implementation.
   */
  bool operator==(const PGOp& other) const;

  /**
   * Checks equality between two instances of the same class.
   * The PGOp object passed as parameter must always be of the same type as
   * this.
   *
   * For base class PGOp, it is sufficient that they have same type
   */
  virtual bool is_equal(const PGOp& other) const = 0;

  /**
   * Performs an efficient and safely under-estimating check of commutation
   * (i.e. returning true means they definitely commute, but returning false
   * means it is unlikely they commute). Checks whether all active_paulis
   * mutually commute between the two PGOps.
   */
  bool commutes_with(const PGOp& other) const;

  /**
   * Returns the size of active_paulis, i.e. a measure of the size of the
   * subspace of the Pauli group on which this operator acts non-trivially.
   */
  virtual unsigned n_paulis() const;

  /**
   * Returns a collection of Pauli operators dictating the subspace on which the
   * op acts non-trivially. The guarantee is that, if another op commutes with
   * all Pauli operators in active_paulis, then it commutes with the PGOp (the
   * converse need not hold, for example Rotation gates with angle 0).
   *
   * The ordering of the Pauli operators may be set by the semantics of the
   * subclass, e.g. the projected stabiliser of a PGReset is the Pauli operator
   * at port 0 and the lost stabiliser is at port 1.
   *
   * SpPauliStabiliser is used to account for phase information in common
   * updates and rewrites (e.g. Clifford reordering rules). Some PGOpTypes won't
   * be phase-sensitive (e.g. Decoherence) and some may double-up on phase
   * information (e.g. CliffordRot(P,3) is the same as CliffordRot(-P,1)), but
   * having just +- phase info on the easily accessible PauliTensors is a
   * reasonable middle ground and the other cases can be easily handled on an
   * ad-hoc basis.
   */
  virtual std::vector<SpPauliStabiliser> active_paulis() const = 0;

  /**
   * Gives direct reference access to the SpPauliStabiliser at index \p p in
   * active_paulis. This is most useful to give immediate, generic access to the
   * active_paulis for rewrites and synthesis without having to inspect the
   * PGOpType and cast to the appropriate subclass.
   */
  virtual SpPauliStabiliser& port(unsigned p) = 0;

  /**
   * The classical bits this PGOp may read from. Generates dependencies between
   * this PGOp and both the last and next PGOp to write to each bit. No
   * dependencies exist when both PGOps just read from the same bit.
   */
  virtual bit_vector_t read_bits() const;

  /**
   * The classical bits this PGOp may write to. Generates dependencies between
   * this PGOp and both the last and next PGOp to read or write to each bit.
   */
  virtual bit_vector_t write_bits() const;

  virtual ~PGOp();

 protected:
  /**
   * Protected constructor subclasses can invoke to set type.
   */
  PGOp(PGOpType type);

  /**
   * An indicator of the subclass used for each object.
   */
  const PGOpType type_;
};

/**
 * PGOp for PGOpType::Rotation, representing a conventional Pauli gadget
 * (exponentiating a Pauli string).
 *
 * Whilst SpSymPauliTensor would completely capture both the string and angle,
 * the generic PGOp interface forces us to split it into a SpPauliStabiliser and
 * a separate angle.
 */
class PGRotation : public PGOp {
 public:
  /**
   * Get the Pauli string about which the rotation occurs. The phase of the
   * coefficient determines the direction of rotation.
   *
   * A const alias for PGRotation::port(0).
   */
  const SpPauliStabiliser& get_tensor() const;

  /**
   * Get the angle of rotation in half-turns.
   */
  const Expr& get_angle() const;

  /**
   * Constructs a rotation corresponding to exp(-i * \p tensor * \p angle *
   * pi/2).
   */
  PGRotation(const SpPauliStabiliser& tensor, const Expr& angle);

  // Overrides from PGOp
  virtual SymSet free_symbols() const override;
  virtual PGOp_ptr symbol_substitution(
      const SymEngine::map_basic_basic& sub_map) const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual std::vector<SpPauliStabiliser> active_paulis() const override;
  virtual SpPauliStabiliser& port(unsigned p) override;

 protected:
  SpPauliStabiliser tensor_;
  Expr angle_;
};

/**
 * PGOp for PGOpType::CliffordRot, representing a Clifford-angled Pauli gadget.
 * The angle of rotation is an integer number of quarter turns.
 */
class PGCliffordRot : public PGOp {
 public:
  /**
   * Get the Pauli string about which the rotation occurs. The phase of the
   * coefficient determines the direction of rotation.
   *
   * A const alias for PGCliffordRot::port(0).
   */
  const SpPauliStabiliser& get_tensor() const;

  /**
   * Get the angle of rotation as an integer number of quarter turns.
   */
  unsigned get_angle() const;

  /**
   * Constructs a Clifford-angled rotation corresponding to exp(-i * \p tensor *
   * \p angle * pi/4).
   */
  PGCliffordRot(const SpPauliStabiliser& tensor, unsigned angle);

  // Overrides from PGOp
  virtual SymSet free_symbols() const override;
  virtual PGOp_ptr symbol_substitution(
      const SymEngine::map_basic_basic& sub_map) const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual std::vector<SpPauliStabiliser> active_paulis() const override;
  virtual SpPauliStabiliser& port(unsigned p) override;

 protected:
  SpPauliStabiliser tensor_;
  unsigned angle_;
};

/**
 * PGOp for PGOpType::Measure, representing a non-destructive measurement of a
 * Pauli observable, writing the result to a given classical bit.
 */
class PGMeasure : public PGOp {
 public:
  /**
   * Get the Pauli observable being measured. The phase of the coefficient
   * determines whether the outcome of the measurement is flipped (i.e. the
   * expected measurement value directly gives the expectation value wrt the
   * phaseful Pauli observable).
   *
   * A const alias for PGMeasure::port(0).
   */
  const SpPauliStabiliser& get_tensor() const;

  /**
   * Get the classical bit to which the measurement result is written.
   */
  const Bit& get_target() const;

  /**
   * Constructs a non-destructive measurement of the phaseful Pauli observable
   * \p tensor which writes the outcome to \p target
   */
  PGMeasure(const SpPauliStabiliser& tensor, const Bit& target);

  // Overrides from PGOp
  virtual SymSet free_symbols() const override;
  virtual PGOp_ptr symbol_substitution(
      const SymEngine::map_basic_basic& sub_map) const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual std::vector<SpPauliStabiliser> active_paulis() const override;
  virtual SpPauliStabiliser& port(unsigned p) override;
  virtual bit_vector_t write_bits() const override;

 protected:
  SpPauliStabiliser tensor_;
  Bit target_;
};

/**
 * PGOp for PGOpType::Decoherence, representing a non-destructive measurement of
 * a Pauli observable where the measurement result is not recorded (i.e. a
 * generalisation of OpTyp::Collapse to an arbitrary Pauli basis).
 */
class PGDecoherence : public PGOp {
 public:
  /**
   * Get the Pauli observable being measured. Since the measurement result is
   * not recorded, the coefficient is irrelevant. This destroys information in
   * any Pauli basis for an anticommuting Pauli tensor.
   *
   * A const alias for PGDecoherence::port(0).
   */
  const SpPauliStabiliser& get_tensor() const;

  /**
   * Constructs a non-destructive measurement of the Pauli observable \p tensor
   * where the outcome is ignored.
   */
  PGDecoherence(const SpPauliStabiliser& tensor);

  // Overrides from PGOp
  virtual SymSet free_symbols() const override;
  virtual PGOp_ptr symbol_substitution(
      const SymEngine::map_basic_basic& sub_map) const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual std::vector<SpPauliStabiliser> active_paulis() const override;
  virtual SpPauliStabiliser& port(unsigned p) override;

 protected:
  SpPauliStabiliser tensor_;
};

/**
 * PGOp for PGOpType::Reset, representing a qubit reset operation (discard and
 * preparation of |0>) conjugated by a Clifford circuit.
 */
class PGReset : public PGOp {
 public:
  /**
   * Get the (phaseful) stabiliser guaranteed by the initialisation of the
   * reset. E.g. a regular reset operation without any Clifford conjugation
   * would guarantee +Z as a stabiliser.
   *
   * A const alias for PGReset::port(0).
   */
  const SpPauliStabiliser& get_stab() const;

  /**
   * Get the (phaseless) destabiliser, i.e. the additional Pauli basis in which
   * information is lost. E.g. a regular reset operation without any Clifford
   * conjugation would remove information in Z (see stab_), as well as X and Y;
   * we may choose either for destab_ as they relate by multiplication by stab_
   * so represent the same operation.
   *
   * A const alias for PGReset::port(1).
   */
  const SpPauliStabiliser& get_destab() const;

  /**
   * Construct a reset operation which removes information in the space spanned
   * by \p stab and \p destab and then instantiates a state to generate \p stab
   * as a stabiliser.
   */
  PGReset(const SpPauliStabiliser& stab, const SpPauliStabiliser& destab);

  // Overrides from PGOp
  virtual SymSet free_symbols() const override;
  virtual PGOp_ptr symbol_substitution(
      const SymEngine::map_basic_basic& sub_map) const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual unsigned n_paulis() const override;
  virtual std::vector<SpPauliStabiliser> active_paulis() const override;
  virtual SpPauliStabiliser& port(unsigned p) override;

 protected:
  SpPauliStabiliser stab_;
  SpPauliStabiliser destab_;
};

/**
 * PGOp for PGOpType::Conditional, wrapping another PGOp and executing it
 * conditional on the state of some classical bits.
 *
 * active_paulis and port defer to the inner op, and the condition bits are
 * added to the end of read_bits.
 */
class PGConditional : public PGOp {
 public:
  /**
   * Get the inner PGOp which is executed if the condition is met.
   */
  PGOp_ptr get_inner_op() const;

  /**
   * Get the classical bits that are checked for the condition.
   */
  const bit_vector_t& get_args() const;

  /**
   * Get the target value the bits need to be in order to execute the inner op.
   */
  unsigned get_value() const;

  /**
   * Construct a conditional operation, executing \p inner if the value of the
   * classical bits \p args is exactly \p value (using a little-endian format,
   * e.g. value 2 (10b) means args[0] must be 0 and args[1] must be 1).
   */
  PGConditional(PGOp_ptr inner, const bit_vector_t& args, unsigned value);

  // Overrides from PGOp
  virtual SymSet free_symbols() const override;
  virtual PGOp_ptr symbol_substitution(
      const SymEngine::map_basic_basic& sub_map) const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual unsigned n_paulis() const override;
  virtual std::vector<SpPauliStabiliser> active_paulis() const override;
  virtual SpPauliStabiliser& port(unsigned p) override;
  virtual bit_vector_t read_bits() const override;
  virtual bit_vector_t write_bits() const override;

 protected:
  PGOp_ptr inner_;
  bit_vector_t args_;
  unsigned value_;
};

/**
 * PGOp for PGOpType::StabAssert, representing a StabiliserAssertionBox,
 * possibly conjugated by a Clifford circuit. A pair of PauliTensors specify the
 * space mapped onto a single qubit to be used as an ancilla - this is reset and
 * the measurement encoded onto it. The result is written to a target bit before
 * the inverse Clifford circuit is applied.
 */
class PGStabAssertion : public PGOp {
 public:
  /**
   * Get the (phaseful) Pauli operator measured by the assertion. Success of the
   * assertion will leave this as a stabiliser of the final state.
   *
   * A const alias for PGStabAssertion::port(0).
   */
  const SpPauliStabiliser& get_stab() const;

  /**
   * Get the (phaseful) Pauli operator mapped into +Z on the ancilla qubit.
   * Success of the assertion will leave this as a stabiliser of the final
   * state.
   *
   * A const alias for PGStabAssertion::port(1).
   */
  const SpPauliStabiliser& get_anc_z() const;

  /**
   * Get the (phaseless) destabiliser wrt the measurement, i.e. a Pauli operator
   * which, along with anc_z_, generates the subspace on which information is
   * lost by the ancilla qubit reset. This is the operator which the conjugating
   * Clifford circuit maps to +X on the ancilla qubit.
   *
   * A const alias for PGStabAssertion::port(2).
   */
  const SpPauliStabiliser& get_anc_x() const;

  /**
   * Get the classical bit to which the measurement outcome is written.
   */
  const Bit& get_target() const;

  /**
   * Construct a stabiliser assertion, reducing the space spanned by \p anc_z
   * and \p anc_x onto a single qubit which is reset (the ancilla for the
   * assertion), then \p stab is loaded onto the ancilla before it is measured
   * and recorded in \p target and the ancilla mapped back into \p anc_z (adding
   * this as a stabiliser on a success) and \p anc_x
   */
  PGStabAssertion(
      const SpPauliStabiliser& stab, const SpPauliStabiliser& anc_z,
      const SpPauliStabiliser& anc_x, const Bit& target);

  // Overrides from PGOp
  virtual SymSet free_symbols() const override;
  virtual PGOp_ptr symbol_substitution(
      const SymEngine::map_basic_basic& sub_map) const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual unsigned n_paulis() const override;
  virtual std::vector<SpPauliStabiliser> active_paulis() const override;
  virtual SpPauliStabiliser& port(unsigned p) override;

 protected:
  SpPauliStabiliser stab_;
  SpPauliStabiliser anc_z_;
  SpPauliStabiliser anc_x_;
  Bit target_;
};

/**
 * PGOp for PGOpType::InputTableau. There should be at most one of these within
 * a PauliGraph, occurring at the start. This represents some ChoiMixTableau at
 * the start of the circuit, describing how any free inputs are mapped into the
 * space for the interior of the PauliGraph and any stabilisers generated by
 * initialisations. The active_paulis are the substrings over the output segment
 * (i.e. the segment relating to the interior of the PauliGraph).
 */
class PGInputTableau : public PGOp {
 public:
  /**
   * Get the tensor of row \p p as from the tableau; first component is for the
   * input segment, second for the output component (the active paulis); RxS
   * means SCR = C.
   */
  const ChoiMixTableau::row_tensor_t& get_full_row(unsigned p) const;

  /**
   * Constructs an input tableau operation from the given tableau.
   */
  PGInputTableau(const ChoiMixTableau& tableau);

  // Overrides from PGOp
  virtual SymSet free_symbols() const override;
  virtual PGOp_ptr symbol_substitution(
      const SymEngine::map_basic_basic& sub_map) const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual unsigned n_paulis() const override;
  virtual std::vector<SpPauliStabiliser> active_paulis() const override;
  virtual SpPauliStabiliser& port(unsigned p) override;

 protected:
  /**
   * Store the rows as SpPauliStabilisers rather than an actual tableau object
   * for easier modification of individual rows in the same way as for rewriting
   * on other PGOps. Specific rewrites making use of the input space (i.e.
   * contextual optimisations making use of initialisations) may wish to convert
   * this back into a tableau to make use of row combinations easier.
   */
  std::vector<ChoiMixTableau::row_tensor_t> rows_;
};

/**
 * PGOp for PGOpType::OutputTableau (dual to PGInputTableau). There should be at
 * most one of these within a PauliGraph, occurring at the end. This represents
 * some ChoiMixTableau at the end of the circuit, describing how Pauli operators
 * in the interior of the PauliGraph are mapped into the output space, and which
 * ones are post-selected or discarded. The active_paulis are the substrings
 * over the input segment (i.e. the segment relating to the interior of the
 * PauliGraph).
 */
class PGOutputTableau : public PGOp {
 public:
  /**
   * Get the tensor of row \p p as from the tableau; first component is for the
   * input segment (the active paulis), second for the output component; RxS
   * means SCR = C.
   */
  const ChoiMixTableau::row_tensor_t& get_full_row(unsigned p) const;

  /**
   * Constructs an output tableau operation from the given tableau.
   */
  PGOutputTableau(const ChoiMixTableau& tableau);

  // Overrides from PGOp
  virtual SymSet free_symbols() const override;
  virtual PGOp_ptr symbol_substitution(
      const SymEngine::map_basic_basic& sub_map) const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual unsigned n_paulis() const override;
  virtual std::vector<SpPauliStabiliser> active_paulis() const override;
  virtual SpPauliStabiliser& port(unsigned p) override;

 protected:
  /**
   * Store the rows as SpPauliStabilisers rather than an actual tableau object
   * for easier modification of individual rows in the same way as for rewriting
   * on other PGOps. Specific rewrites making use of the output space (i.e.
   * contextual optimisations making use of post-selections or discards) may
   * wish to convert this back into a tableau to make use of row combinations
   * easier.
   */
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
 * the interior Pauli strings in the anticommutation matrix. If they are not
 * provided explicitly, they are assumed to be identity circuits.
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
  /**
   * Construct an empty PauliGraph with no Qubits or Bits.
   */
  explicit PauliGraph();

  /**
   * Construct an empty PauliGraph representing the identity over some defined
   * set of Qubits and Bits. This will initially lack any PGInputTableau or
   * PGOutputTableau, so these should be added explicitly if they wish to be
   * used.
   */
  explicit PauliGraph(
      const std::set<Qubit>& qubits, const std::set<Bit>& bits = {});

  /**
   * Writes a graphviz representation of the PauliGraph to a stream. Use this
   * for visualisation. Each vertex in the PauliGraph is represented as a
   * cluster of graphviz vertices (one per active Pauli). Classical dependencies
   * are drawn as edges between clusters and the anti-commutation dependencies
   * between Paulis are drawn as edges between the corresponding vertices.
   */
  void to_graphviz(std::ostream& out) const;

  /**
   * Inserts a new vertex at the end of the PauliGraph. Throws an exception if a
   * PGInitialTableau is inserted after other vertices or if any vertex is
   * inserted after a PGOutputTableau.
   */
  PGVert add_vertex_at_end(PGOp_ptr op);

  /**
   * Verification of validity of the data structure. This is computationally
   * expensive so it is intended for use in debugging and tests, but not live
   * code.
   */
  void verify() const;

  /**
   * Returns all PGOps in a valid topological sort of the diagram. The exact
   * order depends on the internal order of vertices in c_graph_.
   */
  std::list<PGOp_ptr> pgop_sequence() const;

  friend PauliGraph tket::circuit_to_pauli_graph3(const tket::Circuit& circ);
  friend tket::Circuit tket::pauli_graph3_to_circuit_individual(
      const PauliGraph& pg, CXConfigType cx_config);

 private:
  MatrixXb pauli_ac_;
  PGIndex pauli_index_;
  PGClassicalGraph c_graph_;
  std::set<Qubit> qubits_;
  std::set<Bit> bits_;

  // Helper variables for tracking previous reads from and writes to each bit to
  // simplify adding dependencies in add_vertex_at_end.
  std::map<Bit, PGVert> last_writes_;
  std::map<Bit, std::unordered_set<PGVert>> last_reads_;

  std::optional<PGVert> initial_tableau_;
  std::optional<PGVert> final_tableau_;

  /**
   * Replaces the QubitPauliString of row \p target_r with i^{ \p coeff } *
   * source * target and updates pauli_ac_ accordingly.
   */
  void multiply_strings(
      unsigned source_r, unsigned target_r, quarter_turns_t coeff = 0);
};

}  // namespace pg
}  // namespace tket
