// Copyright 2019-2024 Cambridge Quantum Computing
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
#include "tket/Utils/PauliTensor.hpp"

namespace tket {
namespace pg {

// Forwards declarations for internal uses
class PGOp;
// Not const as we wish for these to be updated in-place
typedef std::shared_ptr<PGOp> PGOp_ptr;

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

  // Some other PGOp conditioned on a quantum state
  QControl,

  // A collection of tensors of opaque boxed circuit components, conditioned on
  // different values of a quantum state
  MultiplexedTensoredBox,

  // A collection of rotations in the same basis, conditioned on different
  // values of a quantum state
  MultiplexedRotation,

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
 * The active Paulis of each PGOp can be split into pairs of anti-commuting
 * Pauli strings (reducible to the space of one qubit) and additional Pauli
 * strings that commute with all others (reducible to a qubit with a commuting
 * Pauli operator).
 *
 * This signature indicates the number of qubits used to implement the Op after
 * diagonalisation: one per anti-commuting pair, plus one per additional
 * commuting operator. The PGOp is a valid target for GraySynth when there is
 * exactly one commuting operator which becomes the target "phase", with each
 * anti-commuting pair just acting as an ancilla which comes into play when the
 * Op is ready to be synthesised.
 */
struct PGOp_signature {
  // Pairs of anti-commuting Pauli strings
  std::list<std::pair<SpPauliStabiliser, SpPauliStabiliser>> anti_comm_pairs;
  // Pauli strings which commute with all others within the PGOp
  std::list<SpPauliStabiliser> comm_set;
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
   * Deep copy operation, since PGOp_ptr does not point to a const PGOp
   */
  virtual PGOp_ptr clone() const = 0;

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
   * means it is unlikely they commute). Checks whether all active Paulis
   * mutually commute between the two PGOps.
   */
  bool commutes_with(const PGOp& other) const;

  /**
   * Returns the number of active Paulis, i.e. a measure of the size of the
   * subspace of the Pauli group on which this operator acts non-trivially.
   */
  virtual unsigned n_paulis() const;

  /**
   * Returns a collection of Pauli operators dictating the subspace on which the
   * op acts non-trivially. The guarantee is that, if another op commutes with
   * all active Pauli operators, then it commutes with the PGOp (the
   * converse need not hold, for example Rotation gates with angle 0).
   *
   * SpPauliStabiliser is used to account for phase information in common
   * updates and rewrites (e.g. Clifford reordering rules). Some PGOpTypes won't
   * be phase-sensitive (e.g. Decoherence) and some may double-up on phase
   * information (e.g. CliffordRot(P,3) is the same as CliffordRot(-P,1)), but
   * having just +- phase info on the easily accessible PauliTensors is a
   * reasonable middle ground and the other cases can be easily handled on an
   * ad-hoc basis.
   *
   * This signature groups the active Pauli operators according to their
   * commutativity with each other.
   */
  virtual PGOp_signature pauli_signature() const = 0;

  /**
   * Gives direct reference access to each active Pauli as a SpPauliStabiliser
   * via an index into some fixed ordering set by the semantics of the subclass,
   * e.g. the projected stabiliser of a PGReset is the Pauli operator at port 0
   * and the lost stabiliser is at port 1.
   *
   * This is most useful to give immediate, generic access to the active Paulis
   * for rewrites and synthesis without having to inspect the PGOpType and cast
   * to the appropriate subclass.
   */
  SpPauliStabiliser& port(unsigned p);
  virtual const SpPauliStabiliser& port(unsigned p) const = 0;

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
  virtual PGOp_ptr clone() const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual PGOp_signature pauli_signature() const override;
  virtual const SpPauliStabiliser& port(unsigned p) const override;

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
  virtual PGOp_ptr clone() const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual PGOp_signature pauli_signature() const override;
  virtual const SpPauliStabiliser& port(unsigned p) const override;

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
  virtual PGOp_ptr clone() const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual PGOp_signature pauli_signature() const override;
  virtual const SpPauliStabiliser& port(unsigned p) const override;
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
  virtual PGOp_ptr clone() const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual PGOp_signature pauli_signature() const override;
  virtual const SpPauliStabiliser& port(unsigned p) const override;

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
  virtual PGOp_ptr clone() const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual unsigned n_paulis() const override;
  virtual PGOp_signature pauli_signature() const override;
  virtual const SpPauliStabiliser& port(unsigned p) const override;

 protected:
  SpPauliStabiliser stab_;
  SpPauliStabiliser destab_;
};

/**
 * PGOp for PGOpType::Conditional, wrapping another PGOp and executing it
 * conditional on the state of some classical bits.
 *
 * pauli_signature and port defer to the inner op, and the condition bits are
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
  virtual PGOp_ptr clone() const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual unsigned n_paulis() const override;
  virtual PGOp_signature pauli_signature() const override;
  virtual const SpPauliStabiliser& port(unsigned p) const override;
  virtual bit_vector_t read_bits() const override;
  virtual bit_vector_t write_bits() const override;

 protected:
  PGOp_ptr inner_;
  bit_vector_t args_;
  unsigned value_;
};

/**
 * PGOp for PGOpType::QControl, wrapping another (unitary) PGOp and executing it
 * conditional on the state of some qubits.
 *
 * The first ports give the paulis into which the control qubits are
 * encoded, followed by the active Paulis of the inner op.
 */
class PGQControl : public PGOp {
 public:
  /**
   * Get the inner PGOp which is executed coherently according to the control
   * qubits.
   */
  PGOp_ptr get_inner_op() const;

  /**
   * Get the Pauli strings into which the controls are encoded.
   */
  const std::vector<SpPauliStabiliser>& get_control_paulis() const;

  /**
   * Get the target value the control qubits need to be in order to execute the
   * inner op.
   */
  std::vector<bool> get_value() const;

  /**
   * Construct a quantum-controlled operation, executing \p inner coherently if
   * the value of the \p control_paulis is exactly \p value (e.g. value [false,
   * true] means we apply the inner op on states that are +1 (false) eigenstates
   * of control_paulis[0] and -1 (true) eigenstates of control_paulis[1])
   */
  PGQControl(
      PGOp_ptr inner, const std::vector<SpPauliStabiliser>& control_paulis,
      std::vector<bool> value);

  // Overrides from PGOp
  virtual SymSet free_symbols() const override;
  virtual PGOp_ptr symbol_substitution(
      const SymEngine::map_basic_basic& sub_map) const override;
  virtual PGOp_ptr clone() const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual unsigned n_paulis() const override;
  virtual PGOp_signature pauli_signature() const override;
  virtual const SpPauliStabiliser& port(unsigned p) const override;

 protected:
  PGOp_ptr inner_;
  std::vector<SpPauliStabiliser> control_paulis_;
  std::vector<bool> value_;
};

/**
 * PGOp for PGOpType::MultiplexedRotation, encapsulating rotations of different
 * angles in the same basis conditioned on different values of the state of some
 * qubits.
 *
 * The first ports give the paulis into which the control qubits are
 * encoded, followed by the pauli into which the target rotation is encoded.
 */
class PGMultiplexedRotation : public PGOp {
 public:
  /**
   * Get the map between values of the control qubits and the angle of rotation
   * (in half-turns) that is performed coherently at that value.
   */
  const std::map<std::vector<bool>, Expr>& get_angle_map() const;

  /**
   * Get the Pauli strings into which the controls are encoded.
   */
  const std::vector<SpPauliStabiliser>& get_control_paulis() const;

  /**
   * Get the Pauli string about which the target rotation is applied.
   */
  const SpPauliStabiliser& get_target_pauli() const;

  /**
   * Construct a multiplexed operation where, if the input state's eigenvalues
   * wrt \p control_paulis are the vector ``value`` (e.g. value [false, false,
   * true] means a +1 eigenvalue for control_paulis[0-1] and a -1 eigenvalue for
   * control_paulis[2]), then the rotation exp(-i * \p target_pauli * \p
   * angle_map [value] * pi/4) is applied.
   */
  PGMultiplexedRotation(
      const std::map<std::vector<bool>, Expr>& angle_map,
      const std::vector<SpPauliStabiliser>& control_paulis,
      const SpPauliStabiliser& target_pauli);

  // Overrides from PGOp
  virtual SymSet free_symbols() const override;
  virtual PGOp_ptr symbol_substitution(
      const SymEngine::map_basic_basic& sub_map) const override;
  virtual PGOp_ptr clone() const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual unsigned n_paulis() const override;
  virtual PGOp_signature pauli_signature() const override;
  virtual const SpPauliStabiliser& port(unsigned p) const override;

 protected:
  std::map<std::vector<bool>, Expr> angle_map_;
  std::vector<SpPauliStabiliser> control_paulis_;
  SpPauliStabiliser target_pauli_;
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
  virtual PGOp_ptr clone() const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual unsigned n_paulis() const override;
  virtual PGOp_signature pauli_signature() const override;
  virtual const SpPauliStabiliser& port(unsigned p) const override;
  virtual bit_vector_t write_bits() const override;

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
 * initialisations. The active Paulis are the substrings over the output segment
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
   * Combine all rows back into a ChoiMixTableau object for a complete view of
   * the process.
   */
  ChoiMixTableau to_cm_tableau() const;

  /**
   * Constructs an input tableau operation from the given tableau.
   */
  PGInputTableau(const ChoiMixTableau& tableau);

  // Overrides from PGOp
  virtual SymSet free_symbols() const override;
  virtual PGOp_ptr symbol_substitution(
      const SymEngine::map_basic_basic& sub_map) const override;
  virtual PGOp_ptr clone() const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual unsigned n_paulis() const override;
  // CAUTION: Paulis in signature may not match ports due to gaussian
  // elimination used in determining anti-commuting pairs
  virtual PGOp_signature pauli_signature() const override;
  virtual const SpPauliStabiliser& port(unsigned p) const override;

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
 * ones are post-selected or discarded. The active Paulis are the substrings
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
   * Combine all rows back into a ChoiMixTableau object for a complete view of
   * the process.
   */
  ChoiMixTableau to_cm_tableau() const;

  /**
   * Constructs an output tableau operation from the given tableau.
   */
  PGOutputTableau(const ChoiMixTableau& tableau);

  // Overrides from PGOp
  virtual SymSet free_symbols() const override;
  virtual PGOp_ptr symbol_substitution(
      const SymEngine::map_basic_basic& sub_map) const override;
  virtual PGOp_ptr clone() const override;
  virtual std::string get_name(bool latex = false) const override;
  virtual bool is_equal(const PGOp& other) const override;
  virtual unsigned n_paulis() const override;
  // CAUTION: Paulis in signature may not match ports due to gaussian
  // elimination used in determining anti-commuting pairs
  virtual PGOp_signature pauli_signature() const override;
  virtual const SpPauliStabiliser& port(unsigned p) const override;

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

}  // namespace pg
}  // namespace tket
