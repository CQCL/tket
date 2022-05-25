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

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <memory>

#include "Ops/Op.hpp"
#include "Utils/BiMapHeaders.hpp"
#include "Utils/EigenConfig.hpp"
#include "Utils/Json.hpp"
#include "Utils/MatrixAnalysis.hpp"
#include "Utils/UnitID.hpp"

namespace tket {

class Circuit;

/**
 * Abstract class for an operation from which a circuit can be extracted
 */
class Box : public Op {
 public:
  explicit Box(const OpType &type, const op_signature_t &signature = {})
      : Op(type), signature_(signature), circ_(), id_(idgen()) {
    if (!is_box_type(type)) throw NotValid();
  }

  Box(const Box &other)
      : Op(other.get_type()),
        signature_(other.signature_),
        circ_(other.circ_),
        id_(other.id_) {}

  /** Number of Quantum inputs */
  unsigned n_qubits() const override;

  /** Number of Boolean inputs */
  unsigned n_boolean() const;

  /** Number of Classical inputs */
  unsigned n_classical() const;

  op_signature_t get_signature() const override;

  nlohmann::json serialize() const override;

  static Op_ptr deserialize(const nlohmann::json &j);

  /** Circuit represented by box */
  std::shared_ptr<Circuit> to_circuit() const {
    if (circ_ == nullptr) generate_circuit();
    return circ_;
  };

  /** Unique identifier (preserved on copy) */
  boost::uuids::uuid get_id() const { return id_; }

  template <typename BoxT>
  friend Op_ptr set_box_id(BoxT &b, boost::uuids::uuid newid);

 protected:
  static boost::uuids::uuid idgen() {
    static boost::uuids::random_generator gen = {};

    return gen();
  }
  op_signature_t signature_;
  mutable std::shared_ptr<Circuit> circ_;
  boost::uuids::uuid id_;

  virtual void generate_circuit() const = 0;
};

// json for base Box attributes
nlohmann::json core_box_json(const Box &box);

/**
 * @brief Set explicit ID on a box
 *
 * This is used for deserialization.
 *
 * @param b box
 * @param[in] newid new ID
 *
 * @tparam BoxT concrete box type
 *
 * @return box with desired ID
 */
template <typename BoxT>
Op_ptr set_box_id(BoxT &b, boost::uuids::uuid newid) {
  b.id_ = newid;
  return std::make_shared<BoxT>(b);
}

/**
 * Operation defined as a circuit.
 */
class CircBox : public Box {
 public:
  /**
   * Construct from a given circuit. The circuit must be simple.
   */
  explicit CircBox(const Circuit &circ);

  /**
   * Construct from the empty circuit
   */
  CircBox();

  /**
   * Copy constructor
   */
  CircBox(const CircBox &other);

  ~CircBox() override {}

  bool is_clifford() const override;

  Op_ptr symbol_substitution(
      const SymEngine::map_basic_basic &sub_map) const override;

  SymSet free_symbols() const override;

  /**
   * Equality check between two CircBox instances
   */
  bool is_equal(const Op &op_other) const override {
    const CircBox &other = dynamic_cast<const CircBox &>(op_other);
    return id_ == other.get_id();
  }

  Op_ptr dagger() const override;

  Op_ptr transpose() const override;

  static Op_ptr from_json(const nlohmann::json &j);

  static nlohmann::json to_json(const Op_ptr &op);

 protected:
  void generate_circuit() const override {}  // Already set by constructor
};

/**
 * One-qubit operation defined as a unitary matrix
 */
class Unitary1qBox : public Box {
 public:
  /**
   * Construct from a given 2x2 unitary matrix
   *
   * @param m unitary matrix
   */
  explicit Unitary1qBox(const Eigen::Matrix2cd &m);

  /**
   * Construct from the identity matrix
   */
  Unitary1qBox();

  /**
   * Copy constructor
   */
  Unitary1qBox(const Unitary1qBox &other);

  ~Unitary1qBox() override {}

  Op_ptr symbol_substitution(
      const SymEngine::map_basic_basic &) const override {
    return Op_ptr();
  }

  SymSet free_symbols() const override { return {}; }

  /**
   * Equality check between two Unitary1qBox instances
   */
  bool is_equal(const Op &op_other) const override {
    const Unitary1qBox &other = dynamic_cast<const Unitary1qBox &>(op_other);
    return id_ == other.get_id();
  }

  /** Get the unitary matrix correspnding to this operation */
  Eigen::Matrix2cd get_matrix() const { return m_; }

  Eigen::MatrixXcd get_unitary() const override { return m_; }

  Op_ptr dagger() const override;

  Op_ptr transpose() const override;

  bool is_clifford() const override;

  static Op_ptr from_json(const nlohmann::json &j);

  static nlohmann::json to_json(const Op_ptr &op);

 protected:
  void generate_circuit() const override;

 private:
  const Eigen::Matrix2cd m_;
};

/**
 * Two-qubit operation defined as a unitary matrix (ILO-BE)
 */
class Unitary2qBox : public Box {
 public:
  /**
   * Construct from a given 4x4 unitary matrix
   *
   * @param m unitary matrix
   * @param basis basis order convention for matrix
   */
  explicit Unitary2qBox(
      const Eigen::Matrix4cd &m, BasisOrder basis = BasisOrder::ilo);

  /**
   * Construct from the identity matrix
   */
  Unitary2qBox();

  /**
   * Copy constructor
   */
  Unitary2qBox(const Unitary2qBox &other);

  ~Unitary2qBox() override {}

  Op_ptr symbol_substitution(
      const SymEngine::map_basic_basic &) const override {
    return Op_ptr();
  }

  SymSet free_symbols() const override { return {}; }

  /**
   * Equality check between two Unitary2qBox instances
   */
  bool is_equal(const Op &op_other) const override {
    const Unitary2qBox &other = dynamic_cast<const Unitary2qBox &>(op_other);
    return id_ == other.get_id();
  }

  /** Get the unitary matrix correspnding to this operation */
  Eigen::Matrix4cd get_matrix() const { return m_; }

  Eigen::MatrixXcd get_unitary() const override { return m_; }

  Op_ptr dagger() const override;

  Op_ptr transpose() const override;

  static Op_ptr from_json(const nlohmann::json &j);

  static nlohmann::json to_json(const Op_ptr &op);

 protected:
  void generate_circuit() const override;

 private:
  const Eigen::Matrix4cd m_;
};

/**
 * Three-qubit operation defined as a unitary matrix (ILO-BE)
 */
class Unitary3qBox : public Box {
 public:
  /**
   * Construct from a given 8x8 unitary matrix
   *
   * @param m unitary matrix
   * @param basis basis order convention for matrix
   */
  explicit Unitary3qBox(const Matrix8cd &m, BasisOrder basis = BasisOrder::ilo);

  /**
   * Construct from the identity matrix
   */
  Unitary3qBox();

  /**
   * Copy constructor
   */
  Unitary3qBox(const Unitary3qBox &other);

  ~Unitary3qBox() override {}

  Op_ptr symbol_substitution(
      const SymEngine::map_basic_basic &) const override {
    return Op_ptr();
  }

  SymSet free_symbols() const override { return {}; }

  /**
   * Equality check between two Unitary3qBox instances
   */
  bool is_equal(const Op &op_other) const override {
    const Unitary3qBox &other = dynamic_cast<const Unitary3qBox &>(op_other);
    return id_ == other.get_id();
  }

  /** Get the unitary matrix correspnding to this operation */
  Matrix8cd get_matrix() const { return m_; }

  Eigen::MatrixXcd get_unitary() const override { return m_; }

  Op_ptr dagger() const override;

  Op_ptr transpose() const override;

  static Op_ptr from_json(const nlohmann::json &j);

  static nlohmann::json to_json(const Op_ptr &op);

 protected:
  void generate_circuit() const override;

 private:
  const Matrix8cd m_;
};

/**
 * Two-qubit operation defined in terms of a hermitian matrix and a phase.
 *
 * The unitary corresponding to the matrix A and phase t is exp(i*t*A).
 * Matrix A is stored in ILO-BE form.
 */
class ExpBox : public Box {
 public:
  /**
   * Construct from a given 4x4 hermitian matrix and optional phase.
   *
   * @param A hermitian matrix
   * @param t phase to apply
   * @param basis basis order convention for matrix
   */
  ExpBox(
      const Eigen::Matrix4cd &A, double t = 1.,
      BasisOrder basis = BasisOrder::ilo);

  /**
   * Construct from the zero matrix (resulting in the identity)
   */
  ExpBox();

  /**
   * Copy constructor
   */
  ExpBox(const ExpBox &other);

  ~ExpBox() override {}

  Op_ptr symbol_substitution(
      const SymEngine::map_basic_basic &) const override {
    return Op_ptr();
  }

  SymSet free_symbols() const override { return {}; }

  /**
   * Equality check between two ExpBox instances
   */
  bool is_equal(const Op &op_other) const override {
    const ExpBox &other = dynamic_cast<const ExpBox &>(op_other);
    return id_ == other.get_id();
  }

  /** Get the hermitian matrix and phase parameter */
  std::pair<Eigen::Matrix4cd, double> get_matrix_and_phase() const {
    return std::make_pair(A_, t_);
  }

  Op_ptr dagger() const override;

  Op_ptr transpose() const override;

  static Op_ptr from_json(const nlohmann::json &j);

  static nlohmann::json to_json(const Op_ptr &op);

 protected:
  void generate_circuit() const override;

 private:
  const Eigen::Matrix4cd A_;
  double t_;
};

/**
 * Operation defined as the exponential of a tensor of Pauli operators
 */
class PauliExpBox : public Box {
 public:
  /**
   * The operation implements the unitary operator
   * \f$ e^{-\frac12 i \pi t \sigma_0 \otimes \sigma_1 \otimes \cdots} \f$
   * where \f$ \sigma_i \in \{I,X,Y,Z\} \f$ are the Pauli operators.
   */
  PauliExpBox(const std::vector<Pauli> &paulis, const Expr &t);

  /**
   * Construct from the empty vector
   */
  PauliExpBox();

  /**
   * Copy constructor
   */
  PauliExpBox(const PauliExpBox &other);

  ~PauliExpBox() override {}

  bool is_clifford() const override;

  SymSet free_symbols() const override;

  /**
   * Equality check between two PauliExpBox instances
   */
  bool is_equal(const Op &op_other) const override {
    const PauliExpBox &other = dynamic_cast<const PauliExpBox &>(op_other);
    return id_ == other.get_id();
  }

  /** Get the Pauli string */
  std::vector<Pauli> get_paulis() const { return paulis_; }

  /** Get the phase parameter */
  Expr get_phase() const { return t_; }

  Op_ptr dagger() const override;

  Op_ptr transpose() const override;

  Op_ptr symbol_substitution(
      const SymEngine::map_basic_basic &sub_map) const override;

  static Op_ptr from_json(const nlohmann::json &j);

  static nlohmann::json to_json(const Op_ptr &op);

 protected:
  void generate_circuit() const override;

 private:
  std::vector<Pauli> paulis_;
  Expr t_;
};

class CompositeGateDef;
typedef std::shared_ptr<CompositeGateDef> composite_def_ptr_t;

// CompositeGateDef
JSON_DECL(composite_def_ptr_t)

class CompositeGateDef : public std::enable_shared_from_this<CompositeGateDef> {
 public:
  CompositeGateDef(
      const std::string &name, const Circuit &def,
      const std::vector<Sym> &args);

  static composite_def_ptr_t define_gate(
      const std::string &name, const Circuit &def,
      const std::vector<Sym> &args);

  Circuit instance(const std::vector<Expr> &params) const;

  std::string get_name() const { return name_; }
  std::vector<Sym> get_args() const { return args_; }
  std::shared_ptr<Circuit> get_def() const { return def_; }
  unsigned n_args() const { return args_.size(); }
  op_signature_t signature() const;

  bool operator==(const CompositeGateDef &other) const;

 private:
  std::string name_;
  std::shared_ptr<Circuit> def_;
  std::vector<Sym> args_;

  CompositeGateDef() {}
};

class CustomGate : public Box {
 public:
  CustomGate(const composite_def_ptr_t &gate, const std::vector<Expr> &params);
  CustomGate(const CustomGate &other);

  SymSet free_symbols() const override;

  composite_def_ptr_t get_gate() const { return gate_; }
  std::vector<Expr> get_params() const override { return params_; }
  std::string get_name(bool latex = false) const override;

  Op_ptr symbol_substitution(
      const SymEngine::map_basic_basic &sub_map) const override;

  /**
   * Equality check between two CustomGate instances.
   * This does more than simply checking id_.
   */
  bool is_equal(const Op &op_other) const override;

  static Op_ptr from_json(const nlohmann::json &j);

  static nlohmann::json to_json(const Op_ptr &op);

  bool is_clifford() const override;

 protected:
  void generate_circuit() const override;
  CustomGate() : Box(OpType::CustomGate), gate_(), params_() {}

 private:
  composite_def_ptr_t gate_;
  const std::vector<Expr> params_;
};

/**
 * Wraps another quantum op, adding control qubits
 */
class QControlBox : public Box {
 public:
  /**
   * Construct from a given op and number of controls
   *
   * @param op op to control
   * @param n_controls number of qubit controls to add
   */
  explicit QControlBox(const Op_ptr &op, unsigned n_controls = 1);

  /**
   * Copy constructor
   */
  QControlBox(const QControlBox &other);

  ~QControlBox() override {}

  Op_ptr symbol_substitution(
      const SymEngine::map_basic_basic &sub_map) const override;

  SymSet free_symbols() const override;

  /**
   * Equality check between two QControlBox instances
   */
  bool is_equal(const Op &op_other) const override {
    const QControlBox &other = dynamic_cast<const QControlBox &>(op_other);
    return id_ == other.get_id();
  }

  std::string get_command_str(const unit_vector_t &args) const override;

  Op_ptr dagger() const override;

  Op_ptr transpose() const override;

  Op_ptr get_op() const { return op_; }
  unsigned get_n_controls() const { return n_controls_; }

  static Op_ptr from_json(const nlohmann::json &j);

  static nlohmann::json to_json(const Op_ptr &op);

 protected:
  void generate_circuit() const override;
  QControlBox()
      : Box(OpType::QControlBox), op_(), n_controls_(0), n_inner_qubits_(0) {}

 private:
  const Op_ptr op_;
  const unsigned n_controls_;
  unsigned n_inner_qubits_;
};

class ProjectorAssertionBox : public Box {
 public:
  /**
   * Construct from a given 2x2, 4x4 or 8x8 projector matrix
   *
   * @param m projector matrix
   * @param basis basis order convention for matrix
   */
  explicit ProjectorAssertionBox(
      const Eigen::MatrixXcd &m, BasisOrder basis = BasisOrder::ilo);

  /**
   * Copy constructor
   */
  ProjectorAssertionBox(const ProjectorAssertionBox &other);

  ~ProjectorAssertionBox() override {}

  Op_ptr symbol_substitution(
      const SymEngine::map_basic_basic &) const override {
    return Op_ptr();
  }

  SymSet free_symbols() const override { return {}; }

  /**
   * Equality check between two ProjectorAssertionBox instances
   */
  bool is_equal(const Op &op_other) const override {
    const ProjectorAssertionBox &other =
        dynamic_cast<const ProjectorAssertionBox &>(op_other);
    return id_ == other.get_id();
  }

  /** Get the unitary matrix correspnding to this operation */
  Eigen::MatrixXcd get_matrix() const { return m_; }
  std::vector<bool> get_expected_readouts() const { return expected_readouts_; }

  Op_ptr dagger() const override;

  Op_ptr transpose() const override;

  op_signature_t get_signature() const override;

  static Op_ptr from_json(const nlohmann::json &j);

  static nlohmann::json to_json(const Op_ptr &op);

 protected:
  void generate_circuit() const override;

 private:
  const Eigen::MatrixXcd m_;
  // expected readouts the debug bits
  // false -> 0
  // true -> 1
  mutable std::vector<bool> expected_readouts_;
};

class StabiliserAssertionBox : public Box {
 public:
  /**
   * Construct from a set of stabiliser Pauli strings
   *
   * @param paulis a set of stabiliser Pauli strings
   */
  explicit StabiliserAssertionBox(const PauliStabiliserList &paulis);

  /**
   * Copy constructor
   */
  StabiliserAssertionBox(const StabiliserAssertionBox &other);

  ~StabiliserAssertionBox() override {}
  Op_ptr symbol_substitution(
      const SymEngine::map_basic_basic &) const override {
    return Op_ptr();
  }

  SymSet free_symbols() const override { return {}; }

  /**
   * Equality check between two StabiliserAssertionBox instances
   */
  bool is_equal(const Op &op_other) const override {
    const StabiliserAssertionBox &other =
        dynamic_cast<const StabiliserAssertionBox &>(op_other);
    return id_ == other.get_id();
  }

  /** Get the pauli stabilisers */
  PauliStabiliserList get_stabilisers() const { return paulis_; }
  std::vector<bool> get_expected_readouts() const { return expected_readouts_; }

  Op_ptr dagger() const override;

  Op_ptr transpose() const override;

  static Op_ptr from_json(const nlohmann::json &j);

  static nlohmann::json to_json(const Op_ptr &op);

  op_signature_t get_signature() const override;

 protected:
  void generate_circuit() const override;

 private:
  const PauliStabiliserList paulis_;
  // expected readouts the debug bits
  // false -> 0
  // true -> 1
  mutable std::vector<bool> expected_readouts_;
};

}  // namespace tket
