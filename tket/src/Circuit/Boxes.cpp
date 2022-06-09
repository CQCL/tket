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

#include "Boxes.hpp"

#include <memory>
#include <numeric>

#include "CircUtils.hpp"
#include "Circuit/AssertionSynthesis.hpp"
#include "Command.hpp"
#include "Gate/Rotation.hpp"
#include "Ops/OpJsonFactory.hpp"
#include "Ops/OpPtr.hpp"
#include "ThreeQubitConversion.hpp"
#include "Utils/Assert.hpp"
#include "Utils/EigenConfig.hpp"
#include "Utils/Expression.hpp"
#include "Utils/Json.hpp"
#include "Utils/PauliStrings.hpp"

namespace tket {

unsigned Box::n_qubits() const {
  op_signature_t sig = get_signature();
  return std::count(sig.begin(), sig.end(), EdgeType::Quantum);
}

unsigned Box::n_boolean() const {
  op_signature_t sig = get_signature();
  return std::count(sig.begin(), sig.end(), EdgeType::Boolean);
}

unsigned Box::n_classical() const {
  op_signature_t sig = get_signature();
  return std::count(sig.begin(), sig.end(), EdgeType::Classical);
}

op_signature_t Box::get_signature() const {
  std::optional<op_signature_t> sig = desc_.signature();
  if (sig)
    return *sig;
  else
    return signature_;
}

nlohmann::json Box::serialize() const {
  nlohmann::json j;
  j["type"] = get_type();
  j["box"] = OpJsonFactory::to_json(shared_from_this());
  return j;
}

Op_ptr Box::deserialize(const nlohmann::json &j) {
  return OpJsonFactory::from_json(j.at("box"));
}

CircBox::CircBox(const Circuit &circ) : Box(OpType::CircBox) {
  if (!circ.is_simple()) throw SimpleOnly();
  signature_ = op_signature_t(circ.n_qubits(), EdgeType::Quantum);
  op_signature_t bits(circ.n_bits(), EdgeType::Classical);
  signature_.insert(signature_.end(), bits.begin(), bits.end());
  circ_ = std::make_shared<Circuit>(circ);
}

CircBox::CircBox(const CircBox &other) : Box(other) {}

CircBox::CircBox() : Box(OpType::CircBox) {
  circ_ = std::make_shared<Circuit>();
}

bool CircBox::is_clifford() const {
  BGL_FORALL_VERTICES(v, circ_->dag, DAG) {
    if (!circ_->get_Op_ptr_from_Vertex(v)->is_clifford()) return false;
  }
  return true;
}

Op_ptr CircBox::symbol_substitution(
    const SymEngine::map_basic_basic &sub_map) const {
  Circuit new_circ(*to_circuit());
  new_circ.symbol_substitution(sub_map);
  return std::make_shared<CircBox>(new_circ);
}

SymSet CircBox::free_symbols() const { return to_circuit()->free_symbols(); }

Op_ptr CircBox::dagger() const {
  return std::make_shared<CircBox>(circ_->dagger());
}

Op_ptr CircBox::transpose() const {
  return std::make_shared<CircBox>(circ_->transpose());
}

Unitary1qBox::Unitary1qBox(const Eigen::Matrix2cd &m)
    : Box(OpType::Unitary1qBox), m_(m) {
  if (!is_unitary(m)) {
    throw CircuitInvalidity("Matrix for Unitary1qBox must be unitary");
  }
}

Unitary1qBox::Unitary1qBox(const Unitary1qBox &other)
    : Box(other), m_(other.m_) {}

Unitary1qBox::Unitary1qBox() : Unitary1qBox(Eigen::Matrix2cd::Identity()) {}

Op_ptr Unitary1qBox::dagger() const {
  return std::make_shared<Unitary1qBox>(m_.conjugate().transpose());
}

Op_ptr Unitary1qBox::transpose() const {
  return std::make_shared<Unitary1qBox>(m_.transpose());
}

bool Unitary1qBox::is_clifford() const {
  std::vector<Command> cmds = to_circuit()->get_commands();
  TKET_ASSERT(cmds.size() == 1);
  return cmds[0].get_op_ptr()->is_clifford();
}

void Unitary1qBox::generate_circuit() const {
  std::vector<double> tk1_params = tk1_angles_from_unitary(m_);
  Circuit temp_circ(1);
  temp_circ.add_op<unsigned>(
      OpType::TK1, {tk1_params[0], tk1_params[1], tk1_params[2]}, {0});
  circ_ = std::make_shared<Circuit>(temp_circ);
  circ_->add_phase(tk1_params[3]);
}

Unitary2qBox::Unitary2qBox(const Eigen::Matrix4cd &m, BasisOrder basis)
    : Box(OpType::Unitary2qBox),
      m_(basis == BasisOrder::ilo ? m : reverse_indexing(m)) {
  if (!is_unitary(m)) {
    throw CircuitInvalidity("Matrix for Unitary2qBox must be unitary");
  }
}

Unitary2qBox::Unitary2qBox(const Unitary2qBox &other)
    : Box(other), m_(other.m_) {}

Unitary2qBox::Unitary2qBox() : Unitary2qBox(Eigen::Matrix4cd::Identity()) {}

Op_ptr Unitary2qBox::dagger() const {
  return std::make_shared<Unitary2qBox>(m_.conjugate().transpose());
}

Op_ptr Unitary2qBox::transpose() const {
  return std::make_shared<Unitary2qBox>(m_.transpose());
}

void Unitary2qBox::generate_circuit() const {
  circ_ = std::make_shared<Circuit>(two_qubit_canonical(m_));
}

Unitary3qBox::Unitary3qBox(const Matrix8cd &m, BasisOrder basis)
    : Box(OpType::Unitary3qBox),
      m_(basis == BasisOrder::ilo ? m : reverse_indexing(m)) {}

Unitary3qBox::Unitary3qBox(const Unitary3qBox &other)
    : Box(other), m_(other.m_) {}

Unitary3qBox::Unitary3qBox() : Unitary3qBox(Matrix8cd::Identity()) {}

Op_ptr Unitary3qBox::dagger() const {
  return std::make_shared<Unitary3qBox>(m_.adjoint());
}

Op_ptr Unitary3qBox::transpose() const {
  return std::make_shared<Unitary3qBox>(m_.transpose());
}

void Unitary3qBox::generate_circuit() const {
  circ_ = std::make_shared<Circuit>(three_qubit_synthesis(m_));
}

ExpBox::ExpBox(const Eigen::Matrix4cd &A, double t, BasisOrder basis)
    : Box(OpType::ExpBox),
      A_(basis == BasisOrder::ilo ? A : reverse_indexing(A)),
      t_(t) {
  if (!A.isApprox(A.adjoint())) {
    throw CircuitInvalidity("Matrix for ExpBox must be Hermitian");
  }
}

ExpBox::ExpBox(const ExpBox &other) : Box(other), A_(other.A_), t_(other.t_) {}

ExpBox::ExpBox() : ExpBox(Eigen::Matrix4cd::Zero(), 1.) {}

Op_ptr ExpBox::dagger() const { return std::make_shared<ExpBox>(A_, -t_); }

Op_ptr ExpBox::transpose() const {
  return std::make_shared<ExpBox>(A_.transpose(), t_);
}

void ExpBox::generate_circuit() const {
  circ_ = std::make_shared<Circuit>(two_qubit_canonical((i_ * t_ * A_).exp()));
}

PauliExpBox::PauliExpBox(const std::vector<Pauli> &paulis, const Expr &t)
    : Box(OpType::PauliExpBox,
          op_signature_t(paulis.size(), EdgeType::Quantum)),
      paulis_(paulis),
      t_(t) {}

PauliExpBox::PauliExpBox(const PauliExpBox &other)
    : Box(other), paulis_(other.paulis_), t_(other.t_) {}

PauliExpBox::PauliExpBox() : PauliExpBox({}, 0.) {}

bool PauliExpBox::is_clifford() const {
  return equiv_0(4 * t_) || paulis_.empty();
}

SymSet PauliExpBox::free_symbols() const { return expr_free_symbols(t_); }

Op_ptr PauliExpBox::dagger() const {
  return std::make_shared<PauliExpBox>(paulis_, -t_);
}

Op_ptr PauliExpBox::transpose() const {
  std::vector<Pauli> paulis = get_paulis();
  int y_pauli_counter = std::count(paulis.begin(), paulis.end(), Pauli::Y);

  // Negate the parameter if odd
  if (y_pauli_counter % 2 == 0) {
    return std::make_shared<PauliExpBox>(paulis_, t_);
  } else {
    return std::make_shared<PauliExpBox>(paulis_, -t_);
  };
}

Op_ptr PauliExpBox::symbol_substitution(
    const SymEngine::map_basic_basic &sub_map) const {
  return std::make_shared<PauliExpBox>(this->paulis_, this->t_.subs(sub_map));
}

void PauliExpBox::generate_circuit() const {
  Circuit circ = pauli_gadget(paulis_, t_);
  circ_ = std::make_shared<Circuit>(circ);
}

composite_def_ptr_t CompositeGateDef::define_gate(
    const std::string &name, const Circuit &def, const std::vector<Sym> &args) {
  return std::make_shared<CompositeGateDef>(name, def, args);
}

CompositeGateDef::CompositeGateDef(
    const std::string &name, const Circuit &def, const std::vector<Sym> &args)
    : name_(name), def_(std::make_shared<Circuit>(def)), args_(args) {}

Circuit CompositeGateDef::instance(const std::vector<Expr> &params) const {
  Circuit circ = *def_;
  symbol_map_t symbol_map;
  for (unsigned i = 0; i < params.size(); i++) {
    symbol_map.insert({args_.at(i), params.at(i)});
  }
  circ.symbol_substitution(symbol_map);
  return circ;
}

op_signature_t CompositeGateDef::signature() const {
  op_signature_t qubs(def_->n_qubits(), EdgeType::Quantum);
  op_signature_t bs(def_->n_bits(), EdgeType::Classical);
  qubs.insert(qubs.end(), bs.begin(), bs.end());
  return qubs;
}

bool CompositeGateDef::operator==(const CompositeGateDef &other) const {
  if (this->get_name() != other.get_name()) return false;
  std::vector<Expr> this_args = {this->args_.begin(), this->args_.end()};
  std::vector<Expr> other_args = {other.args_.begin(), other.args_.end()};
  if (this_args != other_args) return false;
  return this->get_def()->circuit_equality(*other.get_def(), {}, false);
}

CustomGate::CustomGate(
    const composite_def_ptr_t &gate, const std::vector<Expr> &params)
    : Box(OpType::CustomGate), gate_(gate), params_(params) {
  if (!gate) {
    throw std::runtime_error(
        "Null CompositeGateDef pointer passed to CustomGate");
  }
  signature_ = gate->signature();

  if (params_.size() != gate_->n_args()) throw InvalidParameterCount();
}

CustomGate::CustomGate(const CustomGate &other)
    : Box(other), gate_(other.gate_), params_(other.params_) {}

bool CustomGate::is_equal(const Op &op_other) const {
  const CustomGate &other = dynamic_cast<const CustomGate &>(op_other);
  if (this->id_ == other.id_) {
    return true;
  }
  TKET_ASSERT(gate_ && other.gate_);
  return params_ == other.params_ && *gate_ == *other.gate_;
}

Op_ptr CustomGate::symbol_substitution(
    const SymEngine::map_basic_basic &sub_map) const {
  std::vector<Expr> new_params;
  for (const Expr &p : this->params_) {
    new_params.push_back(p.subs(sub_map));
  }
  return std::make_shared<CustomGate>(this->gate_, new_params);
}

void CustomGate::generate_circuit() const {
  circ_ = std::make_shared<Circuit>(gate_->instance(params_));
}

SymSet CustomGate::free_symbols() const { return to_circuit()->free_symbols(); }

std::string CustomGate::get_name(bool) const {
  std::stringstream s;
  s << gate_->get_name();
  if (!params_.empty()) {
    s << "(";
    std::string sep = "";
    for (const Expr &e : params_) {
      s << sep << e;
      sep = ",";
    }
    s << ")";
  }
  return s.str();
}

bool CustomGate::is_clifford() const {
  std::shared_ptr<Circuit> circ = to_circuit();
  BGL_FORALL_VERTICES(v, circ->dag, DAG) {
    if (!circ->get_Op_ptr_from_Vertex(v)->is_clifford()) return false;
  }
  return true;
}

QControlBox::QControlBox(const Op_ptr &op, unsigned n_controls)
    : Box(OpType::QControlBox), op_(op), n_controls_(n_controls) {
  op_signature_t inner_sig = op_->get_signature();
  n_inner_qubits_ = inner_sig.size();
  if (std::count(inner_sig.begin(), inner_sig.end(), EdgeType::Quantum) !=
      n_inner_qubits_) {
    throw NotImplemented("Quantum control of classical wires not supported");
  }
  signature_ = op_signature_t(n_controls + n_inner_qubits_, EdgeType::Quantum);
}

QControlBox::QControlBox(const QControlBox &other)
    : Box(other),
      op_(other.op_),
      n_controls_(other.n_controls_),
      n_inner_qubits_(other.n_inner_qubits_) {}

Op_ptr QControlBox::symbol_substitution(
    const SymEngine::map_basic_basic &sub_map) const {
  return std::make_shared<QControlBox>(
      op_->symbol_substitution(sub_map), n_controls_);
}

SymSet QControlBox::free_symbols() const { return op_->free_symbols(); }

std::string QControlBox::get_command_str(const unit_vector_t &args) const {
  std::stringstream out;
  out << "qif (";
  if (n_controls_ > 0) {
    out << args.at(0).repr();
    for (unsigned i = 1; i < n_controls_; ++i) {
      out << ", " << args.at(i).repr();
    }
  }
  unit_vector_t inner_args(args.begin() + n_controls_, args.end());
  out << ") " << op_->get_command_str(inner_args);
  return out.str();
}

void QControlBox::generate_circuit() const {
  Circuit c(n_inner_qubits_);
  std::vector<unsigned> qbs(n_inner_qubits_);
  std::iota(qbs.begin(), qbs.end(), 0);
  c.add_op(op_, qbs);
  c.decompose_boxes_recursively();
  c = with_controls(c, n_controls_);
  circ_ = std::make_shared<Circuit>(c);
}

Op_ptr QControlBox::dagger() const {
  const Op_ptr inner_dagger = op_->dagger();
  return std::make_shared<QControlBox>(inner_dagger, n_controls_);
}

Op_ptr QControlBox::transpose() const {
  const Op_ptr inner_transpose = op_->transpose();
  return std::make_shared<QControlBox>(inner_transpose, n_controls_);
}

ProjectorAssertionBox::ProjectorAssertionBox(
    const Eigen::MatrixXcd &m, BasisOrder basis)
    : Box(OpType::ProjectorAssertionBox),
      m_(basis == BasisOrder::ilo ? m : reverse_indexing(m)),
      expected_readouts_({}) {
  if ((m.rows() != 2 && m.rows() != 4 && m.rows() != 8) || !is_projector(m)) {
    throw CircuitInvalidity(
        "Matrix for ProjectorAssertionBox must be a 2x2, 4x4, or 8x8 "
        "projector");
  }
  generate_circuit();
}

ProjectorAssertionBox::ProjectorAssertionBox(const ProjectorAssertionBox &other)
    : Box(other), m_(other.m_), expected_readouts_(other.expected_readouts_) {}

Op_ptr ProjectorAssertionBox::dagger() const {
  return std::make_shared<ProjectorAssertionBox>(m_.adjoint());
}

Op_ptr ProjectorAssertionBox::transpose() const {
  return std::make_shared<ProjectorAssertionBox>(m_.transpose());
}

op_signature_t ProjectorAssertionBox::get_signature() const {
  auto circ_ptr = to_circuit();
  op_signature_t qubs(circ_ptr->n_qubits(), EdgeType::Quantum);
  op_signature_t bs(circ_ptr->n_bits(), EdgeType::Classical);
  qubs.insert(qubs.end(), bs.begin(), bs.end());
  return qubs;
}

void ProjectorAssertionBox::generate_circuit() const {
  Circuit c;
  std::tie(c, expected_readouts_) = projector_assertion_synthesis(m_);
  c.decompose_boxes_recursively();
  circ_ = std::make_shared<Circuit>(c);
}

StabiliserAssertionBox::StabiliserAssertionBox(
    const PauliStabiliserList &paulis)
    : Box(OpType::StabiliserAssertionBox),
      paulis_(paulis),
      expected_readouts_({}) {
  generate_circuit();
}

StabiliserAssertionBox::StabiliserAssertionBox(
    const StabiliserAssertionBox &other)
    : Box(other),
      paulis_(other.paulis_),
      expected_readouts_(other.expected_readouts_) {}

Op_ptr StabiliserAssertionBox::dagger() const {
  return std::make_shared<StabiliserAssertionBox>(paulis_);
}

Op_ptr StabiliserAssertionBox::transpose() const {
  PauliStabiliserList new_pauli_list;
  for (auto &pauli : paulis_) {
    int y_pauli_counter =
        std::count(pauli.string.begin(), pauli.string.end(), Pauli::Y);
    if (y_pauli_counter % 2 == 0) {
      new_pauli_list.push_back(PauliStabiliser(pauli.string, pauli.coeff));
    } else {
      new_pauli_list.push_back(PauliStabiliser(pauli.string, !pauli.coeff));
    };
  }
  return std::make_shared<StabiliserAssertionBox>(new_pauli_list);
}

void StabiliserAssertionBox::generate_circuit() const {
  Circuit c;
  std::tie(c, expected_readouts_) = stabiliser_assertion_synthesis(paulis_);
  c.decompose_boxes_recursively();
  circ_ = std::make_shared<Circuit>(c);
}

op_signature_t StabiliserAssertionBox::get_signature() const {
  auto circ_ptr = to_circuit();
  op_signature_t qubs(circ_ptr->n_qubits(), EdgeType::Quantum);
  op_signature_t bs(circ_ptr->n_bits(), EdgeType::Classical);
  qubs.insert(qubs.end(), bs.begin(), bs.end());
  return qubs;
}

nlohmann::json core_box_json(const Box &box) {
  nlohmann::json j;
  j["type"] = box.get_type();
  j["id"] = boost::lexical_cast<std::string>(box.get_id());
  return j;
}

nlohmann::json CircBox::to_json(const Op_ptr &op) {
  const auto &box = static_cast<const CircBox &>(*op);
  nlohmann::json j = core_box_json(box);
  j["circuit"] = *(box.to_circuit());
  return j;
}

Op_ptr CircBox::from_json(const nlohmann::json &j) {
  CircBox box = CircBox(j.at("circuit").get<Circuit>());
  return set_box_id(
      box,
      boost::lexical_cast<boost::uuids::uuid>(j.at("id").get<std::string>()));
}

nlohmann::json Unitary1qBox::to_json(const Op_ptr &op) {
  const auto &box = static_cast<const Unitary1qBox &>(*op);
  nlohmann::json j = core_box_json(box);
  j["matrix"] = box.get_matrix();
  return j;
}

Op_ptr Unitary1qBox::from_json(const nlohmann::json &j) {
  Unitary1qBox box = Unitary1qBox(j.at("matrix").get<Eigen::Matrix2cd>());
  return set_box_id(
      box,
      boost::lexical_cast<boost::uuids::uuid>(j.at("id").get<std::string>()));
}

nlohmann::json Unitary2qBox::to_json(const Op_ptr &op) {
  const auto &box = static_cast<const Unitary2qBox &>(*op);
  nlohmann::json j = core_box_json(box);
  j["matrix"] = box.get_matrix();
  return j;
}

Op_ptr Unitary2qBox::from_json(const nlohmann::json &j) {
  Unitary2qBox box = Unitary2qBox(j.at("matrix").get<Eigen::Matrix4cd>());
  return set_box_id(
      box,
      boost::lexical_cast<boost::uuids::uuid>(j.at("id").get<std::string>()));
}
nlohmann::json Unitary3qBox::to_json(const Op_ptr &op) {
  const auto &box = static_cast<const Unitary3qBox &>(*op);
  nlohmann::json j = core_box_json(box);
  j["matrix"] = box.get_matrix();
  return j;
}

Op_ptr Unitary3qBox::from_json(const nlohmann::json &j) {
  Unitary3qBox box = Unitary3qBox(j.at("matrix").get<Matrix8cd>());
  return set_box_id(
      box,
      boost::lexical_cast<boost::uuids::uuid>(j.at("id").get<std::string>()));
}

nlohmann::json ExpBox::to_json(const Op_ptr &op) {
  const auto &box = static_cast<const ExpBox &>(*op);
  nlohmann::json j = core_box_json(box);
  const auto &matrix_phase = box.get_matrix_and_phase();
  j["matrix"] = matrix_phase.first;
  j["phase"] = matrix_phase.second;
  return j;
}

Op_ptr ExpBox::from_json(const nlohmann::json &j) {
  ExpBox box = ExpBox(
      j.at("matrix").get<Eigen::Matrix4cd>(), j.at("phase").get<double>());
  return set_box_id(
      box,
      boost::lexical_cast<boost::uuids::uuid>(j.at("id").get<std::string>()));
}

nlohmann::json PauliExpBox::to_json(const Op_ptr &op) {
  const auto &box = static_cast<const PauliExpBox &>(*op);
  nlohmann::json j = core_box_json(box);
  j["paulis"] = box.get_paulis();
  j["phase"] = box.get_phase();
  return j;
}

Op_ptr PauliExpBox::from_json(const nlohmann::json &j) {
  PauliExpBox box = PauliExpBox(
      j.at("paulis").get<std::vector<Pauli>>(), j.at("phase").get<Expr>());
  return set_box_id(
      box,
      boost::lexical_cast<boost::uuids::uuid>(j.at("id").get<std::string>()));
}

void to_json(nlohmann::json &j, const composite_def_ptr_t &cdef) {
  j["name"] = cdef->get_name();
  j["definition"] = *cdef->get_def();
  j["args"] = cdef->get_args();
}

void from_json(const nlohmann::json &j, composite_def_ptr_t &cdef) {
  cdef = CompositeGateDef::define_gate(
      j.at("name").get<std::string>(), j.at("definition").get<Circuit>(),
      j.at("args").get<std::vector<Sym>>());
}

nlohmann::json CustomGate::to_json(const Op_ptr &op) {
  const auto &box = static_cast<const CustomGate &>(*op);
  nlohmann::json j = core_box_json(box);
  j["gate"] = box.get_gate();
  j["params"] = box.get_params();
  return j;
}

Op_ptr CustomGate::from_json(const nlohmann::json &j) {
  CustomGate box = CustomGate(
      j.at("gate").get<composite_def_ptr_t>(),
      j.at("params").get<std::vector<Expr>>());
  return set_box_id(
      box,
      boost::lexical_cast<boost::uuids::uuid>(j.at("id").get<std::string>()));
}

nlohmann::json QControlBox::to_json(const Op_ptr &op) {
  const auto &box = static_cast<const QControlBox &>(*op);
  nlohmann::json j = core_box_json(box);
  j["n_controls"] = box.get_n_controls();
  j["op"] = box.get_op();
  return j;
}

Op_ptr QControlBox::from_json(const nlohmann::json &j) {
  QControlBox box =
      QControlBox(j.at("op").get<Op_ptr>(), j.at("n_controls").get<unsigned>());
  return set_box_id(
      box,
      boost::lexical_cast<boost::uuids::uuid>(j.at("id").get<std::string>()));
}

nlohmann::json ProjectorAssertionBox::to_json(const Op_ptr &op) {
  const auto &box = static_cast<const ProjectorAssertionBox &>(*op);
  nlohmann::json j = core_box_json(box);
  j["matrix"] = box.get_matrix();
  return j;
}

Op_ptr ProjectorAssertionBox::from_json(const nlohmann::json &j) {
  ProjectorAssertionBox box =
      ProjectorAssertionBox(j.at("matrix").get<Eigen::MatrixXcd>());
  return set_box_id(
      box,
      boost::lexical_cast<boost::uuids::uuid>(j.at("id").get<std::string>()));
}

nlohmann::json StabiliserAssertionBox::to_json(const Op_ptr &op) {
  const auto &box = static_cast<const StabiliserAssertionBox &>(*op);
  nlohmann::json j = core_box_json(box);
  j["stabilisers"] = box.get_stabilisers();
  return j;
}

Op_ptr StabiliserAssertionBox::from_json(const nlohmann::json &j) {
  StabiliserAssertionBox box =
      StabiliserAssertionBox(j.at("stabilisers").get<PauliStabiliserList>());
  return set_box_id(
      box,
      boost::lexical_cast<boost::uuids::uuid>(j.at("id").get<std::string>()));
}

// use macro to register converters defined in this file with OpJsonFactory
REGISTER_OPFACTORY(CircBox, CircBox)
REGISTER_OPFACTORY(Unitary1qBox, Unitary1qBox)
REGISTER_OPFACTORY(Unitary2qBox, Unitary2qBox)
REGISTER_OPFACTORY(Unitary3qBox, Unitary3qBox)
REGISTER_OPFACTORY(ExpBox, ExpBox)
REGISTER_OPFACTORY(PauliExpBox, PauliExpBox)
REGISTER_OPFACTORY(CustomGate, CustomGate)
REGISTER_OPFACTORY(QControlBox, QControlBox)
REGISTER_OPFACTORY(ProjectorAssertionBox, ProjectorAssertionBox)
REGISTER_OPFACTORY(StabiliserAssertionBox, StabiliserAssertionBox)
}  // namespace tket
