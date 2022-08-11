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
#include <tkassert/Assert.hpp>

#include "CircUtils.hpp"
#include "Circuit/AssertionSynthesis.hpp"
#include "Command.hpp"
#include "Gate/Rotation.hpp"
#include "OpType/OpTypeInfo.hpp"
#include "Ops/OpJsonFactory.hpp"
#include "Ops/OpPtr.hpp"
#include "ThreeQubitConversion.hpp"
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
  circ_ = std::make_shared<Circuit>(three_qubit_tk_synthesis(m_));
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
    throw BadOpType(
        "Quantum control of classical wires not supported", op_->get_type());
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

unsigned get_hamming_distance(const vector<bool> &a, const vector<bool> &b) {
  std::vector<bool> diff = a ^ b;
  return std::accumulate(diff.begin(), diff.end(), 0);
}

std::vector<ToffoliBox::transposition_t> ToffoliBox::cycle_to_transposition(
    const cycle_t &cycle) const {
  // use summation of the hamming distance to decide between two basic options
  // std::cout << "Cycle: " << std::endl;
  // for (auto x : cycle) {
  //   for (auto b : x) {
  //     std::cout << b;
  //   }
  //   std::cout << std::endl;
  // }

  // iterate through cycle_t object pairwise and produce two cheap
  // transpositions options
  cycle_t::const_iterator it = cycle.begin();
  std::vector<bool> first = *it;
  ++it;
  cycle_t::const_iterator jt = cycle.begin();
  // also cost transposition options as they are produced
  unsigned accumulated_hamming_distance0, accumulated_hamming_distance1;
  std::vector<ToffoliBox::transposition_t> transpositions0, transpositions1;
  while (it != cycle.end()) {
    transpositions0.push_back({first, *it});
    transpositions1.push_back({*jt, *it});

    accumulated_hamming_distance0 += get_hamming_distance(first, *it);
    accumulated_hamming_distance1 += get_hamming_distance(*jt, *it);

    ++it;
    ++jt;
  }
  if (accumulated_hamming_distance0 > accumulated_hamming_distance1) {
    return transpositions1;
  }
  return transpositions0;
}

std::set<std::vector<ToffoliBox::transposition_t>>
ToffoliBox::get_transpositions() const {
  std::set<std::vector<ToffoliBox::transposition_t>> transpositions;
  for (const cycle_t &cycle : this->cycles_) {
    transpositions.insert(this->cycle_to_transposition);
  }
  return transpositions;
}

// TODO: check every vector is the same size
ToffoliBox::ToffoliBox(
    std::map<std::vector<bool>, std::vector<bool>> &_permutation)
    : Box(OpType::ToffoliBox) {
  // Convert passed permutation to cycles
  while (!_permutation.empty()) {
    auto it = _permutation.begin();
    // std::cout << "{";
    // for(auto b : it->first){
    //   std::cout << b << " ";
    // }
    // std::cout << "}" << std::endl;
    cycle_t cycle = {it->first};
    it = _permutation.find(it->second);
    while (it->first != cycle[0]) {
      // std::cout << "While Loop" << std::endl;
      // std::cout << "{";
      // for(auto b : it->first){
      //   std::cout << b << " ";
      // }
      // std::cout << "}" << std::endl;
      cycle.push_back(it->first);
      it = _permutation.find(it->second);
      if (it == _permutation.end()) {
        throw std::invalid_argument("Permutation is not complete.");
      }
    }
    // std::cout << "Produced Cycle!" << std::endl;
    // for(auto x : cycle){
    //   std::cout << "{";
    //   for(auto b : x){
    //     std::cout << b << " ";
    //   }
    //   std::cout << "}" << std::endl;
    // }
    // std::cout << cycle.size() << std::endl;
    if (cycle.size() > 1) {
      this->cycles_.push_back(cycle);
    }
    // TODO: what's quicker, deleting entries or adding to set for look up?
    for (const std::vector<bool> &bitstring : cycle) {
      _permutation.erase(bitstring);
    }
  }
}

Circuit ToffoliBox::get_bitstring_circuit(
    const std::vector<bool> &bitstring, const unsigned &target) const {
  // flip qubits that need to be state 0
  unsigned n_qubits = bitstring.size();
  Circuit x_circuit(n_qubits);
  std::vector<unsigned> controls;
  for (unsigned i = 0; i < n_qubits; i++) {
    if (i != target) {
      if (!bitstring[i]) {
        x_circuit.add_op<unsigned>(OpType::X, {i});
      }
      controls.push_back(i);
    }
  }
  controls.push_back(target);
  TKET_ASSERT(controls.size() == n_qubits);

  Circuit return_circuit(n_qubits);
  return_circuit.append(x_circuit);
  return_circuit.add_op<unsigned>(OpType::CnX, controls);
  return_circuit.append(x_circuit);
  return return_circuit;
}

gray_code_t ToffoliBox::order_cycle_graycodes(
    const std::set<std::vector<
        std::tuple<std::vector<bool>, std::vector<bool>, gray_code_t>>>
        &transposition_gray_codes) const {
  // First note that the transpositions have been produced to attempt to reduce
  // hamming distance in cycle Gray code have then been found for each
  // transposition Two rules:
  // 1) for two sequential bitstrings b0 and b1 with
  // target indices i0 and i1 can optionally do b0[i0] = !b0[i0] and b0[i1] =
  // !b0[i1]
  // 2) cycles can be ordered s.t. the last transposition_t in a vector
  // matches the first in the next, meaning they can be cancelled out

  // Reduction is done with these two rules in mind by first looking at all
  // graycodes for a single cycle and attempt to switch graycode entries to
  // reduce overlap at merge points, and then ordering cycles s.t. graycode
  // entries can cancel (with swapping if needed)

  // first order cycles so as to encourage cancellation
  std::vector<std::tuple<std::vector<bool>, std::vector<bool>, gray_code_t>>
      sequenced_graycodes;
  for (const std::vector<
           std::tuple<std::vector<bool>, std::vector<bool>, gray_code_t>>
           &cycle : transposition_gray_codes) {
    // a cycle can start with one of two bitstrings + replacement target and end
    // with one of two bitstrings + replacement target
    std::set<std::pair<std::vector<bool, unsigned>>> initial_options;
    // first find these
    auto it = cycle.begin();
    gray_code_t gray_code_front = std::get<2>(*it);
    TKET_ASSERT(!gray_code_front.empty());
    // if transpositions had hamming distance 1 then only one option
    initial_options.insert(gray_code_front[0]);
    if (gray_code_front.size() > 1) {
      std::pair<std::vector<bool>, unsigned> first = gray_code_front[0];
      unsigned second_target = gray_code_front[1].second;
      first[second_target] = !first[second_target];
      first[first.second] = !first[first.second];
      initial_options.insert(first);
    }

    // repeat for final options
    std::set<std::pair<std::vector<bool, unsigned>>> final_options;
    it = cycle.rbegin();
    gray_code_t gray_code_back = std::get<2>(*it);
    TKET_ASSERT(!gray_code_back.empty());
    final_options.insert(gray_code_back[0]);
    if (gray_code_back.size() > 1) {
      std::pair<std::vector<bool>, unsigned> first = gray_code_back[0];
      unsigned second_target = gray_code_back[1].second;
      first[second_target] = !first[second_target];
      first[first.second] = !first[first.second];
      final_options.insert(first);
    }
  }
}

std::set<std::pair<std::vector<bool, unsigned>>> gen_gray_code_step(
    const std::vector<bool> &a, const std::vector<bool> &b) {
  if (a.size() != b.size()) {
    throw std::invalid_argument("Bitstrings must have identical length.");
  }
  std::set < std::pair<std::vector<bool, unsigned>> output;
  for (unsigned i = 0; i < a.size(); i++) {
    if (a[i] != b[i]) {
      std::vector<bool> entry = a;
      entry[i] = !entry[i];
      output.insert({entry, i});
    }
  }
  return output;
}

// attempting to merge b to back of a
std::pair<std::vector<bool>, std::vector<bool>> merge_transpositions(
    const std::vector<bool> &a_last, const std::vector<bool> &a_middle,
    const std::vector<bool> &b_first, const std::vector<bool> &b_middle) {
  //
  std::set<std::pair<std::vector<bool, unsigned>>> transposition_a_graycodes =
      gen_gray_code_step(a_last, b_middle);
  std::set<std::pair<std::vector<bool, unsigned>>> transposition_b_graycodes =
      gen_gray_code_step(b_first, a_middle);

  std::set<std::pair<std::vector<bool, unsigned>>> intersection;
  std::set_intersection(
      transposition_a_graycodes.begin(), transposition_a_graycodes.end(),
      transposition_b_graycodes.begin(), transposition_b_graycodes.end(),
      std::back_inserter(intersection));

  if (!intersections.empty()) {
    while (true) {
      std::set<std::pair<std::vector<bool, unsigned>>> new_intersections = {};
      for (const std::pair<std::vector<bool, unsigned>> &step : intersection) {
        // number of returned steps is empty if the target has already been
        // reached all instances of transposition_a_graycodes or
        // transposition_b_graycodes will be empty if one is
        transposition_a_graycodes = gen_gray_code_step(step.first, a_middle);
        transposition_b_graycodes = gen_gray_code_step(step.first, b_middle);
        std::set_intersection(
            transposition_a_graycodes.begin(), transposition_a_graycodes.end(),
            transposition_b_graycodes.begin(), transposition_b_graycodes.end(),
            std::back_inserter(new_intersections));
      }
      if (new_intersections.empty()) {
        break;
      }
      intersections = new_intersections;
    }
    // if matches, just get first in intersection
    std::pair<std::vector<bool, unsigned>> match = *intersection.begin();
    return {match, match};
  } else {
    return {a_last, b_first};
  }
}

cycle_transposition_t ToffoliBox::merge_cycles(
    std::set<cycle_transposition_t> &cycle_transpositions) const {
  // For comparison, use bool from construction to choose how cycle is made

  // First note that the transpositions have been produced to attempt to reduce
  // hamming distance in cycle.
  // There are multiple gray code possible between transpositions.
  // A transposition is defined as a {first, middle, last} struct

  // Taking each transposition in the a cycle pairwise, produce an option of
  // first step of gray code based on a last->middle of the first entry in the
  // pair, and first->middle of the second entry in the pair for all that match,
  // produce a new option of match->middle of the first and match->middle of the
  // second repeat until there is one match left (or take one of the "equal"
  // matches) set the first transposition to {first, middle, match}, set the
  // second transposition to {match, middle, last}
  for (const cycle_transposition_t &cycle : cycle_transpositions) {
    auto it = cycle.begin();
    auto jt = it;
    ++jt;
    while (jt != cycle.end()) {
      std::pair<std::vector<bool>, std::vector<bool>> update =
          merge_transpositions(it->last, it->middle, jt->first, jt->middle);
      it->last = update.first;
      jt->first = update.second;
      ++it;
      ++jt;
    }
  }

  // Next repeat this process with the whole cycles.
  // A cycle is defined by the {first, middle} of its first transposition and
  // then {middle, last} of its final transposition Go through permutations (can
  // be capped) of ordered pairs pt0, pt1 of cycles. Generate range of first
  // step of gray code for pt0 last->middle and pt1 first->middle If any match,
  // produce options as before, repeat until best chosen set pt0 last
  // transposition last to match, set pt1 first transposition to match combine
  // into new cycle, repeat from top

  auto it = cycle_transpositions.begin();
  while (it != cycle_transpositions.end()) {
    // it and jt are cycle_transposition_t objects
    auto jt = it;

    while (jt != cycle_transpositions.end()) {
      // try adding cycle it to back of cycle jt
      // TODO: instead of trying one then the other, try both and pick ones with
      // best hamming distance to indivdual middles
      std::pair<std::vector<bool>, std::vector<bool>> output =
          merge_transpositions(jt->last, jt->middle, it->first, it->middle);
      if (output.first != jt->last) {
        jt->last = pair.first;
        it->first = pair.second;
        // combine cycle transpositions
        jt->insert(jt->end(), it->begin(), it->end());
        // n.b. jt after it, so it now == jt != cycle_transpositions.end()
        // but we increment jt immediately after
        it = cycle_transpositions.erase(it);
      } else {
        // try adding cycle jt to back of cycle it
        output =
            merge_transpositions(it->last, it->middle, jt->first, jt->middle);
        if (output.first != it->last) {
          it->last = pair.first;
          jt->first = pair.second;
          // combine cycle transpositions
          // such new information is held in jt
          jt->insert(jt->begin(), it->rbegin(), it->rend());
          // n.b. jt after it, so now it == jt != cycle_transpositions.end()
          // but we increment jt immediately after
          it = cycle_transpositions.erase(it);
        }
      }
      ++jt;
    }
    ++it;
  }

  // std::set<std::pair<std::vector<bool, unsigned>>>
  // transposition_a_graycodes = gen_gray_code_step(it->first, it->middle);
  // std::set<std::pair<std::vector<bool, unsigned>>>
  //     transposition_b_graycodes = gen_gray_code_step(jt->last, jt->middle);

  // std::set<std::pair<std::vector<bool, unsigned>>> intersection;
  // std::set_intersection(
  //     transposition_a_graycodes.begin(), transposition_a_graycodes.end(),
  //     transposition_b_graycodes.begin(), transposition_b_graycodes.end(),
  //     std::back_inserter(intersection));

  // if (!intersections.empty()) {
  //   // find best gray code match
  //   while (true) {
  //     std::set<std::pair<std::vector<bool, unsigned>>> new_intersections =
  //         {};
  //     for (const std::pair<std::vector<bool, unsigned>> &step :
  //          intersection) {
  //       // number of returned steps is empty if the target has already been
  //       // reached all instances of transposition_a_graycodes or
  //       // transposition_b_graycodes will be empty if one is
  //       transposition_a_graycodes =
  //           gen_gray_code_step(step.first, it->middle);
  //       transposition_b_graycodes =
  //           gen_gray_code_step(step.first, jt->middle);
  //       std::set_intersection(
  //           transposition_a_graycodes.begin(),
  //           transposition_a_graycodes.end(),
  //           transposition_b_graycodes.begin(),
  //           transposition_b_graycodes.end(),
  //           std::back_inserter(new_intersections));
  //     }
  //     if (new_intersections.empty()) {
  //       break;
  //     }
  //     intersections = new_intersections;
  //   }
  //   // if matches, just get first in intersection
  //   std::pair<std::vector<bool, unsigned>> match =
  //       *intersection.begin() it->first = match;
  //   jt->last = match;
  //   // combine cycle transpositions
  //   jt->insert(jt->end(), it->begin(), it->end());
  //   // n.b. jt after it, so it now == jt != cycle_transpositions.end()
  //   // but we increment jt immediately after
  //   it = cycle_transpositions.erase(it);
  //     } else {
  //       // try with back
  //       transposition_a_graycodes = gen_gray_code_step(it->last, it->middle);
  //       transposition_b_graycodes = gen_gray_code_step(jt->first,
  //       jt->middle); intersection.clear(); std::set_intersection(
  //           transposition_a_graycodes.begin(),
  //           transposition_a_graycodes.end(),
  //           transposition_b_graycodes.begin(),
  //           transposition_b_graycodes.end(),
  //           std::back_inserter(intersection));
  //       if (!intersections.empty()) {
  //         while (true) {
  //           std::set <
  //               std::pair<std::vector<bool, unsigned>> new_intersections =
  //               {};
  //           for (const std::pair < std::vector < bool,
  //                unsigned >>> &step : intersection) {
  //             // number of returned steps is empty if the target has already
  //             // been reached all instances of transposition_a_graycodes or
  //             // transposition_b_graycodes will be empty if one is
  //             transposition_a_graycodes =
  //                 gen_gray_code_step(step.first, it->middle);
  //             transposition_b_graycodes =
  //                 gen_gray_code_step(step.first, jt->middle);
  //             std::set_intersection(
  //                 transposition_a_graycodes.begin(),
  //                 transposition_a_graycodes.end(),
  //                 transposition_b_graycodes.begin(),
  //                 transposition_b_graycodes.end(),
  //                 std::back_inserter(new_intersections));
  //           }
  //           if (new_intersections.empty()) {
  //             break;
  //           }
  //           intersections = new_intersections;
  //         }

  //         // if matches, just get first in intersection
  //         std::pair<std::vector<bool, unsigned>> match =
  //             *intersection.begin() it->last = match;
  //         jt->first = match;
  //         // combine cycle transpositions
  //         // such new information is held in jt
  //         jt->insert(jt->begin(), it->rbegin(), it->rend());
  //         // n.b. jt after it, so now it == jt != cycle_transpositions.end()
  //         // but we increment jt immediately after
  //         it = cycle_transpositions.erase(it);
  //       }
  //     }
  //     ++jt;
  //   }
  //   ++it;
  // }
}

void ToffoliBox::generate_circuit() const {
  // This decomposition is as described on page 191, section 4.5.2 "Single qubit
  // and CNOT gates are universal" of Nielsen & Chuang
  std::set<std::vector<ToffoliBox::transposition_t>> cycle_transpositions =
      this->get_transpositions();

  std::vector<ToffoliBox::transposition_t> ordered_transpositions =
      this->reorder_transpositions(cycle_transpositions);

  std::cout << "Generate circuit, post get_transpositions:" << std::endl;
  std::cout << "n of sets of transpositions: " << transpositions.size()
            << std::endl;
  if (transpositions.empty()) {
    this->circ_ = std::make_shared<Circuit>();
    return;
  }
  std::set<std::vector<gray_code_t>> all_cycle_gray_codes;
  for (const std::vector<ToffoliBox::transposition_t>
           &single_cycle_transpositions : transpositions) {
    TKET_ASSERT(single_cycle_transpositions.size() > 0);
    size_t n_qubits = single_cycle_transpositions.begin()->first.size();
    std::vector<gray_code_t> single_cycle_gray_codes;
    // each transposition is converted into a gray code
    for (const ToffoliBox::transposition_t &transposition :
         single_cycle_transpositions) {
      TKET_ASSERT(transposition.first.size() == n_qubits);
      TKET_ASSERT(transposition.second.size() == n_qubits);
      // find a sequence of bitstrings one flip away from eachother
      gray_code_t bitstrings;
      std::vector<bool> last_bitstring = transposition.first;
      for (unsigned i = 0; i < transposition.first.size(); i++) {
        // if different, update last bitstring with entry
        if (transposition.first[i] != transposition.second[i]) {
          last_bitstring[i] = !last_bitstring[i];
          bitstrings.push_back({last_bitstring, i});
        }
      }

      gray_code_t gray_code = {{bitstrings[0]}, bitstrings};
      // => first bitstring has alternative, which could lead to later
      // optimisations
      if (bitstrings.size() > 1) {
        std::pair<std::vector<bool>, unsigned> first = bitstrings[0];
        unsigned second_target = bitstrings[1].second;
        first[second_target] = !first[second_target];
        first[first.second] = !first[first.second];
        gray_code.first.insert(first);
      }

      single_cycle_gray_codes.push_back(gray_code);
    }
    all_cycle_gray_codes.push_back(single_cycle_gray_codes);
  }

  gray_code_t concatenated_gray_codes =
      this->order_cycle_graycodes(all_cycle_gray_codes);

  TKET_ASSERT(!concatenated_gray_codes.empty());
  this->circuit_ =
      std::make_shared<Circuit>(concatenated_gray_codes[0].first.size());
  for (const std::pair<std::vector<bool>, unsigned> &bitstring :
       concatenated_gray_codes) {
    this->circ_->append(
        this->get_bitstring_circuit(bitstring.first, bitstring.second));
  }

  // std::vector<std::vector<std::pair<std::vector<bool>, unsigned>>>
  // reordered_bitstrings =  this->reorder_bitstrings(all_bitstrings); for(const
  // std::vector<std::pair<std::vector<bool>, unsigned>>& bitstrings :
  // reordered_bitstrings){
  // }
  // if (!transpositions.empty()) {
  //   size_t n_qubits = transpositions.begin()->first.size();
  //   this->circ_ = std::make_shared<Circuit>(n_qubits);
  //   for (const ToffoliBox::transposition_t &transposition : transpositions) {
  //     TKET_ASSERT(transposition.first.size() == n_qubits);
  //     TKET_ASSERT(transposition.second.size() == n_qubits);
  //     std::cout << "Transposition to convert: " << std::endl;
  //     for(auto b : transposition.first){
  //       std::cout << b;
  //     }
  //     std::cout << " | ";
  //     for(auto b : transposition.second){
  //       std::cout << b;
  //     }
  //     std::cout << std::endl;
  //     // find a sequence of bitstrings one flip away from eachother
  //     std::vector<std::pair<std::vector<bool>, unsigned>> bitstrings;
  //     std::vector<bool> last_bitstring = transposition.first;
  //     for (unsigned i = 0; i < transposition.first.size(); i++) {
  //       // if different, update last bitstring with entry
  //       if (transposition.first[i] != transposition.second[i]) {
  //         last_bitstring[i] = !last_bitstring[i];
  //         bitstrings.push_back({last_bitstring, i});
  //       }
  //     }
  //     std::cout << "Bitstrings produced: " << std::endl;
  //     for(auto bs : bitstrings){
  //       for(auto bo : bs.first){
  //         std::cout << bo;
  //       }
  //       std::cout << " | " << bs.second << std::endl;
  //     }
  //     std::vector<std::pair<std::vector<bool>, unsigned>>::const_iterator it
  //     =
  //         bitstrings.begin();
  //     while (it != bitstrings.end()) {
  //       this->circ_->append(this->get_bitstring_circuit(it->first,
  //       it->second));
  //       ++it;
  //     }
  //     // as one past end & don't want to undo gate
  //     std::vector<
  //         std::pair<std::vector<bool>, unsigned>>::const_reverse_iterator jt
  //         = bitstrings.rbegin();
  //     std::advance(jt, 1);
  //     // uncompute basis flips
  //     while (jt != bitstrings.rend()) {
  //       this->circ_->append(this->get_bitstring_circuit(jt->first,
  //       jt->second));
  //       ++jt;
  //     }
  //   }
  //   std::cout << "Finished making circuit! " << std::endl;
  // } else {
  //   this->circ_ = std::make_shared<Circuit>();
  // }
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
