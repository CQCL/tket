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

#include "PauliGraph.hpp"

namespace tket {
namespace pg {

/**
 * PGOp Abstract Implementation
 */

PGOpType PGOp::get_type() const { return type_; }

bool PGOp::operator==(const PGOp& other) const {
  return (type_ == other.type_) && is_equal(other);
}

PGOp::~PGOp() {}

PGOp::PGOp(PGOpType type) : type_(type) {}

bool PGOp::commutes_with(const PGOp& other) const {
  for (const QubitPauliTensor& t : active_paulis()) {
    for (const QubitPauliTensor& ot : other.active_paulis()) {
      if (!t.commutes_with(ot)) return false;
    }
  }
  for (const Bit& b : write_bits()) {
    for (const Bit& ob : other.write_bits()) {
      if (b == ob) return false;
    }
    for (const Bit& ob : other.read_bits()) {
      if (b == ob) return false;
    }
  }
  for (const Bit& b : read_bits()) {
    for (const Bit& ob : other.write_bits()) {
      if (b == ob) return false;
    }
  }
  return true;
}

unsigned PGOp::n_paulis() const { return 1; }

bit_vector_t PGOp::read_bits() const { return {}; }

bit_vector_t PGOp::write_bits() const { return {}; }

/**
 * PGRotation Implementation
 */

const QubitPauliTensor& PGRotation::get_tensor() const { return tensor_; }

const Expr& PGRotation::get_angle() const { return angle_; }

PGRotation::PGRotation(const QubitPauliTensor& tensor, const Expr& angle)
    : PGOp(PGOpType::Rotation), tensor_(tensor), angle_(angle) {
  if (std::abs(tensor.coeff + 1.) < EPS) {
    angle_ *= -1.;
    tensor_.coeff = 1.;
  } else if (std::abs(tensor.coeff - 1.) >= EPS)
    throw PGError("Invalid coefficient in tensor for PauliGraph rotation");
}

SymSet PGRotation::free_symbols() const { return expr_free_symbols(angle_); }

PGOp_ptr PGRotation::symbol_substitution(
    const SymEngine::map_basic_basic& sub_map) const {
  return std::make_shared<PGRotation>(tensor_, angle_.subs(sub_map));
}

std::string PGRotation::get_name(bool) const {
  std::stringstream str;
  str << "Rot(" << tensor_.to_str() << "; " << angle_ << ")";
  return str.str();
}

bool PGRotation::is_equal(const PGOp& op_other) const {
  const PGRotation& other = dynamic_cast<const PGRotation&>(op_other);
  return (tensor_ == other.tensor_) && equiv_expr(angle_, other.angle_, 2.);
}

std::vector<QubitPauliTensor> PGRotation::active_paulis() const {
  return {tensor_};
}

QubitPauliTensor& PGRotation::port(unsigned p) {
  if (p != 0)
    throw PGError(
        "Cannot dereference port on PGRotation: " + std::to_string(p));
  return tensor_;
}

/**
 * PGCliffordRot Implementation
 */

const QubitPauliTensor& PGCliffordRot::get_tensor() const { return tensor_; }

unsigned PGCliffordRot::get_angle() const { return angle_; }

PGCliffordRot::PGCliffordRot(const QubitPauliTensor& tensor, unsigned angle)
    : PGOp(PGOpType::CliffordRot), tensor_(tensor), angle_(angle) {
  if (tensor_.coeff == -1.) {
    angle_ = (4 - angle) % 4;
    tensor_.coeff = 1.;
  } else if (tensor_.coeff != 1.)
    throw PGError(
        "Invalid coefficient in tensor for PauliGraph Clifford rotation");
}

SymSet PGCliffordRot::free_symbols() const { return {}; }

PGOp_ptr PGCliffordRot::symbol_substitution(
    const SymEngine::map_basic_basic&) const {
  return PGOp_ptr();
}

std::string PGCliffordRot::get_name(bool) const {
  std::stringstream str;
  str << "ClfRot(" << tensor_.to_str() << "; " << (angle_ * 0.5) << ")";
  return str.str();
}

bool PGCliffordRot::is_equal(const PGOp& op_other) const {
  const PGCliffordRot& other = dynamic_cast<const PGCliffordRot&>(op_other);
  return (tensor_ == other.tensor_) && (angle_ % 4 == other.angle_ % 4);
}

std::vector<QubitPauliTensor> PGCliffordRot::active_paulis() const {
  return {tensor_};
}

QubitPauliTensor& PGCliffordRot::port(unsigned p) {
  if (p != 0)
    throw PGError(
        "Cannot dereference port on PGCliffordRot: " + std::to_string(p));
  return tensor_;
}

/**
 * PGMeasure Implementation
 */

const QubitPauliTensor& PGMeasure::get_tensor() const { return tensor_; }

const Bit& PGMeasure::get_target() const { return target_; }

PGMeasure::PGMeasure(const QubitPauliTensor& tensor, const Bit& target)
    : PGOp(PGOpType::Measure), tensor_(tensor), target_(target) {
  if (tensor_.coeff != 1. && tensor_.coeff != -1.)
    throw PGError("Invalid coefficient in tensor for PauliGraph measurement");
}

SymSet PGMeasure::free_symbols() const { return {}; }

PGOp_ptr PGMeasure::symbol_substitution(
    const SymEngine::map_basic_basic&) const {
  return PGOp_ptr();
}

std::string PGMeasure::get_name(bool) const {
  std::stringstream str;
  str << "Meas(" << tensor_.to_str() << " -> " << target_.repr() << ")";
  return str.str();
}

bool PGMeasure::is_equal(const PGOp& op_other) const {
  const PGMeasure& other = dynamic_cast<const PGMeasure&>(op_other);
  return (tensor_ == other.tensor_) && (target_ == other.target_);
}

std::vector<QubitPauliTensor> PGMeasure::active_paulis() const {
  return {tensor_};
}

QubitPauliTensor& PGMeasure::port(unsigned p) {
  if (p != 0)
    throw PGError("Cannot dereference port on PGMeasure: " + std::to_string(p));
  return tensor_;
}

bit_vector_t PGMeasure::write_bits() const { return {target_}; }

/**
 * PGDecoherence Implementation
 */

const QubitPauliTensor& PGDecoherence::get_tensor() const { return tensor_; }

PGDecoherence::PGDecoherence(const QubitPauliTensor& tensor)
    : PGOp(PGOpType::Decoherence), tensor_(tensor) {
  if (tensor_.coeff == -1.)
    tensor_.coeff = 1.;
  else if (tensor_.coeff != 1.)
    throw PGError("Invalid coefficient in tensor for PauliGraph decoherence");
}

SymSet PGDecoherence::free_symbols() const { return {}; }

PGOp_ptr PGDecoherence::symbol_substitution(
    const SymEngine::map_basic_basic&) const {
  return PGOp_ptr();
}

std::string PGDecoherence::get_name(bool) const {
  std::stringstream str;
  str << "Deco(" << tensor_.to_str() << ")";
  return str.str();
}

bool PGDecoherence::is_equal(const PGOp& op_other) const {
  const PGDecoherence& other = dynamic_cast<const PGDecoherence&>(op_other);
  return tensor_ == other.tensor_;
}

std::vector<QubitPauliTensor> PGDecoherence::active_paulis() const {
  return {tensor_};
}

QubitPauliTensor& PGDecoherence::port(unsigned p) {
  if (p != 0)
    throw PGError(
        "Cannot dereference port on PGDecoherence: " + std::to_string(p));
  return tensor_;
}

/**
 * PGReset Implementation
 */

const QubitPauliTensor& PGReset::get_stab() const { return stab_; }

const QubitPauliTensor& PGReset::get_destab() const { return destab_; }

PGReset::PGReset(const QubitPauliTensor& stab, const QubitPauliTensor& destab)
    : PGOp(PGOpType::Reset), stab_(stab), destab_(destab) {
  if (destab_.coeff == -1.)
    destab_.coeff = 1.;
  else if (destab_.coeff != 1. || (stab_.coeff != 1. && stab_.coeff != -1.))
    throw PGError("Invalid coefficient in tensor for PauliGraph reset");
}

SymSet PGReset::free_symbols() const { return {}; }

PGOp_ptr PGReset::symbol_substitution(const SymEngine::map_basic_basic&) const {
  return PGOp_ptr();
}

std::string PGReset::get_name(bool) const {
  std::stringstream str;
  str << "Reset(" << stab_.to_str() << "; " << destab_.to_str() << ")";
  return str.str();
}

bool PGReset::is_equal(const PGOp& op_other) const {
  const PGReset& other = dynamic_cast<const PGReset&>(op_other);
  return (stab_ == other.stab_) && (destab_ == other.destab_);
}

unsigned PGReset::n_paulis() const { return 2; }

std::vector<QubitPauliTensor> PGReset::active_paulis() const {
  return {stab_, destab_};
}

QubitPauliTensor& PGReset::port(unsigned p) {
  if (p == 0)
    return stab_;
  else if (p == 1)
    return destab_;
  else
    throw PGError("Cannot dereference port on PGReset: " + std::to_string(p));
}

/**
 * PGConditional Implementation
 */

PGOp_ptr PGConditional::get_inner_op() const { return inner_; }

const bit_vector_t& PGConditional::get_args() const { return args_; }

unsigned PGConditional::get_value() const { return value_; }

PGConditional::PGConditional(
    PGOp_ptr inner, const bit_vector_t& args, unsigned value)
    : PGOp(PGOpType::Conditional), inner_(inner), args_(args), value_(value) {}

SymSet PGConditional::free_symbols() const { return inner_->free_symbols(); }

PGOp_ptr PGConditional::symbol_substitution(
    const SymEngine::map_basic_basic& sub_map) const {
  PGOp_ptr inner_sub = inner_->symbol_substitution(sub_map);
  if (inner_sub)
    return std::make_shared<PGConditional>(inner_sub, args_, value_);
  else
    return PGOp_ptr();
}

std::string PGConditional::get_name(bool latex) const {
  std::stringstream str;
  str << "[";
  for (const Bit& b : args_) {
    str << b.repr() << ", ";
  }
  str << "\b\b] == " << value_ << " ? " << inner_->get_name(latex);
  return str.str();
}

bool PGConditional::is_equal(const PGOp& op_other) const {
  const PGConditional& other = dynamic_cast<const PGConditional&>(op_other);
  return (value_ == other.value_) && (args_ == other.args_) &&
         (*inner_ == *other.inner_);
}

unsigned PGConditional::n_paulis() const { return inner_->n_paulis(); }

std::vector<QubitPauliTensor> PGConditional::active_paulis() const {
  return inner_->active_paulis();
}

QubitPauliTensor& PGConditional::port(unsigned p) { return inner_->port(p); }

bit_vector_t PGConditional::read_bits() const {
  bit_vector_t bits = inner_->read_bits();
  bits.insert(bits.end(), args_.begin(), args_.end());
  return bits;
}

bit_vector_t PGConditional::write_bits() const { return inner_->write_bits(); }

/**
 * PGStabilizer Implementation
 */

const QubitPauliTensor& PGStabilizer::get_stab() const { return stab_; }

SymSet PGStabilizer::free_symbols() const { return {}; }

PGOp_ptr PGStabilizer::symbol_substitution(
    const SymEngine::map_basic_basic&) const {
  return PGOp_ptr();
}

std::string PGStabilizer::get_name(bool) const {
  std::stringstream str;
  str << "Stab(" << stab_.to_str() << ")";
  return str.str();
}

bool PGStabilizer::is_equal(const PGOp& op_other) const {
  const PGStabilizer& other = dynamic_cast<const PGStabilizer&>(op_other);
  return (stab_ == other.stab_);
}

std::vector<QubitPauliTensor> PGStabilizer::active_paulis() const {
  return {stab_};
}

QubitPauliTensor& PGStabilizer::port(unsigned p) {
  if (p != 0)
    throw PGError(
        "Cannot dereference port on PGStabilizer: " + std::to_string(p));
  return stab_;
}

/**
 * PGInputRow Implementation
 */

const QubitPauliTensor& PGInputRow::get_tensor() const { return tensor_; }

const QubitPauliTensor& PGInputRow::get_input_tensor() const {
  return input_tensor_;
}

SymSet PGInputRow::free_symbols() const { return {}; }

PGOp_ptr PGInputRow::symbol_substitution(
    const SymEngine::map_basic_basic&) const {
  return PGOp_ptr();
}

std::string PGInputRow::get_name(bool) const {
  std::stringstream str;
  str << "Input(" << input_tensor_.to_str() << " -> " << tensor_.to_str()
      << ")";
  return str.str();
}

bool PGInputRow::is_equal(const PGOp& op_other) const {
  const PGInputRow& other = dynamic_cast<const PGInputRow&>(op_other);
  return (input_tensor_ == other.input_tensor_) && (tensor_ == other.tensor_);
}

std::vector<QubitPauliTensor> PGInputRow::active_paulis() const {
  return {tensor_};
}

QubitPauliTensor& PGInputRow::port(unsigned p) {
  if (p != 0)
    throw PGError(
        "Cannot dereference port on PGInputRow: " + std::to_string(p));
  return tensor_;
}

/**
 * PGOutputRow Implementation
 */

const QubitPauliTensor& PGOutputRow::get_tensor() const { return tensor_; }

const QubitPauliTensor& PGOutputRow::get_output_tensor() const {
  return output_tensor_;
}

SymSet PGOutputRow::free_symbols() const { return {}; }

PGOp_ptr PGOutputRow::symbol_substitution(
    const SymEngine::map_basic_basic&) const {
  return PGOp_ptr();
}

std::string PGOutputRow::get_name(bool) const {
  std::stringstream str;
  str << "Output(" << tensor_.to_str() << " -> " << output_tensor_.to_str()
      << ")";
  return str.str();
}

bool PGOutputRow::is_equal(const PGOp& op_other) const {
  const PGOutputRow& other = dynamic_cast<const PGOutputRow&>(op_other);
  return (output_tensor_ == other.output_tensor_) && (tensor_ == other.tensor_);
}

std::vector<QubitPauliTensor> PGOutputRow::active_paulis() const {
  return {tensor_};
}

QubitPauliTensor& PGOutputRow::port(unsigned p) {
  if (p != 0)
    throw PGError(
        "Cannot dereference port on PGOutputRow: " + std::to_string(p));
  return tensor_;
}

/**
 * PauliGraph Implementation
 */

PauliGraph::PauliGraph()
    : pauli_ac_(0, 0),
      pauli_index_(),
      c_graph_(),
      qubits_(),
      bits_(),
      last_writes_(),
      last_reads_() {}

void PauliGraph::to_graphviz(std::ostream& out) const {
  out << "digraph G {\ncompound = true;\n";

  std::map<PGVert, unsigned> v_map;
  unsigned i = 0;
  BGL_FORALL_VERTICES(v, c_graph_, PGClassicalGraph) {
    v_map.insert({v, i});
    PGOp_ptr op = c_graph_[v];
    out << "subgraph v" << i << "{\nlabel = \"" << op->get_name() << "\";\n";
    if (op->n_paulis() == 0) {
      out << i << "_0;\n";
    } else {
      for (unsigned j = 0; j < op->n_paulis(); ++j) {
        out << i << "_" << j << ";\n";
      }
    }
    out << "}\n";
    ++i;
  }

  BGL_FORALL_EDGES(e, c_graph_, PGClassicalGraph) {
    PGVert vs = boost::source(e, c_graph_);
    PGVert vt = boost::target(e, c_graph_);
    unsigned vsi = v_map.at(vs);
    unsigned vti = v_map.at(vt);
    out << vsi << "_0 -> " << vti << "_0 [ltail=v" << vsi << ",lhead=v" << vti
        << "];\n";
  }

  for (const PGPauli& r_pauli : pauli_index_.get<TagID>()) {
    for (const PGPauli& c_pauli : pauli_index_.get<TagID>()) {
      if (pauli_ac_(r_pauli.index, c_pauli.index)) {
        out << v_map.at(c_pauli.vert) << "_" << c_pauli.port << " -> "
            << v_map.at(r_pauli.vert) << "_" << r_pauli.port << ";\n";
      }
    }
  }

  out << "}";
}

PGVert PauliGraph::add_vertex_at_end(PGOp_ptr op) {
  // Add vertex
  PGVert v = boost::add_vertex(c_graph_);
  c_graph_[v] = op;
  unsigned mat_offset = pauli_ac_.rows();
  pauli_ac_.conservativeResize(
      mat_offset + op->n_paulis(), mat_offset + op->n_paulis());
  for (unsigned i = 0; i < op->n_paulis(); ++i) {
    pauli_index_.insert({mat_offset + i, v, i});
    pauli_ac_.col(mat_offset + i).setZero();
    for (const std::pair<const Qubit, Pauli>& qp : op->port(i).string.map)
      qubits_.insert(qp.first);
  }
  // Find ancestors in the anticommutation matrix
  if (op->get_type() == PGOpType::InputRow) {
    const PGInputRow& in = dynamic_cast<const PGInputRow&>(*op);
    QubitPauliTensor tensor = in.get_tensor();
    QubitPauliTensor input_tensor = in.get_input_tensor();
    for (const PGPauli& prev_pauli : pauli_index_.get<TagID>()) {
      PGOp_ptr other_op = c_graph_[prev_pauli.vert];
      if (other_op->get_type() == PGOpType::InputRow) {
        const PGInputRow& other_in = dynamic_cast<const PGInputRow&>(*other_op);
        if (tensor.commutes_with(other_in.get_tensor()) ^
            input_tensor.commutes_with(other_in.get_input_tensor()))
          throw PGError(
              "Cannot insert InputRow into PauliGraph - does not commute with "
              "existing InputRow");
      } else if (!tensor.commutes_with(other_op->port(prev_pauli.port)))
        throw PGError(
            "Cannot insert InputRow into PauliGraph - does not commute with "
            "existing operations");
      pauli_ac_(mat_offset, prev_pauli.index) = false;
    }
  } else if (op->get_type() == PGOpType::OutputRow) {
    const PGOutputRow& out = dynamic_cast<const PGOutputRow&>(*op);
    QubitPauliTensor tensor = out.get_tensor();
    QubitPauliTensor output_tensor = out.get_output_tensor();
    for (const PGPauli& prev_pauli : pauli_index_.get<TagID>()) {
      PGOp_ptr other_op = c_graph_[prev_pauli.vert];
      if (other_op->get_type() == PGOpType::OutputRow) {
        const PGOutputRow& other_out =
            dynamic_cast<const PGOutputRow&>(*other_op);
        if (tensor.commutes_with(other_out.get_tensor()) ^
            output_tensor.commutes_with(other_out.get_output_tensor()))
          throw PGError(
              "Cannot insert OutputRow into PauliGraph - does not commute with "
              "existing OutputRow");
        pauli_ac_(mat_offset, prev_pauli.index) = false;
      } else {
        pauli_ac_(mat_offset, prev_pauli.index) =
            !tensor.commutes_with(other_op->port(prev_pauli.port));
      }
    }
  } else {
    std::vector<QubitPauliTensor> active = op->active_paulis();
    for (unsigned i = 0; i < active.size(); ++i) {
      for (const PGPauli& prev_pauli : pauli_index_.get<TagID>()) {
        PGOp_ptr other_op = c_graph_[prev_pauli.vert];
        QubitPauliTensor other_pauli = other_op->port(prev_pauli.port);
        if (other_op->get_type() == PGOpType::OutputRow) {
          if (active.at(i).commutes_with(other_pauli))
            throw PGError(
                "Cannot insert operation into PauliGraph - does not commute "
                "with existing OutputRow");
          pauli_ac_(mat_offset + i, prev_pauli.index) = false;
        } else {
          pauli_ac_(mat_offset + i, prev_pauli.index) =
              !active.at(i).commutes_with(other_pauli);
        }
      }
    }
  }
  // Find classical predecessors
  std::unordered_set<PGVert> c_preds;
  for (const Bit& b : op->read_bits()) {
    auto prev = last_writes_.find(b);
    if (prev != last_writes_.end()) c_preds.insert(prev->second);
    last_reads_[b].insert(v);
  }
  for (const Bit& b : op->write_bits()) {
    auto prev = last_writes_.find(b);
    if (prev != last_writes_.end()) c_preds.insert(prev->second);
    last_writes_[b] = v;
    auto reads = last_reads_.find(b);
    if (reads != last_reads_.end()) {
      for (const PGVert& read : reads->second) c_preds.insert(read);
      last_reads_.erase(reads);
    }
  }
  // If op both reads and writes to the same bit, it may have been inserted into
  // c_preds
  auto self_loop = c_preds.find(v);
  if (self_loop != c_preds.end()) c_preds.erase(self_loop);
  for (const PGVert& pred : c_preds) boost::add_edge(pred, v, c_graph_);
  return v;
}

void PauliGraph::multiply_strings(
    unsigned source_r, unsigned target_r, Complex coeff) {
  PGPauli source_pgp = *pauli_index_.get<TagID>().find(source_r);
  PGOp_ptr source_op = c_graph_[source_pgp.vert];
  PGPauli target_pgp = *pauli_index_.get<TagID>().find(target_r);
  PGOp_ptr target_op = c_graph_[target_pgp.vert];
  // Update strings in PGOps
  target_op->port(target_pgp.port) = coeff * source_op->port(source_pgp.port) *
                                     target_op->port(target_pgp.port);
  // Update anticommutation matrix
  for (unsigned i = 0; i < pauli_ac_.rows(); ++i) {
    pauli_ac_(i, target_pgp.index) =
        pauli_ac_(i, source_pgp.index) ^ pauli_ac_(i, target_pgp.index);
    pauli_ac_(target_pgp.index) =
        pauli_ac_(source_pgp.index, i) ^ pauli_ac_(target_pgp.index, i);
  }
  pauli_ac_(target_pgp.index, target_pgp.index) = false;
}

void PauliGraph::verify() const {
  // Check validity of the graphs by finding a mutual topological ordering
  std::unordered_set<PGVert> consumed;
  std::map<Bit, PGVert> previous_write;
  std::map<Bit, std::unordered_set<PGVert>> previous_reads;
  std::list<PGOp_ptr> inputs;
  std::list<PGOp_ptr> outputs;
  bool found_more = true;
  while (found_more) {
    found_more = false;
    BGL_FORALL_VERTICES(v, c_graph_, PGClassicalGraph) {
      if (consumed.find(v) != consumed.end()) continue;
      bool initial = true;
      auto in_edge_range = boost::in_edges(v, c_graph_);
      for (auto it = in_edge_range.first; it != in_edge_range.second; ++it) {
        if (consumed.find(boost::source(*it, c_graph_)) == consumed.end()) {
          initial = false;
          break;
        }
      }
      if (!initial) continue;
      auto range = pauli_index_.get<TagOp>().equal_range(v);
      for (auto it = range.first; it != range.second; ++it) {
        for (const PGPauli& c_pauli : pauli_index_.get<TagID>()) {
          if (pauli_ac_(it->index, c_pauli.index) &&
              (consumed.find(c_pauli.vert) == consumed.end())) {
            initial = false;
            break;
          }
        }
      }
      if (!initial) continue;
      // Found a valid next vertex, so check its relative validity
      found_more = true;
      consumed.insert(v);
      // Check Classical history contains exactly predecessor hazards for all
      // classical bits and all active bits are registered
      PGOp_ptr op = c_graph_[v];
      std::unordered_set<PGVert> justified_preds;
      for (const Bit& b : op->write_bits()) {
        if (bits_.find(b) == bits_.end())
          throw PGError("PGOp writes to unregistered bit: " + op->get_name());
        auto w_found = previous_write.find(b);
        if (w_found != previous_write.end()) {
          if (!boost::edge(w_found->second, v, c_graph_).second)
            throw PGError(
                "No edge in PGClassicalGraph for RAW dependency on bit " +
                b.repr() + " between " + c_graph_[w_found->second]->get_name() +
                " and " + op->get_name());
          justified_preds.insert(w_found->second);
        }
        previous_reads[b].insert(v);
      }
      for (const Bit& b : op->write_bits()) {
        if (bits_.find(b) == bits_.end())
          throw PGError("PGOp writes to unregistered bit: " + op->get_name());
        auto w_found = previous_write.find(b);
        if (w_found != previous_write.end()) {
          if (!boost::edge(w_found->second, v, c_graph_).second)
            throw PGError(
                "No edge in PGClassicalGraph for WAW dependency on bit " +
                b.repr() + " between " + c_graph_[w_found->second]->get_name() +
                " and " + op->get_name());
          justified_preds.insert(w_found->second);
        }
        auto r_found = previous_reads.find(b);
        if (r_found != previous_reads.end()) {
          for (const PGVert& u : r_found->second) {
            if (!boost::edge(u, v, c_graph_).second)
              throw PGError(
                  "No edge in PGClassicalGraph for WAR dependency on bit " +
                  b.repr() + " between " + c_graph_[u]->get_name() + " and " +
                  op->get_name());
            justified_preds.insert(u);
          }
          r_found->second.clear();
        }
        previous_write[b] = v;
      }
      for (auto it = in_edge_range.first; it != in_edge_range.second; ++it) {
        PGVert s = boost::source(*it, c_graph_);
        if (justified_preds.find(s) == justified_preds.end())
          throw PGError(
              "Edge in PGClassicalGraph despite no dependency between " +
              c_graph_[s]->get_name() + " and " + op->get_name());
      }
      // Check Pauli history contains all anti-commuting terms and all active
      // qubits are registered
      std::vector<QubitPauliTensor> paulis = op->active_paulis();
      for (auto it = range.first; it != range.second; ++it) {
        const QubitPauliTensor& tensor = paulis.at(it->port);
        for (const std::pair<const Qubit, Pauli>& qp : tensor.string.map) {
          if (qubits_.find(qp.first) == qubits_.end())
            throw PGError(
                "PGOp interacts with unregistered qubit: " + op->get_name());
        }
        if (op->get_type() == PGOpType::InputRow) {
          for (unsigned c = 0; c < pauli_ac_.cols(); ++c) {
            if (pauli_ac_(it->index, c))
              throw PGError(
                  "PauliGraph input row has predecessors in anticommutation "
                  "matrix: " +
                  op->get_name());
          }
          const PGInputRow& in = dynamic_cast<const PGInputRow&>(*op);
          for (const PGOp_ptr& other : inputs) {
            const PGInputRow& other_in =
                dynamic_cast<const PGInputRow&>(*other);
            if (in.get_tensor().commutes_with(other_in.get_tensor()) !=
                in.get_input_tensor().commutes_with(
                    other_in.get_input_tensor()))
              throw PGError(
                  "PauliGraph contains anticommuting input rows: " +
                  op->get_name() + " and " + other->get_name());
          }
          inputs.push_back(op);
        } else if (op->get_type() == PGOpType::OutputRow) {
          for (const PGPauli& c_pauli : pauli_index_.get<TagID>()) {
            if ((consumed.find(c_pauli.vert) != consumed.end()) &&
                (pauli_ac_(it->index, c_pauli.index) !=
                 (c_graph_[it->vert]->get_type() != PGOpType::OutputRow &&
                  !tensor.commutes_with(c_graph_[it->vert]->port(it->port)))))
              throw PGError(
                  "PauliGraph anticommutation matrix is missing a link "
                  "between " +
                  c_graph_[it->vert]->get_name() + " and " + op->get_name());
          }
          for (unsigned r = 0; r < pauli_ac_.rows(); ++r) {
            if (pauli_ac_(r, it->index))
              throw PGError(
                  "PauliGraph output row has successors in anticommutation "
                  "matrix: " +
                  op->get_name());
          }
          const PGOutputRow& out = dynamic_cast<const PGOutputRow&>(*op);
          for (const PGOp_ptr& other : outputs) {
            const PGOutputRow& other_out =
                dynamic_cast<const PGOutputRow&>(*other);
            if (out.get_tensor().commutes_with(other_out.get_tensor()) !=
                out.get_output_tensor().commutes_with(
                    other_out.get_output_tensor()))
              throw PGError(
                  "PauliGraph contains anticommuting output rows: " +
                  op->get_name() + " and " + other->get_name());
          }
          outputs.push_back(op);
        } else {
          for (const PGPauli& c_pauli : pauli_index_.get<TagID>()) {
            if ((consumed.find(c_pauli.vert) != consumed.end()) &&
                (pauli_ac_(it->index, c_pauli.index) ==
                 tensor.commutes_with(
                     c_graph_[c_pauli.vert]->port(c_pauli.port))))
              throw PGError(
                  "PauliGraph anticommutation matrix is missing a link "
                  "between " +
                  c_graph_[c_pauli.vert]->get_name() + " and " +
                  op->get_name());
          }
        }
      }
    }
  }
}

}  // namespace pg
}  // namespace tket
