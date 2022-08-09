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

bool PGOp::is_equal(const PGOp&) const { return true; }

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

bit_vector_t PGOp::read_bits() const { return {}; }

bit_vector_t PGOp::write_bits() const { return {}; }

/**
 * PGRotation Implementation
 */

const QubitPauliTensor& PGRotation::get_tensor() const { return tensor_; }

const Expr& PGRotation::get_angle() const { return angle_; }

PGRotation::PGRotation(const QubitPauliTensor& tensor, const Expr& angle) : PGOp(PGOpType::Rotation), tensor_(tensor), angle_(angle) {
  if (tensor_.coeff == -1.) {
    angle_ *= -1.;
    tensor_.coeff = 1.;
  }
  else if (tensor_.coeff != 1.) throw PGError("Invalid coefficient in tensor for PauliGraph rotation");
}

SymSet PGRotation::free_symbols() const { return expr_free_symbols(angle_); }

PGOp_ptr PGRotation::symbol_substitution(const SymEngine::map_basic_basic& sub_map) const {
  return std::make_shared<const PGRotation>(tensor_, angle_.subs(sub_map));
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

std::vector<QubitPauliTensor> PGRotation::active_paulis() const { return {tensor_}; }

/**
 * PGCliffordRot Implementation
 */

const QubitPauliTensor& PGCliffordRot::get_tensor() const { return tensor_; }

unsigned PGCliffordRot::get_angle() const { return angle_; }

PGCliffordRot::PGCliffordRot(const QubitPauliTensor& tensor, unsigned angle) : PGOp(PGOpType::CliffordRot), tensor_(tensor), angle_(angle) {
  if (tensor_.coeff == -1.) {
    angle_ = (4 - angle) % 4;
    tensor_.coeff = 1.;
  }
  else if (tensor_.coeff != 1.) throw PGError("Invalid coefficient in tensor for PauliGraph Clifford rotation");
}

SymSet PGCliffordRot::free_symbols() const { return {}; }

PGOp_ptr PGCliffordRot::symbol_substitution(const SymEngine::map_basic_basic&) const { return PGOp_ptr(); }

std::string PGCliffordRot::get_name(bool) const {
  std::stringstream str;
  str << "ClfRot(" << tensor_.to_str() << "; " << (angle_ * 0.5) << ")";
  return str.str();
}

bool PGCliffordRot::is_equal(const PGOp& op_other) const {
  const PGCliffordRot& other = dynamic_cast<const PGCliffordRot&>(op_other);
  return (tensor_ == other.tensor_) && (angle_ % 4 == other.angle_ % 4);
}

std::vector<QubitPauliTensor> PGCliffordRot::active_paulis() const { return {tensor_}; }

/**
 * PGMeasure Implementation
 */

const QubitPauliTensor& PGMeasure::get_tensor() const { return tensor_; }

const Bit& PGMeasure::get_target() const { return target_; }

PGMeasure::PGMeasure(const QubitPauliTensor& tensor, const Bit& target) : PGOp(PGOpType::Measure), tensor_(tensor), target_(target) {
  if (tensor_.coeff != 1. && tensor_.coeff != -1.) throw PGError("Invalid coefficient in tensor for PauliGraph measurement");
}

SymSet PGMeasure::free_symbols() const { return {}; }

PGOp_ptr PGMeasure::symbol_substitution(const SymEngine::map_basic_basic&) const { return PGOp_ptr(); }

std::string PGMeasure::get_name(bool) const {
  std::stringstream str;
  str << "Meas(" << tensor_.to_str() << " -> " << target_.repr() << ")";
  return str.str();
}

bool PGMeasure::is_equal(const PGOp& op_other) const {
  const PGMeasure& other = dynamic_cast<const PGMeasure&>(op_other);
  return (tensor_ == other.tensor_) && (target_ == other.target_);
}

std::vector<QubitPauliTensor> PGMeasure::active_paulis() const { return {tensor_}; }

bit_vector_t PGMeasure::write_bits() const { return {target_}; }

/**
 * PGDecoherence Implementation
 */

const QubitPauliTensor& PGDecoherence::get_tensor() const { return tensor_; }

PGDecoherence::PGDecoherence(const QubitPauliTensor& tensor) : PGOp(PGOpType::Decoherence), tensor_(tensor) {
  if (tensor_.coeff == -1.) tensor_.coeff = 1.;
  else if (tensor_.coeff != 1.) throw PGError("Invalid coefficient in tensor for PauliGraph decoherence");
}

SymSet PGDecoherence::free_symbols() const { return {}; }

PGOp_ptr PGDecoherence::symbol_substitution(const SymEngine::map_basic_basic&) const { return PGOp_ptr(); }

std::string PGDecoherence::get_name(bool) const {
  std::stringstream str;
  str << "Deco(" << tensor_.to_str() << ")";
  return str.str();
}

bool PGDecoherence::is_equal(const PGOp& op_other) const {
  const PGDecoherence& other = dynamic_cast<const PGDecoherence&>(op_other);
  return tensor_ == other.tensor_;
}

std::vector<QubitPauliTensor> PGDecoherence::active_paulis() const { return {tensor_}; }

/**
 * PGReset Implementation
 */

const QubitPauliTensor& PGReset::get_stab() const { return stab_; }

const QubitPauliTensor& PGReset::get_destab() const { return destab_; }

PGReset::PGReset(const QubitPauliTensor& stab, const QubitPauliTensor& destab) : PGOp(PGOpType::Reset), stab_(stab), destab_(destab) {
  if (destab_.coeff == -1.) destab_.coeff = 1.;
  else if (destab_.coeff != 1. || (stab_.coeff != 1. && stab_.coeff != -1.)) throw PGError("Invalid coefficient in tensor for PauliGraph reset");
}

SymSet PGReset::free_symbols() const { return {}; }

PGOp_ptr PGReset::symbol_substitution(const SymEngine::map_basic_basic&) const { return PGOp_ptr(); }

std::string PGReset::get_name(bool) const {
  std::stringstream str;
  str << "Reset(" << stab_.to_str() << "; " << destab_.to_str() << ")";
  return str.str();
}

bool PGReset::is_equal(const PGOp& op_other) const {
  const PGReset& other = dynamic_cast<const PGReset&>(op_other);
  return (stab_ == other.stab_) && (destab_ == other.destab_);
}

std::vector<QubitPauliTensor> PGReset::active_paulis() const { return {stab_, destab_}; }

/**
 * PGConditional Implementation
 */

PGOp_ptr PGConditional::get_inner_op() const { return inner_; }

const bit_vector_t& PGConditional::get_args() const { return args_; }

unsigned PGConditional::get_value() const { return value_; }

PGConditional::PGConditional(PGOp_ptr inner, const bit_vector_t& args, unsigned value) : PGOp(PGOpType::Conditional), inner_(inner), args_(args), value_(value) {}

SymSet PGConditional::free_symbols() const { return inner_->free_symbols(); }

PGOp_ptr PGConditional::symbol_substitution(const SymEngine::map_basic_basic& sub_map) const {
  PGOp_ptr inner_sub = inner_->symbol_substitution(sub_map);
  if (inner_sub) return std::make_shared<const PGConditional>(inner_sub, args_, value_);
  else return PGOp_ptr();
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
  return (value_ == other.value_) && (args_ == other.args_) && (*inner_ == *other.inner_);
}

std::vector<QubitPauliTensor> PGConditional::active_paulis() const { return inner_->active_paulis(); }

bit_vector_t PGConditional::read_bits() const {
  bit_vector_t bits = inner_->read_bits();
  bits.insert(bits.end(), args_.begin(), args_.end());
  return bits;
}

bit_vector_t PGConditional::write_bits() const { return inner_->write_bits(); }

/**
 * PauliGraph Implementation
 */

PauliGraph::PauliGraph(unsigned n) : initial_({}), final_({}) {
  MatrixXb initial_xmat = MatrixXb::Zero(2*n, n);
  MatrixXb initial_zmat = MatrixXb::Zero(2*n, n);
  MatrixXb final_xmat = MatrixXb::Zero(2*n, n);
  MatrixXb final_zmat = MatrixXb::Zero(2*n, n);
  for (unsigned i = 0; i < n; ++i) {
    initial_cols_.insert({Qubit(i), i});
    initial_rows_.insert({{Qubit(i), TableauRowType::ZRow}, 2*i});
    initial_zmat(2*i, i) = true;
    initial_rows_.insert({{Qubit(i), TableauRowType::XRow}, 2*i + 1});
    initial_xmat(2*i + 1, i) = true;
    final_cols_.insert({Qubit(i), i});
    final_rows_.insert({{Qubit(i), TableauRowType::ZRow}, 2*i});
    final_zmat(2*i, i) = true;
    final_rows_.insert({{Qubit(i), TableauRowType::XRow}, 2*i + 1});
    final_xmat(2*i + 1, i) = true;
  }
  initial_ = SymplecticTableau(initial_xmat, initial_zmat, VectorXb::Zero(2 * n));
  final_ = SymplecticTableau(final_xmat, final_zmat, VectorXb::Zero(2 * n));
}

PauliGraph::PauliGraph(const qubit_vector_t& qbs, const bit_vector_t& bits) : initial_({}), final_({}), bits_(bits) {
  unsigned n = qbs.size();
  MatrixXb initial_xmat = MatrixXb::Zero(2*n, n);
  MatrixXb initial_zmat = MatrixXb::Zero(2*n, n);
  MatrixXb final_xmat = MatrixXb::Zero(2*n, n);
  MatrixXb final_zmat = MatrixXb::Zero(2*n, n);
  for (unsigned i = 0; i < n; ++i) {
    initial_cols_.insert({qbs.at(i), i});
    initial_rows_.insert({{qbs.at(i), TableauRowType::ZRow}, 2*i});
    initial_zmat(2*i, i) = true;
    initial_rows_.insert({{qbs.at(i), TableauRowType::XRow}, 2*i + 1});
    initial_xmat(2*i + 1, i) = true;
    final_cols_.insert({qbs.at(i), i});
    final_rows_.insert({{qbs.at(i), TableauRowType::ZRow}, 2*i});
    final_zmat(2*i, i) = true;
    final_rows_.insert({{qbs.at(i), TableauRowType::XRow}, 2*i + 1});
    final_xmat(2*i + 1, i) = true;
  }
  initial_ = SymplecticTableau(initial_xmat, initial_zmat, VectorXb::Zero(2 * n));
  final_ = SymplecticTableau(final_xmat, final_zmat, VectorXb::Zero(2 * n));
}

void PauliGraph::to_graphviz(std::ostream &out) const {
  out << "digraph G {\n";

  std::map<PGVert, unsigned> index_map;
  out << "0 [label = \"" << initial_ << "\"];\n";
  out << "1 [label = \"" << final_ << "\"];\n";
  unsigned i = 2;
  BGL_FORALL_VERTICES(v, graph_, PGDAG) {
    index_map.insert({v, i});
    out << i << " [label = \"" << graph_[v].op_->get_name() << "\"];\n";
    ++i;
  }

  for (const PGVert& sv : start_line_.get<TagSeq>()) {
    out << 0 << " --> " << index_map.at(sv) << ";\n";
  }
  for (const PGVert& ev : end_line_.get<TagSeq>()) {
    out << index_map.at(ev) << " --> 1;\n";
  }
  if (start_line_.empty()) out << "0 --> 1;\n";
  BGL_FORALL_EDGES(e, graph_, PGDAG) {
    PGVert vs = source(e);
    PGVert vt = target(e);
    out << index_map.at(vs) << " --> " << index_map.at(vt) << ";\n";
  }

  out << "}";
}

unsigned PauliGraph::n_vertices() const { return boost::num_vertices(graph_); }

PGVertSet PauliGraph::get_successors(const PGVert& vert) const {
  PGVertSet succs;
  for (auto iter = boost::adjacent_vertices(vert, graph_); iter.first != iter.second; ++iter.first) {
    succs.insert(*iter.first);
  }
  return succs;
}

PGVertSet PauliGraph::get_predecessors(const PGVert& vert) const {
  PGVertSet preds;
  for (auto iter = boost::inv_adjacent_vertices(vert, graph_); iter.first != iter.second; ++iter.first) {
    preds.insert(*iter.first);
  }
  return preds;
}

PGEdgeSet PauliGraph::get_in_edges(const PGVert& vert) const {
  PGEdgeSet ins;
  for (auto iter = boost::in_edges(vert, graph_); iter.first != iter.second; ++iter.first) {
    ins.insert(*iter.first);
  }
  return ins;
}

PGEdgeSet PauliGraph::get_out_edges(const PGVert& vert) const {
  PGEdgeSet outs;
  for (auto iter = boost::out_edges(vert, graph_); iter.first != iter.second; ++iter.first) {
    outs.insert(*iter.first);
  }
  return outs;
}

PGVert PauliGraph::source(const PGEdge& edge) const {
  return boost::source(edge, graph_);
}

PGVert PauliGraph::target(const PGEdge& edge) const {
  return boost::target(edge, graph_);
}

void PauliGraph::add_vertex_at_end(PGOp_ptr op) {
  PGVertSet to_search = end_line_;
  PGVertSet commuted;
  PGVert new_vert = boost::add_vertex(graph_);
  graph_[new_vert] = {op};
  while (!to_search.empty()) {
    // Get next candidate parent
    PGVertSet::index<TagSeq>::type& to_search_seq = to_search.get<TagSeq>();
    PGVert to_compare = to_search_seq.front();
    to_search_seq.pop_front();

    // Check that we have already commuted past all of its children
    bool ready = true;
    PGVertSet::index<TagKey>::type& commuted_key = commuted.get<TagKey>();
    for (const PGVert &child : get_successors(to_compare)) {
      if (commuted_key.find(child) == commuted_key.end()) {
        ready = false;
        break;
      }
    }
    if (!ready) continue;

    // Check if we can commute past it
    PGOp_ptr compare_op = graph_[to_compare].op_;
    if (op->commutes_with(*compare_op)) {
      // Commute and continue searching
      PGVertSet preds = get_predecessors(to_compare);
      to_search.insert(preds.begin(), preds.end());
      commuted.insert(to_compare);
    }
    else {
      // Does not commute - add dependency edge
      boost::add_edge(to_compare, new_vert, graph_);
      end_line_.erase(to_compare);
    }
  }
  end_line_.insert(new_vert);
  if (get_predecessors(new_vert).empty()) start_line_.insert(new_vert);
}

PauliGraph::TopSortIterator::TopSortIterator() : pg_(nullptr), current_vert_(boost::graph_traits<PGDAG>::null_vertex()) {}

PauliGraph::TopSortIterator::TopSortIterator(const PauliGraph& pg) {
  if (pg.start_line_.empty()) {
    current_vert_ = boost::graph_traits<PGDAG>::null_vertex();
    return;
  }
  pg_ = &pg;
  search_set_ = pg.start_line_;
  current_vert_ = search_set_.get<TagSeq>().front();
  search_set_.get<TagSeq>().pop_front();
  visited_ = {current_vert_};
  for (const PGVert& s : pg_->get_successors(current_vert_)) search_set_.insert(s);
}

const PGVert& PauliGraph::TopSortIterator::operator*() const { return current_vert_; }

const PGVert* PauliGraph::TopSortIterator::operator->() const { return &current_vert_; }

bool PauliGraph::TopSortIterator::operator==(const TopSortIterator& other) const { return current_vert_ == other.current_vert_; }

bool PauliGraph::TopSortIterator::operator!=(const TopSortIterator& other) const { return !(*this == other); }

PauliGraph::TopSortIterator PauliGraph::TopSortIterator::operator++(int) {
  PauliGraph::TopSortIterator it = *this;
  ++*this;
  return it;
}

PauliGraph::TopSortIterator& PauliGraph::TopSortIterator::operator++() {
  bool found_next = false;
  PGVertSet::index<TagSeq>::type& search_seq = search_set_.get<TagSeq>();
  while (!found_next && !search_set_.empty()) {
    current_vert_ = search_seq.front();
    search_seq.pop_front();

    // Check that we have visited all parents
    found_next = true;
    for (const PGVert& p : pg_->get_predecessors(current_vert_)) {
      if (visited_.find(p) == visited_.end()) {
        found_next = false;
        break;
      }
    }
  }
  if (found_next) {
    visited_.insert(current_vert_);
    for (const PGVert& s : pg_->get_successors(current_vert_)) {
      search_set_.insert(s);
    }
  }
  else {
    *this = TopSortIterator();
  }
  return *this;
}

PauliGraph::TopSortIterator PauliGraph::begin() const { return TopSortIterator(*this); }

PauliGraph::TopSortIterator PauliGraph::end() const { return TopSortIterator(); }

}  // namespace pg
}  // namespace tket
