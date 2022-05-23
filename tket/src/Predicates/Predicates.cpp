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

#include "Predicates.hpp"

#include "Gate/Gate.hpp"
#include "Mapping/Verification.hpp"
#include "OpType/OpTypeFunctions.hpp"
#include "Placement/Placement.hpp"
#include "Utils/UnitID.hpp"

namespace tket {

template <typename T>
static bool auto_implication(const T&, const Predicate& other) {
  try {
    (void)dynamic_cast<const T&>(other);
  } catch (const std::bad_cast&) {
    throw IncorrectPredicate(
        "Cannot compare predicates of different subclasses");
  }
  return true;
}

template <typename T>
static PredicatePtr auto_meet(const T&, const Predicate& other) {
  try {
    (void)dynamic_cast<const T&>(other);
    PredicatePtr pp = std::make_shared<T>();
    return pp;
  } catch (const std::bad_cast&) {
    throw IncorrectPredicate(
        "Cannot compare predicates of different subclasses");
  }
}

template <typename T>
static std::string auto_name(const T&) {
  return predicate_name(typeid(T));
}

const std::string& predicate_name(std::type_index idx) {
  static const std::map<std::type_index, std::string> predicate_names = {
#define SET_PRED_NAME(a) {typeid(a), #a}
      SET_PRED_NAME(CliffordCircuitPredicate),
      SET_PRED_NAME(ConnectivityPredicate),
      SET_PRED_NAME(DefaultRegisterPredicate),
      SET_PRED_NAME(DirectednessPredicate),
      SET_PRED_NAME(GateSetPredicate),
      SET_PRED_NAME(MaxNQubitsPredicate),
      SET_PRED_NAME(MaxTwoQubitGatesPredicate),
      SET_PRED_NAME(NoBarriersPredicate),
      SET_PRED_NAME(NoClassicalBitsPredicate),
      SET_PRED_NAME(NoClassicalControlPredicate),
      SET_PRED_NAME(NoFastFeedforwardPredicate),
      SET_PRED_NAME(NoMidMeasurePredicate),
      SET_PRED_NAME(NoSymbolsPredicate),
      SET_PRED_NAME(GlobalPhasedXPredicate),
      SET_PRED_NAME(NoWireSwapsPredicate),
      SET_PRED_NAME(PlacementPredicate),
      SET_PRED_NAME(UserDefinedPredicate)};
#undef SET_PRED_NAME
  return predicate_names.at(idx);
}

/////////////////////
// PREDICATE METHODS//
/////////////////////

bool GateSetPredicate::verify(const Circuit& circ) const {
  BGL_FORALL_VERTICES(v, circ.dag, DAG) {
    Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
    OpDesc desc = op->get_desc();
    if (desc.is_meta()) continue;
    OpType type = op->get_type();
    if (type == OpType::Conditional) {
      const Conditional& cond = static_cast<const Conditional&>(*op);
      type = cond.get_op()->get_type();
    }
    if (!find_in_set(type, allowed_types_)) return false;
  }
  return true;
}

bool GateSetPredicate::implies(const Predicate& other) const {
  try {
    const GateSetPredicate& other_p =
        dynamic_cast<const GateSetPredicate&>(other);
    for (const OpType& ot : allowed_types_) {
      if (other_p.allowed_types_.find(ot) == other_p.allowed_types_.end())
        return false;
    }
    return true;
  } catch (const std::bad_cast&) {
    throw IncorrectPredicate(
        "Cannot compare predicates of different subclasses");
  }
}

PredicatePtr GateSetPredicate::meet(const Predicate& other) const {
  try {
    const GateSetPredicate& other_p =
        dynamic_cast<const GateSetPredicate&>(other);
    OpTypeSet new_set;
    for (const OpType& ot : allowed_types_) {
      if (other_p.allowed_types_.find(ot) != other_p.allowed_types_.end()) {
        new_set.insert(ot);
      }
    }
    PredicatePtr pp = std::make_shared<GateSetPredicate>(new_set);
    return pp;
  } catch (const std::bad_cast&) {
    throw IncorrectPredicate(
        "Cannot compare predicates of different subclasses");
  }
}

std::string GateSetPredicate::to_string() const {
  std::string str = auto_name(*this) + ":{ ";
  for (const OpType& ot : allowed_types_) {
    str += (optypeinfo().find(ot)->second.name + " ");
  }
  str += "}";
  return str;
}

bool NoClassicalControlPredicate::verify(const Circuit& circ) const {
  BGL_FORALL_VERTICES(v, circ.dag, DAG) {
    Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
    OpType ot = op->get_type();
    if (ot == OpType::Conditional)
      return false;
    else if (ot == OpType::CircBox || ot == OpType::CustomGate) {
      const Box& box = static_cast<const Box&>(*op);
      if (!verify(*box.to_circuit())) return false;
    }
  }
  return true;
}

bool NoClassicalControlPredicate::implies(const Predicate& other) const {
  return auto_implication(*this, other);
}

PredicatePtr NoClassicalControlPredicate::meet(const Predicate& other) const {
  return auto_meet(*this, other);
}

std::string NoClassicalControlPredicate::to_string() const {
  return auto_name(*this);
}

static bool fast_feed_forward_helper(
    const Command& com, unit_set_t& unset_bits) {
  // Allows conditionals from unset_bits
  // Encountering a measurement removed the bits from unset_bits
  // Returns whether or not a feed-forward conditional is found
  // Applies recursively for CircBoxes
  if (com.get_op_ptr()->get_type() == OpType::Conditional) {
    const Conditional& cond =
        static_cast<const Conditional&>(*com.get_op_ptr());
    unit_vector_t all_args = com.get_args();
    unit_vector_t::iterator arg_it = all_args.begin();
    for (unsigned i = 0; i < cond.get_width(); ++i) {
      if (unset_bits.find(*arg_it) == unset_bits.end()) return false;
      ++arg_it;
    }
    unit_vector_t new_args = {arg_it, all_args.end()};
    Command new_com = {cond.get_op(), new_args};
    return fast_feed_forward_helper(new_com, unset_bits);
  } else if (
      com.get_op_ptr()->get_type() == OpType::CircBox ||
      com.get_op_ptr()->get_type() == OpType::CustomGate) {
    const Box& box = static_cast<const Box&>(*com.get_op_ptr());
    unit_map_t interface;
    unit_set_t inner_set;
    unsigned i = 0;
    for (const Bit& b : com.get_bits()) {
      Bit inner_bit(i);
      interface.insert({Bit(i), b});
      if (unset_bits.find(b) != unset_bits.end()) {
        inner_set.insert(inner_bit);
      }
      ++i;
    }
    for (const Command& c : *box.to_circuit()) {
      if (!fast_feed_forward_helper(c, inner_set)) return false;
    }
    for (const std::pair<const UnitID, UnitID>& pair : interface) {
      if (inner_set.find(pair.first) == inner_set.end())
        unset_bits.erase(pair.second);
    }
  } else if (com.get_op_ptr()->get_type() == OpType::Measure) {
    unset_bits.erase(com.get_args().at(1));
  }
  return true;
}

bool NoFastFeedforwardPredicate::verify(const Circuit& circ) const {
  if (circ.n_bits() == 0) return true;
  bit_vector_t all_bits = circ.all_bits();
  unit_set_t unset_bits = {all_bits.begin(), all_bits.end()};
  for (const Command& com : circ) {
    if (!fast_feed_forward_helper(com, unset_bits)) return false;
  }
  return true;
}

bool NoFastFeedforwardPredicate::implies(const Predicate& other) const {
  return auto_implication(*this, other);
}

PredicatePtr NoFastFeedforwardPredicate::meet(const Predicate& other) const {
  return auto_meet(*this, other);
}

std::string NoFastFeedforwardPredicate::to_string() const {
  return auto_name(*this);
}

bool NoClassicalBitsPredicate::verify(const Circuit& circ) const {
  // for classical edges, we currently require classical input
  return (circ.n_bits() == 0);
}

bool NoClassicalBitsPredicate::implies(const Predicate& other) const {
  return auto_implication(*this, other);
}

PredicatePtr NoClassicalBitsPredicate::meet(const Predicate& other) const {
  return auto_meet(*this, other);
}

std::string NoClassicalBitsPredicate::to_string() const {
  return auto_name(*this);
}

bool NoWireSwapsPredicate::verify(const Circuit& circ) const {
  std::map<UnitID, QPathDetailed> paths = circ.all_unit_paths();
  for (const std::pair<const UnitID, QPathDetailed>& path : paths) {
    UnitID out_id = circ.get_id_from_out(path.second.back().first);
    if (path.first != out_id) return false;
  }
  return true;
}

bool NoWireSwapsPredicate::implies(const Predicate& other) const {
  return auto_implication(*this, other);
}

PredicatePtr NoWireSwapsPredicate::meet(const Predicate& other) const {
  return auto_meet(*this, other);
}

std::string NoWireSwapsPredicate::to_string() const { return auto_name(*this); }

bool MaxTwoQubitGatesPredicate::verify(const Circuit& circ) const {
  BGL_FORALL_VERTICES(v, circ.dag, DAG) {
    if (circ.get_OpType_from_Vertex(v) == OpType::Barrier) continue;
    if (circ.n_in_edges_of_type(v, EdgeType::Quantum) > 2) return false;
  }
  return true;
}

bool MaxTwoQubitGatesPredicate::implies(const Predicate& other) const {
  return auto_implication(*this, other);
}

PredicatePtr MaxTwoQubitGatesPredicate::meet(const Predicate& other) const {
  return auto_meet(*this, other);
}

std::string MaxTwoQubitGatesPredicate::to_string() const {
  return auto_name(*this);
}

bool PlacementPredicate::verify(const Circuit& circ) const {
  for (const Qubit& qb : circ.all_qubits()) {
    if (qb.reg_name() == Placement::unplaced_reg()) continue;
    if (nodes_.find(Node(qb)) == nodes_.end()) return false;
  }
  return true;
}

bool PlacementPredicate::implies(const Predicate& other) const {
  const PlacementPredicate other_p =
      dynamic_cast<const PlacementPredicate&>(other);
  for (const Node& node : nodes_) {
    if (other_p.nodes_.find(node) == other_p.nodes_.end()) return false;
  }
  return true;
}

PredicatePtr PlacementPredicate::meet(const Predicate& other) const {
  try {
    const PlacementPredicate& other_c =
        dynamic_cast<const PlacementPredicate&>(other);
    node_set_t nodes;
    for (const Node& node : nodes_) {
      if (other_c.nodes_.find(node) != other_c.nodes_.end()) nodes.insert(node);
    }
    PredicatePtr pp = std::make_shared<PlacementPredicate>(nodes);
    return pp;
  } catch (const std::bad_cast&) {
    throw IncorrectPredicate(
        "Cannot compare predicates of different subclasses");
  }
}

std::string PlacementPredicate::to_string() const {
  std::string str = auto_name(*this) + ":{ ";
  str += ("Nodes: " + std::to_string(nodes_.size()) + " }");
  return str;
}

bool ConnectivityPredicate::verify(const Circuit& circ) const {
  return respects_connectivity_constraints(circ, arch_, false, true);
}

bool ConnectivityPredicate::implies(const Predicate& other) const {
  try {
    const ConnectivityPredicate& other_c =
        dynamic_cast<const ConnectivityPredicate&>(other);
    const Architecture& arc1 = arch_;
    const Architecture& arc2 = other_c.arch_;
    // Check that all nodes in arc1 are in arc2:
    for (const Node& n : arc1.get_all_nodes_vec()) {
      if (!arc2.node_exists(n)) {
        return false;
      }
    }
    // Collect all edges in arc1
    for (auto [n1, n2] : arc1.get_all_edges_vec()) {
      if (!arc2.edge_exists(n1, n2) && !arc2.edge_exists(n2, n1)) {
        return false;  // if not in second architecture, return false
      }
    }
    return true;
  } catch (const std::bad_cast&) {
    throw IncorrectPredicate(
        "Cannot compare predicates of different subclasses");
  }
}

PredicatePtr ConnectivityPredicate::meet(const Predicate& other) const {
  try {
    const ConnectivityPredicate& other_c =
        dynamic_cast<const ConnectivityPredicate&>(other);
    const Architecture& arc1 = arch_;
    const Architecture& arc2 = other_c.arch_;
    std::vector<std::pair<Node, Node>> new_edges;
    // Collect all edges in arc1 which are also in arc2
    for (auto [n1, n2] : arc1.get_all_edges_vec()) {
      if (arc2.edge_exists(n1, n2)) {
        new_edges.push_back({n1, n2});
        new_edges.push_back({n2, n1});
      }
    }
    Architecture arc3(new_edges);
    PredicatePtr pp = std::make_shared<ConnectivityPredicate>(arc3);
    return pp;
  } catch (const std::bad_cast&) {
    throw IncorrectPredicate(
        "Cannot compare predicates of different subclasses");
  }
}

std::string ConnectivityPredicate::to_string() const {
  std::string str = auto_name(*this) + ":{ ";
  str +=
      ("Nodes: " + std::to_string(arch_.n_nodes()) +
       ", Edges: " + std::to_string(arch_.n_connections())) += " }";
  return str;
}

bool DirectednessPredicate::verify(const Circuit& circ) const {
  return respects_connectivity_constraints(circ, arch_, true, false);
}

bool DirectednessPredicate::implies(const Predicate& other) const {
  try {
    const DirectednessPredicate& other_c =
        dynamic_cast<const DirectednessPredicate&>(other);
    const Architecture& arc1 = arch_;
    const Architecture& arc2 = other_c.arch_;
    // Collect all edges in arc1
    for (auto [n1, n2] : arc1.get_all_edges_vec()) {
      // directedness accounted for
      if (!arc2.edge_exists(n1, n2)) {
        return false;  // if not in second architecture, return false
      }
    }
    return true;
  } catch (const std::bad_cast&) {
    throw IncorrectPredicate(
        "Cannot compare predicates of different subclasses");
  }
}

PredicatePtr DirectednessPredicate::meet(const Predicate& other) const {
  try {
    const DirectednessPredicate& other_c =
        dynamic_cast<const DirectednessPredicate&>(other);
    const Architecture& arc1 = arch_;
    const Architecture& arc2 = other_c.arch_;
    std::vector<std::pair<Node, Node>> new_edges;
    // Collect all edges in arc1 which are also in arc2
    for (auto [n1, n2] : arc1.get_all_edges_vec()) {
      // this also accounts for directedness, do we want that?
      if (arc2.edge_exists(n1, n2)) {
        new_edges.push_back({n1, n2});
      }
    }
    Architecture arc3(new_edges);
    PredicatePtr pp = std::make_shared<DirectednessPredicate>(arc3);
    return pp;
  } catch (const std::bad_cast&) {
    throw IncorrectPredicate(
        "Cannot compare predicates of different subclasses");
  }
}

std::string DirectednessPredicate::to_string() const {
  std::string str = auto_name(*this) + ":{ ";
  str +=
      ("Nodes: " + std::to_string(arch_.n_nodes()) +
       ", Edges: " + std::to_string(arch_.n_connections())) += " }";
  return str;
}

bool CliffordCircuitPredicate::verify(const Circuit& circ) const {
  BGL_FORALL_VERTICES(v, circ.dag, DAG) {
    if (!circ.get_Op_ptr_from_Vertex(v)->is_clifford()) return false;
  }
  return true;
}

bool CliffordCircuitPredicate::implies(const Predicate& other) const {
  return auto_implication(*this, other);
}

PredicatePtr CliffordCircuitPredicate::meet(const Predicate& other) const {
  return auto_meet(*this, other);
}

std::string CliffordCircuitPredicate::to_string() const {
  return auto_name(*this);
}

bool UserDefinedPredicate::implies(const Predicate&) const {
  throw IncorrectPredicate(
      "Cannot deduce implication relations of user defined Predicates");
}

PredicatePtr UserDefinedPredicate::meet(const Predicate&) const {
  throw IncorrectPredicate("Cannot find the meet of user defined Predicates");
}

bool UserDefinedPredicate::verify(const Circuit& circ) const {
  return func_(circ);
}

std::string UserDefinedPredicate::to_string() const { return auto_name(*this); }

bool DefaultRegisterPredicate::verify(const Circuit& circ) const {
  return circ.is_simple();
}

bool DefaultRegisterPredicate::implies(const Predicate& other) const {
  return auto_implication(*this, other);
}

PredicatePtr DefaultRegisterPredicate::meet(const Predicate& other) const {
  return auto_meet(*this, other);
}

std::string DefaultRegisterPredicate::to_string() const {
  return auto_name(*this);
}

bool MaxNQubitsPredicate::verify(const Circuit& circ) const {
  return circ.n_qubits() <= n_qubits_;
}

bool MaxNQubitsPredicate::implies(const Predicate& other) const {
  try {
    const MaxNQubitsPredicate& other_p =
        dynamic_cast<const MaxNQubitsPredicate&>(other);
    return n_qubits_ <= other_p.n_qubits_;
  } catch (const std::bad_cast&) {
    throw IncorrectPredicate(
        "Cannot compare predicates of different subclasses");
  }
}

PredicatePtr MaxNQubitsPredicate::meet(const Predicate& other) const {
  try {
    const MaxNQubitsPredicate& other_p =
        dynamic_cast<const MaxNQubitsPredicate&>(other);
    return std::make_shared<MaxNQubitsPredicate>(
        std::min(n_qubits_, other_p.n_qubits_));
  } catch (const std::bad_cast&) {
    throw IncorrectPredicate(
        "Cannot compare predicates of different subclasses");
  }
}

std::string MaxNQubitsPredicate::to_string() const {
  return auto_name(*this) + "(" + std::to_string(n_qubits_) + ")";
}

bool NoBarriersPredicate::verify(const Circuit& circ) const {
  BGL_FORALL_VERTICES(v, circ.dag, DAG) {
    Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
    if (op->get_type() == OpType::Barrier) return false;
    if (op->get_type() == OpType::CircBox ||
        op->get_type() == OpType::CustomGate) {
      const Box& box = static_cast<const Box&>(*op);
      if (!verify(*box.to_circuit())) return false;
    }
  }
  return true;
}

bool NoBarriersPredicate::implies(const Predicate& other) const {
  return auto_implication(*this, other);
}

PredicatePtr NoBarriersPredicate::meet(const Predicate& other) const {
  return auto_meet(*this, other);
}

std::string NoBarriersPredicate::to_string() const { return auto_name(*this); }

static bool mid_measure_helper(const Command& com, unit_set_t& measured_units) {
  // Rejects gates acting on measured_units
  // Encountering a measurement adds the qubit and bit to measured_units
  // Returns whether or not a mid-circuit measurement is found
  // Applies recursively for CircBoxes
  if (com.get_op_ptr()->get_type() == OpType::Conditional) {
    unit_vector_t all_args = com.get_args();
    const Conditional& cond =
        static_cast<const Conditional&>(*com.get_op_ptr());
    unit_vector_t::iterator arg_it = all_args.begin();
    for (unsigned i = 0; i < cond.get_width(); ++i) {
      if (measured_units.find(*arg_it) != measured_units.end()) return false;
      ++arg_it;
    }
    unit_vector_t new_args = {arg_it, all_args.end()};
    Command new_com = {cond.get_op(), new_args};
    return mid_measure_helper(new_com, measured_units);
  } else if (
      com.get_op_ptr()->get_type() == OpType::CircBox ||
      com.get_op_ptr()->get_type() == OpType::CustomGate) {
    const Box& box = static_cast<const Box&>(*com.get_op_ptr());
    unit_map_t interface;
    unit_set_t inner_set;
    unsigned q_count = 0;
    unsigned b_count = 0;
    for (const UnitID& u : com.get_args()) {
      UnitID inner_unit = (u.type() == UnitType::Qubit)
                              ? static_cast<UnitID>(Qubit(q_count++))
                              : static_cast<UnitID>(Bit(b_count++));
      interface.insert({inner_unit, u});
      if (measured_units.find(u) != measured_units.end()) {
        inner_set.insert(inner_unit);
      }
    }
    for (const Command& c : *box.to_circuit()) {
      if (!mid_measure_helper(c, inner_set)) return false;
    }
    for (const UnitID& u : inner_set) {
      measured_units.insert(interface.at(u));
    }
    return true;
  } else if (com.get_op_ptr()->get_type() == OpType::Measure) {
    std::pair<unit_set_t::iterator, bool> q_inserted =
        measured_units.insert(com.get_args().at(0));
    std::pair<unit_set_t::iterator, bool> c_inserted =
        measured_units.insert(com.get_args().at(1));
    return q_inserted.second && c_inserted.second;
  } else {
    for (const UnitID& a : com.get_args()) {
      if (measured_units.find(a) != measured_units.end()) return false;
    }
    return true;
  }
}

bool NoMidMeasurePredicate::verify(const Circuit& circ) const {
  if (circ.n_bits() == 0) return true;
  unit_set_t measured_units;
  for (const Command& com : circ) {
    if (!mid_measure_helper(com, measured_units)) return false;
  }
  return true;
}

bool NoMidMeasurePredicate::implies(const Predicate& other) const {
  return auto_implication(*this, other);
}

PredicatePtr NoMidMeasurePredicate::meet(const Predicate& other) const {
  return auto_meet(*this, other);
}

std::string NoMidMeasurePredicate::to_string() const {
  return auto_name(*this);
}

bool NoSymbolsPredicate::verify(const Circuit& circ) const {
  return !circ.is_symbolic();
}

bool NoSymbolsPredicate::implies(const Predicate& other) const {
  return auto_implication(*this, other);
}

PredicatePtr NoSymbolsPredicate::meet(const Predicate& other) const {
  return auto_meet(*this, other);
}

std::string NoSymbolsPredicate::to_string() const { return auto_name(*this); }

bool GlobalPhasedXPredicate::verify(const Circuit& circ) const {
  BGL_FORALL_VERTICES(v, circ.dag, DAG) {
    if (circ.get_OpType_from_Vertex(v) == OpType::NPhasedX) {
      if (circ.n_in_edges_of_type(v, EdgeType::Quantum) != circ.n_qubits()) {
        return false;
      }
    }
  }
  return true;
}

bool GlobalPhasedXPredicate::implies(const Predicate& other) const {
  return auto_implication(*this, other);
}

PredicatePtr GlobalPhasedXPredicate::meet(const Predicate& other) const {
  return auto_meet(*this, other);
}

std::string GlobalPhasedXPredicate::to_string() const {
  return auto_name(*this);
}

void to_json(nlohmann::json& j, const PredicatePtr& pred_ptr) {
  if (std::shared_ptr<GateSetPredicate> cast_pred =
          std::dynamic_pointer_cast<GateSetPredicate>(pred_ptr)) {
    j["type"] = "GateSetPredicate";
    j["allowed_types"] = cast_pred->get_allowed_types();
    std::sort(j["allowed_types"].begin(), j["allowed_types"].end());
  } else if (
      std::shared_ptr<NoClassicalControlPredicate> cast_pred =
          std::dynamic_pointer_cast<NoClassicalControlPredicate>(pred_ptr)) {
    j["type"] = "NoClassicalControlPredicate";
  } else if (
      std::shared_ptr<NoFastFeedforwardPredicate> cast_pred =
          std::dynamic_pointer_cast<NoFastFeedforwardPredicate>(pred_ptr)) {
    j["type"] = "NoFastFeedforwardPredicate";
  } else if (
      std::shared_ptr<NoClassicalBitsPredicate> cast_pred =
          std::dynamic_pointer_cast<NoClassicalBitsPredicate>(pred_ptr)) {
    j["type"] = "NoClassicalBitsPredicate";
  } else if (
      std::shared_ptr<NoWireSwapsPredicate> cast_pred =
          std::dynamic_pointer_cast<NoWireSwapsPredicate>(pred_ptr)) {
    j["type"] = "NoWireSwapsPredicate";
  } else if (
      std::shared_ptr<MaxTwoQubitGatesPredicate> cast_pred =
          std::dynamic_pointer_cast<MaxTwoQubitGatesPredicate>(pred_ptr)) {
    j["type"] = "MaxTwoQubitGatesPredicate";
  } else if (
      std::shared_ptr<PlacementPredicate> cast_pred =
          std::dynamic_pointer_cast<PlacementPredicate>(pred_ptr)) {
    j["type"] = "PlacementPredicate";
    j["node_set"] = cast_pred->get_nodes();
  } else if (
      std::shared_ptr<ConnectivityPredicate> cast_pred =
          std::dynamic_pointer_cast<ConnectivityPredicate>(pred_ptr)) {
    j["type"] = "ConnectivityPredicate";
    j["architecture"] = cast_pred->get_arch();
  } else if (
      std::shared_ptr<DirectednessPredicate> cast_pred =
          std::dynamic_pointer_cast<DirectednessPredicate>(pred_ptr)) {
    j["type"] = "DirectednessPredicate";
    j["architecture"] = cast_pred->get_arch();
  } else if (
      std::shared_ptr<CliffordCircuitPredicate> cast_pred =
          std::dynamic_pointer_cast<CliffordCircuitPredicate>(pred_ptr)) {
    j["type"] = "CliffordCircuitPredicate";
  } else if (
      std::shared_ptr<UserDefinedPredicate> cast_pred =
          std::dynamic_pointer_cast<UserDefinedPredicate>(pred_ptr)) {
    j["type"] = "UserDefinedPredicate";
    j["custom"] = "SERIALIZATION OF FUNCTIONS IS NOT YET SUPPORTED";
  } else if (
      std::shared_ptr<DefaultRegisterPredicate> cast_pred =
          std::dynamic_pointer_cast<DefaultRegisterPredicate>(pred_ptr)) {
    j["type"] = "DefaultRegisterPredicate";
  } else if (
      std::shared_ptr<MaxNQubitsPredicate> cast_pred =
          std::dynamic_pointer_cast<MaxNQubitsPredicate>(pred_ptr)) {
    j["type"] = "MaxNQubitsPredicate";
    j["n_qubits"] = cast_pred->get_n_qubits();
  } else if (
      std::shared_ptr<NoBarriersPredicate> cast_pred =
          std::dynamic_pointer_cast<NoBarriersPredicate>(pred_ptr)) {
    j["type"] = "NoBarriersPredicate";
  } else if (
      std::shared_ptr<NoMidMeasurePredicate> cast_pred =
          std::dynamic_pointer_cast<NoMidMeasurePredicate>(pred_ptr)) {
    j["type"] = "NoMidMeasurePredicate";
  } else if (
      std::shared_ptr<NoSymbolsPredicate> cast_pred =
          std::dynamic_pointer_cast<NoSymbolsPredicate>(pred_ptr)) {
    j["type"] = "NoSymbolsPredicate";
  } else if (
      std::shared_ptr<GlobalPhasedXPredicate> cast_pred =
          std::dynamic_pointer_cast<GlobalPhasedXPredicate>(pred_ptr)) {
    j["type"] = "GlobalPhasedXPredicate";
  } else {
    throw JsonError("Cannot serialize PredicatePtr of unknown type.");
  }
}

void from_json(const nlohmann::json& j, PredicatePtr& pred_ptr) {
  std::string classname = j.at("type").get<std::string>();
  if (classname == "GateSetPredicate") {
    OpTypeSet allowed_types = j.at("allowed_types").get<OpTypeSet>();
    pred_ptr = std::make_shared<GateSetPredicate>(allowed_types);
  } else if (classname == "NoClassicalControlPredicate") {
    pred_ptr = std::make_shared<NoClassicalControlPredicate>();
  } else if (classname == "NoFastFeedforwardPredicate") {
    pred_ptr = std::make_shared<NoFastFeedforwardPredicate>();
  } else if (classname == "NoClassicalBitsPredicate") {
    pred_ptr = std::make_shared<NoClassicalBitsPredicate>();
  } else if (classname == "NoWireSwapsPredicate") {
    pred_ptr = std::make_shared<NoWireSwapsPredicate>();
  } else if (classname == "MaxTwoQubitGatesPredicate") {
    pred_ptr = std::make_shared<MaxTwoQubitGatesPredicate>();
  } else if (classname == "PlacementPredicate") {
    node_set_t node_set = j.at("node_set").get<node_set_t>();
    pred_ptr = std::make_shared<PlacementPredicate>(node_set);
  } else if (classname == "ConnectivityPredicate") {
    Architecture arch = j.at("architecture").get<Architecture>();
    pred_ptr = std::make_shared<ConnectivityPredicate>(arch);
  } else if (classname == "DirectednessPredicate") {
    Architecture arch = j.at("architecture").get<Architecture>();
    pred_ptr = std::make_shared<DirectednessPredicate>(arch);
  } else if (classname == "CliffordCircuitPredicate") {
    pred_ptr = std::make_shared<CliffordCircuitPredicate>();
  } else if (classname == "UserDefinedPredicate") {
    throw NotImplemented(
        "Deserialization of UserDefinedPredicates not yet implemented.");
  } else if (classname == "DefaultRegisterPredicate") {
    pred_ptr = std::make_shared<DefaultRegisterPredicate>();
  } else if (classname == "MaxNQubitsPredicate") {
    unsigned n_qubits = j.at("n_qubits").get<unsigned>();
    pred_ptr = std::make_shared<MaxNQubitsPredicate>(n_qubits);
  } else if (classname == "NoBarriersPredicate") {
    pred_ptr = std::make_shared<NoBarriersPredicate>();
  } else if (classname == "NoMidMeasurePredicate") {
    pred_ptr = std::make_shared<NoMidMeasurePredicate>();
  } else if (classname == "NoSymbolsPredicate") {
    pred_ptr = std::make_shared<NoSymbolsPredicate>();
  } else if (classname == "GlobalPhasedXPredicate") {
    pred_ptr = std::make_shared<GlobalPhasedXPredicate>();
  } else {
    throw JsonError("Cannot load PredicatePtr of unknown type.");
  }
}

}  // namespace tket
