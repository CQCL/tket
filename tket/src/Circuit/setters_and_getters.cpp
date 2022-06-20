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

////////////////////////////////////////////////////////
// ALL METHODS TO SET AND GET BASIC CIRCUIT INFORMATION//
////////////////////////////////////////////////////////

#include "Circuit.hpp"
#include "DAGDefs.hpp"
#include "DAGProperties.hpp"
#include "OpType/OpDesc.hpp"
#include "OpType/OpType.hpp"
#include "Ops/OpPtr.hpp"
#include "Utils/Assert.hpp"
#include "Utils/Exceptions.hpp"
#include "Utils/GraphHeaders.hpp"
#include "Utils/TketLog.hpp"

namespace tket {

Circuit::Circuit(const std::string &_name) : Circuit() { name = _name; }

Circuit::Circuit(unsigned n, const std::optional<std::string> _name)
    : Circuit() {
  name = _name;
  add_q_register(q_default_reg(), n);
}

Circuit::Circuit(unsigned n, unsigned m, const std::optional<std::string> _name)
    : Circuit(n, _name) {
  add_c_register(c_default_reg(), m);
}

// Copy constructor.
// Makes no assumptions about the graph
Circuit::Circuit(const Circuit &circ) : Circuit() {
  copy_graph(circ);
  phase = circ.get_phase();
  name = circ.name;
}

// copy assignment. Moves boundary pointers.
Circuit &Circuit::operator=(const Circuit &other)  // (1)
{
  dag = DAG();
  boundary = boundary_t();
  copy_graph(other);
  phase = other.get_phase();
  name = other.name;
  return *this;
}

void Circuit::assert_valid() const {  //
  TKET_ASSERT(is_valid(dag));
}

VertexVec Circuit::all_inputs() const {
  VertexVec ins = q_inputs();
  VertexVec c_ins = c_inputs();
  ins.insert(ins.end(), c_ins.begin(), c_ins.end());
  return ins;
}

VertexVec Circuit::q_inputs() const {
  VertexVec ins;
  for (auto [it, end] = boundary.get<TagType>().equal_range(UnitType::Qubit);
       it != end; it++) {
    ins.push_back(it->in_);
  }
  return ins;
}

VertexVec Circuit::c_inputs() const {
  VertexVec ins;
  for (auto [it, end] = boundary.get<TagType>().equal_range(UnitType::Bit);
       it != end; it++) {
    ins.push_back(it->in_);
  }
  return ins;
}

VertexVec Circuit::all_outputs() const {
  VertexVec outs = q_outputs();
  VertexVec c_outs = c_outputs();
  outs.insert(outs.end(), c_outs.begin(), c_outs.end());
  return outs;
}

VertexVec Circuit::q_outputs() const {
  VertexVec outs;
  for (auto [it, end] = boundary.get<TagType>().equal_range(UnitType::Qubit);
       it != end; it++) {
    outs.push_back(it->out_);
  }
  return outs;
}

VertexVec Circuit::c_outputs() const {
  VertexVec outs;
  for (auto [it, end] = boundary.get<TagType>().equal_range(UnitType::Bit);
       it != end; it++) {
    outs.push_back(it->out_);
  }
  return outs;
}

qubit_vector_t Circuit::all_qubits() const {
  qubit_vector_t all_qbs;
  for (auto [it, end] = boundary.get<TagType>().equal_range(UnitType::Qubit);
       it != end; it++) {
    all_qbs.push_back(Qubit(it->id_));
  }
  std::sort(all_qbs.begin(), all_qbs.end());
  return all_qbs;
}

bit_vector_t Circuit::all_bits() const {
  bit_vector_t all_bs;
  for (auto [it, end] = boundary.get<TagType>().equal_range(UnitType::Bit);
       it != end; it++) {
    all_bs.push_back(Bit(it->id_));
  }
  std::sort(all_bs.begin(), all_bs.end());
  return all_bs;
}

unit_vector_t Circuit::all_units() const {
  unit_vector_t all_us;
  for (const BoundaryElement &el : boundary.get<TagID>()) {
    all_us.push_back(el.id_);
  }
  return all_us;
}

std::map<Bit, unsigned> Circuit::bit_readout() const {
  std::map<Bit, unsigned> res;

  // Order bits to generate indices
  bit_vector_t all_bs = all_bits();
  std::sort(all_bs.begin(), all_bs.end());
  unsigned i = 0;
  for (const Bit &b : all_bs) {
    res.insert({b, i});
    i++;
  }

  return res;
}

std::map<Qubit, unsigned> Circuit::qubit_readout() const {
  std::map<Qubit, unsigned> res;
  std::map<Bit, unsigned> bmap = bit_readout();

  // Find measurement map from qubits to index
  for (auto [it, end] = boundary.get<TagType>().equal_range(UnitType::Qubit);
       it != end; it++) {
    Vertex q_out = it->out_;
    Vertex last_gate = source(get_nth_in_edge(q_out, 0));
    if (get_OpType_from_Vertex(last_gate) == OpType::Measure) {
      Vertex possible_c_out = target(get_nth_out_edge(last_gate, 1));
      if (get_OpType_from_Vertex(possible_c_out) == OpType::ClOutput) {
        Bit b(get_id_from_out(possible_c_out));
        res.insert({Qubit(it->id_), bmap.at(b)});
      }
    }
  }

  return res;
}

std::map<Qubit, Bit> Circuit::qubit_to_bit_map() const {
  std::map<Qubit, Bit> res;
  for (auto [it, end] = boundary.get<TagType>().equal_range(UnitType::Qubit);
       it != end; it++) {
    Vertex q_out = it->out_;
    Vertex last_gate = source(get_nth_in_edge(q_out, 0));
    if (get_OpType_from_Vertex(last_gate) == OpType::Measure) {
      Vertex possible_c_out = target(get_nth_out_edge(last_gate, 1));
      if (get_OpType_from_Vertex(possible_c_out) == OpType::ClOutput) {
        Bit b(get_id_from_out(possible_c_out));
        res.insert({Qubit(it->id_), b});
      }
    }
  }

  return res;
}

bool Circuit::contains_unit(const UnitID &id) const {
  return boundary.get<TagID>().find(id) != boundary.get<TagID>().end();
}

Vertex Circuit::get_in(const UnitID &id) const {
  boundary_t::iterator found = boundary.get<TagID>().find(id);
  if (found == boundary.get<TagID>().end())
    throw CircuitInvalidity(
        "Circuit does not contain unit with id: " + id.repr());
  return found->in_;
}

Vertex Circuit::get_out(const UnitID &id) const {
  boundary_t::iterator found = boundary.get<TagID>().find(id);
  if (found == boundary.get<TagID>().end())
    throw CircuitInvalidity(
        "Circuit does not contain unit with id: " + id.repr());
  return found->out_;
}

UnitID Circuit::get_id_from_in(const Vertex &in) const {
  boundary_t::index<TagIn>::type::iterator found =
      boundary.get<TagIn>().find(in);
  if (found == boundary.get<TagIn>().end())
    throw CircuitInvalidity("Input not found in Circuit");
  return found->id_;
}

UnitID Circuit::get_id_from_out(const Vertex &out) const {
  boundary_t::index<TagOut>::type::iterator found =
      boundary.get<TagOut>().find(out);
  if (found == boundary.get<TagOut>().end())
    throw CircuitInvalidity("Output not found in Circuit");
  return found->id_;
}

opt_reg_info_t Circuit::get_reg_info(std::string reg_name) const {
  boundary_t::index<TagReg>::type::iterator found =
      boundary.get<TagReg>().find(reg_name);
  if (found == boundary.get<TagReg>().end())
    return std::nullopt;
  else
    return found->reg_info();
}

register_t Circuit::get_reg(std::string reg_name) const {
  register_t reg;
  for (auto [it, end] = boundary.get<TagReg>().equal_range(reg_name); it != end;
       it++) {
    if (it->id_.reg_dim() != 1)
      throw CircuitInvalidity("Cannot linearise register " + reg_name);
    reg.insert({it->id_.index()[0], it->id_});
  }
  return reg;
}

// returns the total number of vertices in dag
// just a wrapper for a boost method.
// makes no assumptions about graph structure.
unsigned Circuit::n_vertices() const { return boost::num_vertices(this->dag); }

unsigned Circuit::n_qubits() const {
  return boundary.get<TagType>().count(UnitType::Qubit);
}

unsigned Circuit::n_bits() const {
  return boundary.get<TagType>().count(UnitType::Bit);
}

unsigned Circuit::n_units() const { return boundary.size(); }

unsigned Circuit::n_gates() const { return n_vertices() - 2 * n_units(); }

bool Circuit::is_created(const Qubit &id) const {
  return get_OpType_from_Vertex(get_in(id)) == OpType::Create;
}

bool Circuit::is_discarded(const Qubit &id) const {
  return get_OpType_from_Vertex(get_out(id)) == OpType::Discard;
}

// given a vertex, returns a set of all its successor vertices
// this set can be empty, and no warnings are given if it is
// there are no checks to ensure the vertex exists in the graph
VertexVec Circuit::get_successors(const Vertex &vert) const {
  EdgeVec outs = get_all_out_edges(vert);
  VertexVec children;
  std::unordered_set<Vertex> lookup;
  for (const Edge &e : outs) {
    Vertex succ = target(e);
    if (lookup.find(succ) == lookup.end()) {
      children.push_back(succ);
      lookup.insert(succ);
    }
  }
  return children;
}

VertexVec Circuit::get_successors_of_type(
    const Vertex &vert, EdgeType type) const {
  EdgeVec outs = get_out_edges_of_type(vert, type);
  VertexVec children;
  std::unordered_set<Vertex> lookup;
  for (const Edge &e : outs) {
    Vertex succ = target(e);
    if (lookup.find(succ) == lookup.end()) {
      children.push_back(succ);
      lookup.insert(succ);
    }
  }
  return children;
}

// given a vertex, returns a set of all its predecessor vertices
// this set can be empty, and no warnings are given if it is
// there are no checks to ensure the vertex exists in the graph
VertexVec Circuit::get_predecessors(const Vertex &vert) const {
  EdgeVec ins = get_in_edges(vert);
  VertexVec parents;
  std::unordered_set<Vertex> lookup;
  for (const Edge &e : ins) {
    Vertex pred = source(e);
    if (lookup.find(pred) == lookup.end()) {
      parents.push_back(pred);
      lookup.insert(pred);
    }
  }
  return parents;
}

VertexVec Circuit::get_predecessors_of_type(
    const Vertex &vert, EdgeType type) const {
  EdgeVec ins = get_in_edges_of_type(vert, type);
  VertexVec parents;
  std::unordered_set<Vertex> lookup;
  for (const Edge &e : ins) {
    Vertex pred = source(e);
    if (lookup.find(pred) == lookup.end()) {
      parents.push_back(pred);
      lookup.insert(pred);
    }
  }
  return parents;
}

// Just a wrapper for a boost method
unsigned Circuit::n_edges() const { return boost::num_edges(dag); }

unsigned Circuit::n_edges_of_type(const EdgeType &et) const {
  unsigned count = 0;
  BGL_FORALL_EDGES(e, dag, DAG) {
    if (dag[e].type == et) {
      ++count;
    }
  }
  return count;
}

// return the ports corresponding to an edge
std::pair<port_t, port_t> Circuit::get_ports(const Edge &edge) const {
  return dag[edge].ports;
}

port_t Circuit::get_source_port(const Edge &edge) const {
  return dag[edge].ports.first;
}

port_t Circuit::get_target_port(const Edge &edge) const {
  return dag[edge].ports.second;
}

EdgeType Circuit::get_edgetype(const Edge &edge) const {
  return dag[edge].type;
}

EdgeVec Circuit::get_in_edges(const Vertex &vert) const {
  unsigned n = n_in_edges(vert);
  EdgeVec inedges(n);
  std::vector<bool> port_found(n, false);
  BGL_FORALL_INEDGES(vert, e, dag, DAG) {
    port_t port = get_target_port(e);
    if (port >= n) {
      inedges.resize(port + 1);
      port_found.resize(port + 1, false);
    } else if (port_found[port]) {
      throw CircuitInvalidity("Vertex has multiple inputs on the same port");
    }
    port_found[port] = true;
    inedges[port] = e;
  }
  for (unsigned i = 0; i < n; i++)
    if (!port_found[i]) {
      throw CircuitInvalidity("Input ports on Vertex are non-contiguous");
    }
  return inedges;
}

EdgeVec Circuit::get_in_edges_of_type(const Vertex &vert, EdgeType type) const {
  EdgeVec ins = get_in_edges(vert);
  EdgeVec matching;
  for (const Edge &e : ins) {
    if (get_edgetype(e) == type) {
      matching.push_back(e);
    }
  }
  return matching;
}

std::vector<std::optional<Edge>> Circuit::get_linear_out_edges(
    const Vertex &vert) const {
  unsigned n = n_ports(vert);
  std::vector<std::optional<Edge>> outedges(n);
  BGL_FORALL_OUTEDGES(vert, e, dag, DAG) {
    if (get_edgetype(e) == EdgeType::Boolean) continue;
    port_t port = get_source_port(e);
    if (port >= n) {
      throw CircuitInvalidity("Vertex has an output on an unexpected port");
    }
    if (outedges[port]) {
      throw CircuitInvalidity(
          "Vertex has multiple linear outputs on the same port");
    }
    outedges[port] = e;
  }
  return outedges;
}

EdgeVec Circuit::get_all_out_edges(const Vertex &vert) const {
  std::vector<std::optional<Edge>> lin_outs = get_linear_out_edges(vert);
  std::vector<EdgeVec> b_bundles = get_b_out_bundles(vert);
  EdgeVec outs;
  for (unsigned i = 0; i < lin_outs.size(); ++i) {
    std::optional<Edge> &l_out = lin_outs[i];
    if (l_out) {
      outs.push_back(*l_out);
      outs.insert(outs.end(), b_bundles[i].begin(), b_bundles[i].end());
    }
  }
  return outs;
}

EdgeVec Circuit::get_out_edges_of_type(
    const Vertex &vert, EdgeType type) const {
  if (type == EdgeType::Boolean) {
    std::vector<EdgeVec> bundles = get_b_out_bundles(vert);
    EdgeVec outs;
    for (const EdgeVec &b : bundles) {
      outs.insert(outs.end(), b.begin(), b.end());
    }
    return outs;
  } else {
    std::vector<std::optional<Edge>> outs = get_linear_out_edges(vert);
    EdgeVec matching;
    for (const std::optional<Edge> &e : outs) {
      if (e && get_edgetype(*e) == type) {
        matching.push_back(*e);
      }
    }
    return matching;
  }
}

std::vector<EdgeVec> Circuit::get_b_out_bundles(const Vertex &vert) const {
  unsigned n = n_ports(vert);
  std::vector<EdgeVec> bundles(n);
  BGL_FORALL_OUTEDGES(vert, e, dag, DAG) {
    if (get_edgetype(e) == EdgeType::Boolean) {
      port_t port = get_source_port(e);
      if (port > n) {
        throw CircuitInvalidity("Vertex has an output on an unexpected port");
      }
      bundles.at(port).push_back(e);
    }
  }
  return bundles;
}

std::vector<EdgeVec> Circuit::get_b_in_bundles(const Vertex &vert) const {
  unsigned n = n_ports(vert);
  std::vector<EdgeVec> bundles(n);
  BGL_FORALL_INEDGES(vert, e, dag, DAG) {
    if (get_edgetype(e) == EdgeType::Boolean) {
      port_t port = get_target_port(e);
      if (port > n) {
        throw CircuitInvalidity("Vertex has an output on an unexpected port");
      }
      bundles.at(port).push_back(e);
    }
  }
  return bundles;
}

// n represents the port of the edge at vert_from
// there are no checks to ensure the vertex exists in the graph
// will only return Quantum or Classical edges
Edge Circuit::get_nth_out_edge(const Vertex &vert_from, const port_t &n) const {
  BGL_FORALL_OUTEDGES(vert_from, e, dag, DAG) {
    if (get_edgetype(e) != EdgeType::Boolean && get_source_port(e) == n) {
      return e;
    }
  }
  throw MissingEdge();
}

EdgeVec Circuit::get_nth_b_out_bundle(
    const Vertex &vert_from, const port_t &n) const {
  EdgeVec bundle;
  BGL_FORALL_OUTEDGES(vert_from, e, dag, DAG) {
    if (get_edgetype(e) == EdgeType::Boolean && get_source_port(e) == n) {
      bundle.push_back(e);
    }
  }
  return bundle;
}

// n represents the port of the edge at vert_to
// there are no checks to ensure the vertex exists in the graph
Edge Circuit::get_nth_in_edge(const Vertex &vert_to, const port_t &n) const {
  BGL_FORALL_INEDGES(vert_to, e, dag, DAG) {
    if (get_target_port(e) == n) {
      return e;
    }
  }
  throw MissingEdge();
}

// there are no checks to ensure the vertex exists in the graph
unsigned Circuit::n_in_edges(const Vertex &vert) const {
  return boost::in_degree(vert, dag);
}

unsigned Circuit::n_in_edges_of_type(const Vertex &vert, EdgeType et) const {
  unsigned count = 0;
  BGL_FORALL_INEDGES(vert, e, dag, DAG) {
    if (get_edgetype(e) == et) ++count;
  }
  return count;
}

// there are no checks to ensure the vertex exists in the graph
unsigned Circuit::n_out_edges(const Vertex &vert) const {
  return boost::out_degree(vert, dag);
}

unsigned Circuit::n_out_edges_of_type(const Vertex &vert, EdgeType et) const {
  unsigned count = 0;
  BGL_FORALL_OUTEDGES(vert, e, dag, DAG) {
    if (get_edgetype(e) == et) ++count;
  }
  return count;
}

bool Circuit::is_quantum_node(const Vertex &vert) const {
  return n_in_edges_of_type(vert, EdgeType::Classical) == 0 &&
         n_out_edges_of_type(vert, EdgeType::Classical) == 0;
}

bool Circuit::is_classical_node(const Vertex &vert) const {
  return n_in_edges_of_type(vert, EdgeType::Quantum) == 0 &&
         n_out_edges_of_type(vert, EdgeType::Quantum) == 0;
}

unsigned Circuit::n_ports(const Vertex &vert) const {
  return get_Op_signature_from_Vertex(vert).size();
}

// there are no checks to ensure the vertex exists in the graph
// returns a pointer to the op, try not to dereference and do anything with the
// op
const Op_ptr Circuit::get_Op_ptr_from_Vertex(const Vertex &vert) const {
  return this->dag[vert].op;
}

const std::optional<std::string> &Circuit::get_opgroup_from_Vertex(
    const Vertex &vert) const {
  return this->dag[vert].opgroup;
}

const std::unordered_set<std::string> Circuit::get_opgroups() const {
  std::unordered_set<std::string> opgroups;
  BGL_FORALL_VERTICES(v, dag, DAG) {
    std::optional<std::string> v_opgroup = get_opgroup_from_Vertex(v);
    if (v_opgroup) {
      opgroups.insert(v_opgroup.value());
    }
  }
  return opgroups;
}

void Circuit::set_vertex_Op_ptr(const Vertex &vert, const Op_ptr &op) {
  this->dag[vert].op = op;
}

OpDesc Circuit::get_OpDesc_from_Vertex(const Vertex &vert) const {
  return get_Op_ptr_from_Vertex(vert)->get_desc();
}

// there are no checks to ensure the vertex exists in the graph
OpType Circuit::get_OpType_from_Vertex(const Vertex &vert) const {
  return get_Op_ptr_from_Vertex(vert)->get_type();
}

op_signature_t Circuit::get_Op_signature_from_Vertex(const Vertex &vert) const {
  return get_Op_ptr_from_Vertex(vert)->get_signature();
}

// there are no checks to ensure the vertex exists in the graph
// checks that the edge belongs to the vertex
Edge Circuit::get_next_edge(const Vertex &vert, const Edge &in_edge) const {
  if (target(in_edge) != vert)
    throw CircuitInvalidity(
        "Cannot get next edge: Edge is not an in edge to Vertex");
  port_t order = get_target_port(in_edge);
  return get_nth_out_edge(vert, order);
}

// there are no checks to ensure the vertex exists in the graph
// checks that the edge belongs to the vertex
Edge Circuit::get_last_edge(const Vertex &vert, const Edge &out_edge) const {
  if (source(out_edge) != vert)
    throw CircuitInvalidity(
        "Cannot get last edge: Edge is not an out edge from Vertex");
  port_t order = get_source_port(out_edge);
  return get_nth_in_edge(vert, order);
}

// // given a vertex and corresponding in edge, returns next vertex
// and edge
// there are no checks to ensure the vertex exists in the graph
std::pair<Vertex, Edge> Circuit::get_next_pair(
    const Vertex &current_vertex, const Edge &inedge) const {
  Edge new_edge = get_next_edge(current_vertex, inedge);
  Vertex new_vert = target(new_edge);
  if (new_vert == current_vertex) {
    throw CircuitInvalidity("A qubit path is looping");
  }
  return {new_vert, new_edge};
}

// given a vertex and corresponding out edge, returns previous vertex
// and edge pair
// there are no checks to ensure the vertex exists in the graph
std::pair<Vertex, Edge> Circuit::get_prev_pair(
    const Vertex &current_vertex, const Edge &outedge) const {
  Edge last_edge = get_last_edge(current_vertex, outedge);
  Vertex last_vertex = source(last_edge);
  if (last_vertex == current_vertex) {
    throw CircuitInvalidity("A qubit path is looping");
  }

  return {last_vertex, last_edge};
}

bool Circuit::detect_initial_Op(const Vertex &vertex) const {
  OpType type = get_OpType_from_Vertex(vertex);
  return is_initial_q_type(type) || type == OpType::ClInput;
}

bool Circuit::detect_final_Op(const Vertex &vertex) const {
  OpType type = get_OpType_from_Vertex(vertex);
  return is_final_q_type(type) || type == OpType::ClOutput;
}

bool Circuit::detect_boundary_Op(const Vertex &vertex) const {
  OpType type = get_OpType_from_Vertex(vertex);
  return is_boundary_q_type(type) || is_boundary_c_type(type);
}

bool Circuit::detect_singleq_unitary_op(const Vertex &vert) const {
  const OpDesc desc = get_OpDesc_from_Vertex(vert);
  return desc.is_gate() && desc.is_singleq_unitary();
}

unsigned Circuit::qubit_index(
    const Vertex &vert, PortType port_type, port_t port) const {
  const EdgeVec edges = (port_type == PortType::Source)
                            ? get_out_edges_of_type(vert, EdgeType::Quantum)
                            : get_in_edges_of_type(vert, EdgeType::Quantum);
  unsigned n_edges = edges.size();
  for (unsigned i = 0; i < n_edges; i++) {
    const Edge &e = edges[i];
    port_t p = (port_type == PortType::Source) ? get_source_port(e)
                                               : get_target_port(e);
    if (p == port) return i;
  }
  throw NotValid("Invalid port for vertex");
}

std::optional<Pauli> Circuit::commuting_basis(
    const Vertex &vert, PortType port_type, port_t port) const {
  Op_ptr op = get_Op_ptr_from_Vertex(vert);
  if (op->get_type() == OpType::Conditional) {
    op = static_cast<const Conditional &>(*op).get_op();
  }
  return op->commuting_basis(qubit_index(vert, port_type, port));
}

bool Circuit::commutes_with_basis(
    const Vertex &vert, const std::optional<Pauli> &colour, PortType port_type,
    port_t port) const {
  Op_ptr op = get_Op_ptr_from_Vertex(vert);
  if (op->get_type() == OpType::Conditional) {
    op = static_cast<const Conditional &>(*op).get_op();
  }
  return op->commutes_with_basis(colour, qubit_index(vert, port_type, port));
}

}  // namespace tket
