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

#include "Circuit.hpp"

#include <fstream>
#include <numeric>
#include <optional>
#include <set>
#include <string>
#include <utility>

#include "Utils/Expression.hpp"
#include "Utils/GraphHeaders.hpp"
#include "Utils/HelperFunctions.hpp"
#include "Utils/TketLog.hpp"

namespace tket {

////////////////////////////
// Public Circuit Methods //
////////////////////////////

// Out stream of graphviz code describing the circuit from Top-Bottom,
// information apart from qubit path and slices can be seen easily.
// Very useful for debugging and eyeball comparison of circuits.
// out << "\nrankdir=\"LR\"" <--- put this in for Left-Right circuit
void Circuit::to_graphviz(std::ostream &out) const {
  IndexMap im = index_map();

  out << "digraph G {\n";
  out << "{ rank = same\n";
  for (const Vertex &v : all_inputs()) {
    out << im[v] << " ";
  }
  out << "}\n";
  out << "{ rank = same\n";
  for (const Vertex &v : all_outputs()) {
    out << im[v] << " ";
  }
  out << "}\n";

  BGL_FORALL_VERTICES(v, dag, DAG) {
    out << im[v] << " [label = \"" << get_Op_ptr_from_Vertex(v)->get_name()
        << ", " << im[v] << "\"];\n";
  }
  BGL_FORALL_EDGES(e, dag, DAG) {
    Vertex v_so = source(e);
    Vertex v_ta = target(e);
    unsigned v_s = im[v_so];
    unsigned v_t = im[v_ta];
    out << v_s << " -> " << v_t << " [label =  \"" << get_source_port(e) << ", "
        << get_target_port(e) << "\"];\n";
  }
  out << "}";
}

void Circuit::to_graphviz_file(const std::string &filename) const {
  std::ofstream dot_file(filename);
  to_graphviz(dot_file);
}

std::string Circuit::to_graphviz_str() const {
  std::stringstream dot_string;
  to_graphviz(dot_string);
  return dot_string.str();
}

void Circuit::extract_slice_segment(unsigned slice_one, unsigned slice_two) {
  SliceVec slices = get_slices();
  VertexList bin;
  for (unsigned i = 0; i < slice_one - 1; ++i) {
    for (Vertex it : slices[i]) {
      bin.push_back(it);
      remove_vertex(it, GraphRewiring::Yes, VertexDeletion::No);
    }
  }
  for (unsigned i = slice_two; i < slices.size(); ++i) {
    for (Vertex it : slices[i]) {
      bin.push_back(it);
      remove_vertex(it, GraphRewiring::Yes, VertexDeletion::No);
    }
  }
  remove_vertices(bin, GraphRewiring::No, VertexDeletion::Yes);
}

std::vector<Command> Circuit::get_commands() const {
  std::vector<Command> coms;
  for (CommandIterator it = begin(); it != end(); ++it) {
    coms.push_back(*it);
  }
  return coms;
}

void Circuit::index_vertices() /*const*/ {
  VIndex index = boost::get(boost::vertex_index, dag);
  int i = 0;
  BGL_FORALL_VERTICES(v, dag, DAG) { boost::put(index, v, i++); }
}

VertexVec Circuit::vertices_in_order() /*const*/ {
  index_vertices();
  VertexVec vertices;
  boost::topological_sort(dag, std::back_inserter(vertices));
  std::reverse(vertices.begin(), vertices.end());
  return vertices;
}

IndexMap Circuit::index_map() const {
  IndexMap im;
  unsigned i = 0;
  BGL_FORALL_VERTICES(v, dag, DAG) { im[v] = i++; }
  return im;
}

Expr Circuit::get_phase() const {
  std::optional<double> x = eval_expr_mod(phase);
  if (x) {
    return x.value();
  } else
    return phase;
}

void Circuit::add_phase(Expr a) { phase += a; }

void Circuit::symbol_substitution(const symbol_map_t &symbol_map) {
  SymEngine::map_basic_basic sub_map;
  for (const std::pair<const Sym, Expr> &p : symbol_map) {
    ExprPtr s = p.first;
    ExprPtr e = p.second;
    // This is a workaround for a symengine issue: symengine currently has poor
    // handling of symbolic evaluations for atan2. However, this may not catch
    // every such issue, so we should revisit it.
    if (approx_0(e)) {
      sub_map[s] = SymEngine::zero;
    } else {
      sub_map[s] = e;
    }
  }
  symbol_substitution(sub_map);
}

void Circuit::symbol_substitution(
    const std::map<Sym, double, SymEngine::RCPBasicKeyLess> &symbol_map) {
  symbol_map_t s_map;
  for (std::pair<Sym, Expr> p : symbol_map) {
    s_map[p.first] = Expr(p.second);
  }
  symbol_substitution(s_map);
}

void Circuit::symbol_substitution(const SymEngine::map_basic_basic sub_map) {
  BGL_FORALL_VERTICES(v, dag, DAG) {
    Op_ptr new_op = get_Op_ptr_from_Vertex(v)->symbol_substitution(sub_map);
    if (new_op) {
      dag[v] = {new_op};
    }
  }
  phase = phase.subs(sub_map);
}

const SymSet Circuit::free_symbols() const {
  SymSet symbols;
  BGL_FORALL_VERTICES(v, dag, DAG) {
    const SymSet s = get_Op_ptr_from_Vertex(v)->free_symbols();
    symbols.insert(s.begin(), s.end());
  }
  SymSet phase_s = expr_free_symbols(phase);
  symbols.insert(phase_s.begin(), phase_s.end());
  return symbols;
}

bool Circuit::is_symbolic() const { return !free_symbols().empty(); }

// check aspects of circuit for equality, and optionally throw exceptions when
// not met
bool Circuit::circuit_equality(
    const Circuit &other, const std::set<Check> &except,
    bool throw_error) const {
  bool check = true;
  check &= check_iterators_equality(*this, other);
  if (throw_error && !check) {
    throw CircuitInequality(std::string("Circuit operations do not match."));
  }

  if (except.count(Check::Phase) == 0) {
    const Expr thisphase = this->get_phase();
    const Expr othephase = other.get_phase();
    check &= equiv_expr(thisphase, othephase);
    if (throw_error && !check) {
      throw CircuitInequality(
          std::string("Circuit phases do not match: ") +
          ExprPtr(thisphase)->__str__() +
          " != " + ExprPtr(othephase)->__str__());
    }
  }
  if (except.count(Check::Units) == 0) {
    check &= (this->all_qubits() == other.all_qubits());
    if (throw_error && !check) {
      throw CircuitInequality(std::string("Circuit qubits do not match."));
    }

    check &= (this->all_bits() == other.all_bits());
    if (throw_error && !check) {
      throw CircuitInequality(std::string("Circuit bits do not match."));
    }
  }

  if (except.count(Check::ImplicitPermutation) == 0) {
    check &=
        (this->implicit_qubit_permutation() ==
         other.implicit_qubit_permutation());
    if (throw_error && !check) {
      throw CircuitInequality(
          std::string("Circuit implicit permutations do not match."));
    }
  }
  if (except.count(Check::Name) == 0) {
    check &= (this->get_name() == other.get_name());
    if (throw_error && !check) {
      const std::optional<std::string> thisname = this->get_name();
      const std::optional<std::string> othename = other.get_name();
      std::string errormsg = "Circuit names do not match: ";
      errormsg += (thisname ? thisname.value() : "None");
      errormsg += " != ";
      errormsg += (othename ? othename.value() : "None");

      throw CircuitInequality(errormsg);
    }
  }
  return check;
}

// Performs a traversal from the given vertex forwards through the dag, looking
// for something on the target qubit We can prune a path if it reaches the depth
// of the target forward = true returns true if target is in causal future of
// from forward = false checks for causal past (v_to_depth should give reverse
// depth)
// TODO:: rewrite to work with classical boxes
bool Circuit::in_causal_order(
    const Vertex &target, const Vertex &from, bool forward,
    const std::map<Vertex, unsigned> &v_to_depth,
    const std::map<Vertex, unit_set_t> &v_to_units, bool strict) const {
  unsigned target_depth = v_to_depth.at(target);
  if (!strict && from == target) return true;
  if (v_to_depth.at(from) >= target_depth) return false;
  typedef std::function<bool(Vertex, Vertex)> Comp;
  Comp c = [&](Vertex a, Vertex b) {
    unsigned deptha = v_to_depth.at(a);
    unsigned depthb = v_to_depth.at(b);
    if (deptha == depthb) {
      unit_set_t unitsa = v_to_units.at(a);
      unit_set_t unitsb = v_to_units.at(b);
      return unitsa < unitsb;
    }
    return deptha < depthb;
  };
  std::set<Vertex, Comp> to_search(c);
  if (forward) {
    VertexVec succs = get_successors(from);
    for (const Vertex &s : succs) {
      if (v_to_depth.find(s) != v_to_depth.end()) {
        to_search.insert(s);
      }
    }
  } else {
    VertexVec preds = get_predecessors(from);
    to_search.insert(preds.begin(), preds.end());
  }
  unit_set_t lookup_units = v_to_units.at(target);
  while (!to_search.empty()) {
    Vertex v = *to_search.begin();
    to_search.erase(to_search.begin());
    if (v_to_depth.at(v) > target_depth) continue;
    unit_set_t v_units = v_to_units.at(v);
    for (const UnitID &u : lookup_units) {
      if (v_units.find(u) != v_units.end()) {
        return true;
      }
    }
    if (forward) {
      VertexVec succs = get_successors(v);
      for (const Vertex &s : succs) {
        if (v_to_depth.find(s) != v_to_depth.end()) {
          to_search.insert(s);
        }
      }
    } else {
      VertexVec preds = get_predecessors(v);
      to_search.insert(preds.begin(), preds.end());
    }
  }
  return false;
}

}  // namespace tket
