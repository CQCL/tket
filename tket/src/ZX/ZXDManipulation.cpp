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

#include "Utils/Assert.hpp"
#include "Utils/GraphHeaders.hpp"
#include "ZX/ZXDiagram.hpp"

namespace tket {

namespace zx {

ZXVert ZXDiagram::add_vertex(ZXGen_ptr op) {
  ZXVertProperties vp{op};
  return boost::add_vertex(vp, *graph);
}

ZXVert ZXDiagram::add_vertex(ZXType type, QuantumType qtype) {
  ZXGen_ptr op = ZXGen::create_gen(type, qtype);
  return add_vertex(op);
}

ZXVert ZXDiagram::add_vertex(
    ZXType type, const Expr& param, QuantumType qtype) {
  ZXGen_ptr op = ZXGen::create_gen(type, param, qtype);
  return add_vertex(op);
}

Wire ZXDiagram::add_wire(
    const ZXVert& va, const ZXVert& vb, const WireProperties& prop) {
  auto [wire, added] = boost::add_edge(va, vb, prop, *graph);
  // add_edge only returns false if the graph cannot support parallel edges, but
  // we have set it up to allow this
  TKET_ASSERT(added);
  return wire;
}

Wire ZXDiagram::add_wire(
    const ZXVert& va, const ZXVert& vb, ZXWireType type, QuantumType qtype,
    std::optional<unsigned> va_port, std::optional<unsigned> vb_port) {
  return add_wire(va, vb, {type, qtype, va_port, vb_port});
}

void ZXDiagram::remove_vertex(const ZXVert& v) {
  // Remove from boundary if `v` is a boundary vertex
  if (is_boundary_type(get_zxtype(v))) {
    ZXVertVec::iterator v_it = std::find(boundary.begin(), boundary.end(), v);
    if (v_it != boundary.end()) boundary.erase(v_it);
  }

  boost::clear_vertex(v, *graph);
  boost::remove_vertex(v, *graph);
}

void ZXDiagram::remove_wire(const Wire& w) { boost::remove_edge(w, *graph); }

bool ZXDiagram::remove_wire(
    const ZXVert& va, const ZXVert& vb, const WireProperties& prop,
    ZXDiagram::WireSearchOption directed) {
  BGL_FORALL_OUTEDGES(va, w, *graph, ZXGraph) {
    if ((target(w) == vb) && (get_wire_info(w) == prop)) {
      remove_wire(w);
      return true;
    }
  }

  // Recheck with reverse wire direction (including ports on the two ends)
  if (directed == WireSearchOption::UNDIRECTED) {
    WireProperties rev_props = prop;
    rev_props.source_port = prop.target_port;
    rev_props.target_port = prop.source_port;
    return remove_wire(vb, va, rev_props, WireSearchOption::DIRECTED);
  }
  return false;
}

void ZXDiagram::symbol_substitution(const symbol_map_t& symbol_map) {
  SymEngine::map_basic_basic sub_map;
  for (const std::pair<const Sym, Expr>& p : symbol_map) {
    ExprPtr s = p.first;
    ExprPtr e = p.second;
    sub_map[s] = e;
  }
  symbol_substitution(sub_map);
}

void ZXDiagram::symbol_substitution(
    const std::map<Sym, double, SymEngine::RCPBasicKeyLess>& symbol_map) {
  SymEngine::map_basic_basic sub_map;
  for (std::pair<Sym, Expr> p : symbol_map) {
    ExprPtr s = p.first;
    ExprPtr e = Expr(p.second);
    sub_map[s] = e;
  }
  symbol_substitution(sub_map);
}

void ZXDiagram::symbol_substitution(const SymEngine::map_basic_basic& sub_map) {
  scalar = scalar.subs(sub_map);
  BGL_FORALL_VERTICES(v, *graph, ZXGraph) {
    ZXGen_ptr new_op = get_vertex_ZXGen_ptr(v)->symbol_substitution(sub_map);
    if (new_op) set_vertex_ZXGen_ptr(v, new_op);
  }
}

SymSet ZXDiagram::free_symbols() const {
  SymSet symbols = expr_free_symbols(get_scalar());
  BGL_FORALL_VERTICES(v, *graph, ZXGraph) {
    const SymSet s = get_vertex_ZXGen_ptr(v)->free_symbols();
    symbols.insert(s.begin(), s.end());
  }
  return symbols;
}

bool ZXDiagram::is_symbolic() const { return !free_symbols().empty(); }

static void check_valid_wire(
    const std::optional<unsigned>& port, QuantumType qtype,
    const std::optional<unsigned>& n_ports, std::vector<bool>& ports_found,
    ZXGen_ptr op) {
  if (port) {
    if (!n_ports)
      throw ZXError("Wire at a named port of an undirected vertex");
    else if (ports_found.at(*port))
      throw ZXError("Multiple wires on the same port of a vertex");
    else
      ports_found.at(*port) = true;
  } else if (n_ports)
    throw ZXError("Wire at an unnamed port of a directed vertex");
  if (!op->valid_edge(port, qtype))
    throw ZXError("QuantumType of wire is incompatible with the given port");
}

void ZXDiagram::check_validity() const {
  std::set<ZXVert> boundary_lookup;
  for (const ZXVert& b : boundary) {
    if (!is_boundary_type(get_zxtype(b)))
      throw ZXError("Non-boundary vertex type in boundary");
    if (!boundary_lookup.insert(b).second)
      throw ZXError("Vertex appears in boundary multiple times");
  }
  BGL_FORALL_VERTICES(v, *graph, ZXGraph) {
    ZXGen_ptr op = get_vertex_ZXGen_ptr(v);
    ZXType type = op->get_type();
    if (is_boundary_type(type)) {
      if (degree(v) != 1)
        throw ZXError("Boundary vertex does not have degree 1");
      if (boundary_lookup.find(v) == boundary_lookup.end())
        throw ZXError("Vertex of boundary type is not in the boundary");
    }
    std::optional<unsigned> n_ports = std::nullopt;
    if (is_directed_type(type)) {
      const ZXDirected& dir = static_cast<const ZXDirected&>(*op);
      n_ports = dir.n_ports();
    }
    std::vector<bool> ports_found(n_ports ? *n_ports : 0, false);
    BGL_FORALL_OUTEDGES(v, w, *graph, ZXGraph) {
      check_valid_wire(source_port(w), get_qtype(w), n_ports, ports_found, op);
    }
    BGL_FORALL_INEDGES(v, w, *graph, ZXGraph) {
      check_valid_wire(target_port(w), get_qtype(w), n_ports, ports_found, op);
    }
    if (n_ports && !std::all_of(
                       ports_found.begin(), ports_found.end(),
                       [](const bool& b) { return b; }))
      throw ZXError("Not all ports of a directed vertex have wires connected");
  }
}

}  // namespace zx

}  // namespace tket
