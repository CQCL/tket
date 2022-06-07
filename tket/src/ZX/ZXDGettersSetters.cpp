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

#include "Utils/GraphHeaders.hpp"
#include "ZX/ZXDiagram.hpp"

namespace tket {

namespace zx {

bool WireProperties::operator==(const WireProperties& other) const {
  return (this->type == other.type) && (this->qtype == other.qtype) &&
         (this->source_port == other.source_port) &&
         (this->target_port == other.target_port);
}

ZXVertVec ZXDiagram::get_boundary(
    std::optional<ZXType> type, std::optional<QuantumType> qtype) const {
  if (!type && !qtype) return boundary;

  ZXVertVec sub_boundary;
  for (const ZXVert& b : boundary) {
    if ((!type || get_zxtype(b) == *type) &&
        (!qtype || *get_qtype(b) == *qtype))
      sub_boundary.push_back(b);
  }
  return sub_boundary;
}

std::unique_ptr<ZXGraph>& ZXDiagram::get_graph() { return graph; }

void ZXDiagram::add_boundary(ZXVert& vert) { boundary.push_back(vert); }
const Expr& ZXDiagram::get_scalar() const { return scalar; }

void ZXDiagram::multiply_scalar(const Expr& sc) { scalar *= sc; }

unsigned ZXDiagram::n_vertices() const { return boost::num_vertices(*graph); }

unsigned ZXDiagram::n_wires() const { return boost::num_edges(*graph); }

unsigned ZXDiagram::count_vertices(ZXType type) const {
  unsigned count = 0;
  BGL_FORALL_VERTICES(v, *graph, ZXGraph) {
    if (get_zxtype(v) == type) ++count;
  }
  return count;
}

unsigned ZXDiagram::count_vertices(ZXType zxtype, QuantumType qtype) const {
  unsigned count = 0;
  BGL_FORALL_VERTICES(v, *graph, ZXGraph) {
    if (get_zxtype(v) == zxtype && get_qtype(v) == qtype) ++count;
  }
  return count;
}

unsigned ZXDiagram::count_wires(ZXWireType type) const {
  unsigned count = 0;
  BGL_FORALL_EDGES(w, *graph, ZXGraph) {
    if (get_wire_type(w) == type) ++count;
  }
  return count;
}

unsigned ZXDiagram::degree(const ZXVert& v) const {
  return boost::out_degree(v, *graph) + boost::in_degree(v, *graph);
}

ZXVertVec ZXDiagram::neighbours(const ZXVert& v) const {
  ZXVertSeqSet neis;
  BGL_FORALL_OUTEDGES(v, w, *graph, ZXGraph) neis.insert(target(w));
  BGL_FORALL_INEDGES(v, w, *graph, ZXGraph) neis.insert(source(w));
  auto& ordered = neis.get<TagSeq>();
  ZXVertVec neighbours(ordered.begin(), ordered.end());
  return neighbours;
}

WireVec ZXDiagram::adj_wires(const ZXVert& v) const {
  WireVec adj;
  BGL_FORALL_OUTEDGES(v, w, *graph, ZXGraph) adj.push_back(w);
  BGL_FORALL_INEDGES(v, w, *graph, ZXGraph) {
    // Only add self loops once (ignore the in_edge wires for loops)
    if (source(w) != v) adj.push_back(w);
  }
  return adj;
}

WireVec ZXDiagram::wires_between(const ZXVert& u, const ZXVert& v) const {
  WireVec wires;
  for (const Wire& w : adj_wires(u)) {
    ZXVert other = other_end(w, u);
    if (other == v) wires.push_back(w);
  }
  return wires;
}

std::optional<Wire> ZXDiagram::wire_between(
    const ZXVert& va, const ZXVert& vb,
    ZXDiagram::WireSearchOption directed) const {
  const auto& [wire, exists] = boost::edge(va, vb, *graph);
  if (exists)
    return wire;
  else if (directed == WireSearchOption::UNDIRECTED)
    return wire_between(vb, va, WireSearchOption::DIRECTED);
  else
    return std::nullopt;
}

Wire ZXDiagram::wire_at_port(
    const ZXVert& v, std::optional<unsigned> port) const {
  Wire w_found;
  unsigned n_found = 0;
  BGL_FORALL_OUTEDGES(v, w, *graph, ZXGraph) {
    if (source_port(w) == port) {
      w_found = w;
      ++n_found;
    }
  }
  BGL_FORALL_INEDGES(v, w, *graph, ZXGraph) {
    if (target_port(w) == port) {
      w_found = w;
      ++n_found;
    }
  }
  if (n_found != 1)
    throw ZXError(
        "Expected only one wire at port, found " + std::to_string(n_found));
  return w_found;
}

ZXGen_ptr ZXDiagram::get_vertex_ZXGen_ptr(const ZXVert& v) const {
  return (*graph)[v].op;
}

std::string ZXDiagram::get_name(const ZXVert& v) const {
  return (*graph)[v].op->get_name();
}

ZXType ZXDiagram::get_zxtype(const ZXVert& v) const {
  return (*graph)[v].op->get_type();
}

std::optional<QuantumType> ZXDiagram::get_qtype(const ZXVert& v) const {
  return (*graph)[v].op->get_qtype();
}

void ZXDiagram::set_vertex_ZXGen_ptr(const ZXVert& v, const ZXGen_ptr& op) {
  (*graph)[v].op = op;
}

WireProperties ZXDiagram::get_wire_info(const Wire& w) const {
  return (*graph)[w];
}

QuantumType ZXDiagram::get_qtype(const Wire& w) const {
  return (*graph)[w].qtype;
}

ZXWireType ZXDiagram::get_wire_type(const Wire& w) const {
  return (*graph)[w].type;
}

ZXVert ZXDiagram::source(const Wire& w) const {
  return boost::source(w, *graph);
}

ZXVert ZXDiagram::target(const Wire& w) const {
  return boost::target(w, *graph);
}

std::optional<unsigned> ZXDiagram::source_port(const Wire& w) const {
  return (*graph)[w].source_port;
}

std::optional<unsigned> ZXDiagram::target_port(const Wire& w) const {
  return (*graph)[w].target_port;
}

ZXVert ZXDiagram::other_end(const Wire& w, const ZXVert& u) const {
  ZXVert s = source(w);
  ZXVert t = target(w);
  if (s == u) {
    return t;
  } else if (t == u) {
    return s;
  } else {
    throw ZXError("In other_end(w, u), u is not adjacent to w.");
  }
}

ZXVert ZXDiagram::vertex_at_end(const Wire& w, WireEnd we) const {
  if (we == WireEnd::Source)
    return source(w);
  else
    return target(w);
}

WireEnd ZXDiagram::end_of(const Wire& w, const ZXVert& u) const {
  if (source(w) == u) {
    return WireEnd::Source;
  } else if (target(w) == u) {
    return WireEnd::Target;
  } else {
    throw ZXError("In end_of(w, u), u is not adjacent to w.");
  }
}

void ZXDiagram::set_wire_info(const Wire& w, const WireProperties& wp) {
  (*graph)[w] = wp;
}

void ZXDiagram::set_wire_qtype(const Wire& w, QuantumType qtype) {
  (*graph)[w].qtype = qtype;
}

void ZXDiagram::set_wire_type(const Wire& w, ZXWireType type) {
  (*graph)[w].type = type;
}

bool ZXDiagram::is_pauli_spider(const ZXVert& v) const {
  ZXGen_ptr op = get_vertex_ZXGen_ptr(v);
  if (!is_spider_type(op->get_type())) return false;
  const PhasedGen& bg = static_cast<const PhasedGen&>(*op);
  std::optional<unsigned> pi2_mult = equiv_Clifford(bg.get_param());
  return (pi2_mult && ((*pi2_mult % 2) == 0));
}

bool ZXDiagram::is_proper_clifford_spider(const ZXVert& v) const {
  ZXGen_ptr op = get_vertex_ZXGen_ptr(v);
  if (!is_spider_type(op->get_type())) return false;
  const PhasedGen& bg = static_cast<const PhasedGen&>(*op);
  std::optional<unsigned> pi2_mult = equiv_Clifford(bg.get_param());
  return (pi2_mult && ((*pi2_mult % 2) == 1));
}

static std::string graphviz_vertex_props(ZXGen_ptr op) {
  std::stringstream ss;

  // Tooltips (rollover text) contains the get_name information
  ss << "tooltip=\"" + op->get_name() + "\" ";

  // Classical nodes are drawn thinner
  if (op->get_qtype() == QuantumType::Classical) ss << "penwidth=1 ";

  // Modify node drawing properties based on this type
  ZXType type = op->get_type();

  switch (type) {
    case ZXType::Input:
    case ZXType::Output:
    case ZXType::Open: {
      ss << "style=\"filled, dashed\" fillcolor=\"white\" shape=circle "
            "label=\""
         << op->get_name() << "\"";
      break;
    }
    case ZXType::ZSpider:
    case ZXType::XSpider: {
      const PhasedGen& bg = static_cast<const PhasedGen&>(*op);
      Expr p = bg.get_param();
      std::string colour = (type == ZXType::ZSpider) ? "green" : "red";
      ss << "fillcolor=\"" << colour << "\" shape=circle label=\"";
      if (!equiv_0(p)) ss << p;
      ss << "\"";
      break;
    }
    case ZXType::Hbox: {
      const PhasedGen& bg = static_cast<const PhasedGen&>(*op);
      Expr p = bg.get_param();
      std::optional<Complex> ev = eval_expr_c(p);
      ss << "fillcolor=\"gold\" shape=square label=\"";
      if (!ev || (std::abs(*ev + 1.) >= EPS)) ss << p;
      ss << "\"";
      break;
    }
    case ZXType::XY:
    case ZXType::XZ:
    case ZXType::YZ:
    case ZXType::PX:
    case ZXType::PY:
    case ZXType::PZ: {
      ss << "shape=circle width=0.1 fixedsize=shape label=\"" << op->get_name()
         << "\\n\\n\"";
      break;
    }
    case ZXType::Triangle: {
      ss << "fillcolor=\"gold\" shape=triangle";
      break;
    }
    case ZXType::ZXBox: {
      ss << "shape=box3d penwidth=2 label=\"Box\"";
      break;
    }
    default: {
      throw ZXError(
          "Trying to render vertex in a ZXDiagram with unknown "
          "ZXType");
    }
  }
  return ss.str();
}

static std::string graphviz_wire_props(const WireProperties& wp) {
  std::stringstream ss;

  // (default assumption is that qtype==Quantum)
  if (wp.qtype == QuantumType::Classical) ss << "penwidth=1 ";

  if (wp.type == ZXWireType::H) ss << "style=dashed color=\"blue\" ";

  // port information:
  // in graphviz, the 'head' of an edge refers to the `Target` end
  // and similarly, the 'tail' of an edge refers to the `Source` end
  if (wp.source_port) ss << " taillabel=\"" << *wp.source_port << "\"";
  if (wp.target_port) ss << " headlabel=\"" << *wp.target_port << "\"";

  return ss.str();
}

std::string ZXDiagram::to_graphviz_str(
    const std::set<ZXVert>& highlights) const {
  std::stringstream out;

  // Construct ZXVert index map (used as vertex IDs by graphviz)
  std::map<ZXVert, unsigned> idm;
  unsigned x = 0;
  BGL_FORALL_VERTICES(v, *graph, ZXGraph) { idm.insert({v, x++}); }

  out << "graph G {\n";

  /**
   * Draw the vertices.
   * By default vertices are assumed to be of type `QuantumType::Quantum`.
   * That is, they will be drawn thick: `penwidth=3` sets this.
   * If a node is classical, we will draw them thinner: `penwidth=1`.
   **/
  out << "node [penwidth=3 style=filled]\n";
  BGL_FORALL_VERTICES(v, *graph, ZXGraph) {
    // Specifying node ID we are defining properties on
    out << idm.at(v);

    // Defining the properties bracket
    out << " [";

    // Define vertex properties based on ZXGen
    out << graphviz_vertex_props(get_vertex_ZXGen_ptr(v));

    // Additional visual properties on node:
    // exterior labels for the node ID information
    out << " xlabel=<<font color=\"grey\">" << idm.at(v) << "</font>>";

    // Highlight the vertices in `highlights` by red border
    auto loc = highlights.find(v);
    if (loc != highlights.end()) out << " color=\"red3\"";

    // Close off the properties list for the vertices
    out << "];\n";
  }

  /**
   * Draw the edges.
   * By default, the edges are assumed to be quantum edges `penwidth=3`,
   * that is, they are drawn 'thick'. Classical edges will be made thinner:
   * that is, `penwidth=1`.
   **/
  out << "edge [penwidth=3]\n";
  BGL_FORALL_EDGES(w, *graph, ZXGraph) {
    // Draw an edge between the Wire's two ends
    ZXVert s = source(w);
    ZXVert t = target(w);
    out << idm.at(s) << " -- " << idm.at(t);

    // Defining the properties bracket
    out << " [";
    out << graphviz_wire_props(get_wire_info(w));
    out << "]\n";
  }

  /**
   * Invisible nodes & connections to force the same ordering for the
   * vertices within the same rank - such that inputs / outputs will be
   * at the same level and with it, have a fixed ordering.
   **/
  out << "rankdir = LR;\n"
         "input_rank [style=invisible];\n"
         "output_rank [style=invisible];\n"
         "input_rank -- output_rank [style=invis];\n";

  out << "{ rank = same\n"
         "input_rank";
  for (const ZXVert& v : get_boundary(ZXType::Input))
    out << " -- " << idm.at(v);
  out << " [style=invis]; }\n";

  out << "{ rank = same\n"
         "output_rank";
  for (const ZXVert& v : get_boundary(ZXType::Output))
    out << " -- " << idm.at(v);
  out << " [style=invis]; }\n";

  out << "}\n";
  return out.str();
}

}  // namespace zx

}  // namespace tket
