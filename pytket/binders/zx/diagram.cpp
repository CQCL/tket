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

#include "Circuit/Circuit.hpp"
#include "Converters/Converters.hpp"
#include "Utils/GraphHeaders.hpp"
#include "ZX/ZXDiagram.hpp"
#include "typecast.hpp"

namespace py = pybind11;

namespace tket {
namespace zx {

void init_rewrite(py::module& m);

class ZXVertWrapper {
  // ZXVert is a void*
  // pybind11 already has void* bound to PyCapsule
  // PyCapsule equality checks aren't correct on ZXVerts
  // Cannot make void* into an opaque type
  // Define a small wrapper class to safely expose vertices
  // Methods on ZXDiagram may require small lambda wrappers to convert
 private:
  ZXVert v_;

 public:
  ZXVertWrapper() : v_() {}
  ZXVertWrapper(const ZXVert& v) : v_(v) {}
  bool operator==(const ZXVertWrapper& other) const {
    return this->v_ == other.v_;
  }
  std::string to_string() const {
    std::stringstream st;
    st << v_;
    return st.str();
  }
  operator const ZXVert&() const { return v_; }
};

std::pair<ZXDiagram, std::map<UnitID, std::pair<ZXVertWrapper, ZXVertWrapper>>>
wrapped_circuit_to_zx(const Circuit& circ) {
  ZXDiagram zxd;
  boost::bimap<ZXVert, Vertex> bmap;
  std::tie(zxd, bmap) = circuit_to_zx(circ);
  std::map<UnitID, std::pair<ZXVertWrapper, ZXVertWrapper>> ret_map;
  for (auto it = bmap.left.begin(); it != bmap.left.end(); it++) {
    OpType io_type = circ.get_OpType_from_Vertex(it->second);
    if (io_type == OpType::Input || io_type == OpType::ClInput) {
      UnitID uid = circ.get_id_from_in(it->second);
      auto found = ret_map.find(uid);
      if (found == ret_map.end())
        ret_map.insert({uid, {ZXVertWrapper(it->first), ZXVertWrapper()}});
      else
        found->second.first = ZXVertWrapper(it->first);
    } else {
      UnitID uid = circ.get_id_from_out(it->second);
      auto found = ret_map.find(uid);
      if (found == ret_map.end())
        ret_map.insert({uid, {ZXVertWrapper(), ZXVertWrapper(it->first)}});
      else
        found->second.second = ZXVertWrapper(it->first);
    }
  }
  return {std::move(zxd), std::move(ret_map)};
}

std::pair<Circuit, std::map<ZXVertWrapper, UnitID>> wrapped_zx_to_circuit(
    const ZXDiagram& diag) {
  Circuit circ = zx_to_circuit(diag);
  std::map<ZXVertWrapper, UnitID> ret_map;
  ZXVertVec ins = diag.get_boundary(ZXType::Input);
  for (unsigned i = 0; i < ins.size(); ++i) {
    ret_map.insert({ins[i], Qubit(i)});
  }
  ZXVertVec outs = diag.get_boundary(ZXType::Output);
  for (unsigned i = 0; i < outs.size(); ++i) {
    ret_map.insert({outs[i], Qubit(i)});
  }
  return {std::move(circ), std::move(ret_map)};
}

class ZXDiagramPybind {
 public:
  static void init_zxdiagram(py::module& m);
};

void ZXDiagramPybind::init_zxdiagram(py::module& m) {
  py::class_<ZXDiagram, std::shared_ptr<ZXDiagram>>(
      m, "ZXDiagram",
      "Undirected graphs for mixed process ZX diagrams. The boundary is an "
      "ordered list which may mix inputs, outputs, and \"open\" vertices (not "
      "specified to be inputs or outputs). Directed vertices (e.g. Boxes, "
      "Triangles, etc.) have numbered ports to distinguish different incident "
      "edges. The content of each vertex is given by a :py:class:`ZXGen` "
      "generator, describing the :py:class:`ZXType` (e.g. XSpider, Input, "
      "Triangle), the QuantumType for single/doubled versions of typical "
      "generators, and any parameters such as phase. Wires are undirected and "
      "have a :py:class:`ZXWireType` (e.g. Basic, Hadamard) and "
      ":py:class:`QuantumType` (a single wire or a doubled pair for a quantum "
      "system).")
      .def(py::init<>(), "Constructs an empty ZX diagram.")
      .def(
          py::init<unsigned, unsigned, unsigned, unsigned>(),
          "Constructs an empty ZX diagram with a given number of unconnected "
          "boundary vertices.\n\n"
          ":param in: Number of quantum inputs.\n"
          ":param out: Number of quantum outputs.\n"
          ":param classical_in: Number of classical inputs.\n"
          ":param classical_out: Number of classical outputs.",
          py::arg("inputs"), py::arg("outputs"), py::arg("classical_inputs"),
          py::arg("classical_outputs"))
      .def(
          py::init<const ZXDiagram&>(),
          "Constructs a copy of an existing ZX diagram.\n\n"
          ":param other: ZX diagram to copy.",
          py::arg("other"))
      .def(
          "get_boundary",
          [](const ZXDiagram& diag, std::optional<ZXType> type,
             std::optional<QuantumType> qtype) {
            ZXVertVec b = diag.get_boundary(type, qtype);
            return std::vector<ZXVertWrapper>(b.begin(), b.end());
          },
          "Returns handles to boundary vertices in order. Optionally filter by "
          "type of boundary vertex.\n\n"
          ":param type: :py:class:`ZXType` to filter by, from "
          "{:py:meth:`ZXType.Input`, :py:meth:`ZXType.Output`, "
          ":py:meth:`ZXType.Open`, None}. Defaults to None."
          ":param qtype: :py:class:`QuantumType` to filter by, from "
          "{:py:meth:`QuantumType.Quantum`, :py:meth:`QuantumType.Classical`, "
          "None}. Defaults to None.",
          py::arg("type") = std::nullopt, py::arg("qtype") = std::nullopt)
      .def_property_readonly(
          "scalar", &ZXDiagram::get_scalar,
          "Returns the global scalar stored numerically. This may be a "
          "symbolic expression.")
      .def(
          "multiply_scalar", &ZXDiagram::multiply_scalar,
          "Multiplies the global scalar by a numerical (possibly symbolic) "
          "constant.",
          py::arg("scalar"))
      .def_property_readonly(
          "vertices",
          [](const ZXDiagram& diag) {
            std::list<ZXVertWrapper> verts;
            BGL_FORALL_VERTICES(v, *diag.graph, ZXGraph) { verts.push_back(v); }
            return verts;
          },
          "Returns a list of handles to all vertices in the diagram. The order "
          "of vertices may not be semantically relevant.")
      .def_property_readonly(
          "wires",
          [](const ZXDiagram& diag) {
            std::list<Wire> wires;
            BGL_FORALL_EDGES(w, *diag.graph, ZXGraph) { wires.push_back(w); }
            return wires;
          },
          "Returns a list of handles to all wires in the diagram. The order of "
          "wires may not be semantically relevant.")
      .def_property_readonly(
          "n_vertices", &ZXDiagram::n_vertices,
          "Counts the number of vertices in the diagram. Includes boundary "
          "vertices and disconnected vertices.")
      .def_property_readonly(
          "n_wires", &ZXDiagram::n_wires,
          "Counts the number of edges in the diagram.")
      .def(
          "count_vertices",
          (unsigned(ZXDiagram::*)(ZXType) const) & ZXDiagram::count_vertices,
          "Counts the number of vertices of a given :py:class:`ZXType` in the "
          "diagram.",
          py::arg("type"))
      .def(
          "count_wires", &ZXDiagram::count_wires,
          "Counts the number of wired of a given :py:class:`ZXWireType` in the "
          "diagram.",
          py::arg("type"))
      .def(
          "degree",
          [](const ZXDiagram& diag, const ZXVertWrapper& v) {
            return diag.degree(v);
          },
          "Returns the degree of the given vertex.", py::arg("v"))
      .def(
          "neighbours",
          [](const ZXDiagram& diag, const ZXVertWrapper& v) {
            ZXVertVec ns = diag.neighbours(v);
            return std::vector<ZXVertWrapper>(ns.begin(), ns.end());
          },
          "Given a vertex, returns a list of all vertices neighbouring it. "
          "Each neighbour will only appear in the list once regardless of how "
          "many shared edges there are. The order of the neighbour list may "
          "not be semantically relevant.",
          py::arg("v"))
      .def(
          "adj_wires",
          [](const ZXDiagram& diag, const ZXVertWrapper& v) {
            return diag.adj_wires(v);
          },
          "Given a vertex, returns a list of all incident wires. Self-loops "
          "will only appear once in the list. The order of the wire list may "
          "not be semantically relevant.",
          py::arg("v"))
      .def(
          "wires_between",
          [](const ZXDiagram& diag, const ZXVertWrapper& u,
             const ZXVertWrapper& v) { return diag.wires_between(u, v); },
          "Given two vertices, returns a list of all wires between them. The "
          "order of the wire list may not be semantically relevant.",
          py::arg("u"), py::arg("v"))
      .def(
          "wire_between",
          [](const ZXDiagram& diag, const ZXVertWrapper& u,
             const ZXVertWrapper& v) { return diag.wire_between(u, v); },
          "Given two vertices, returns either an arbitrary edge between them "
          "if one exists or None if they are not adjacent.",
          py::arg("u"), py::arg("v"))
      .def(
          "wire_at_port",
          [](const ZXDiagram& diag, const ZXVertWrapper& v, unsigned port) {
            return diag.wire_at_port(v, port);
          },
          "Given a vertex, returns the unique wire at the given port number. "
          "Raises an exception if multiple wires are found at the given port.",
          py::arg("v"), py::arg("port"))
      .def(
          "get_vertex_ZXGen",
          [](const ZXDiagram& diag, const ZXVertWrapper& v) {
            return diag.get_vertex_ZXGen_ptr(v);
          },
          "Returns the content of a given vertex as a :py:class:`ZXGen`.",
          py::arg("v"))
      .def(
          "get_name",
          [](const ZXDiagram& diag, const ZXVertWrapper& v) {
            return diag.get_name(v);
          },
          "Returns the readable string description of a given vertex",
          py::arg("v"))
      .def(
          "get_zxtype",
          [](const ZXDiagram& diag, const ZXVertWrapper& v) {
            return diag.get_zxtype(v);
          },
          "Returns the :py:class:`ZXType` of the given vertex.", py::arg("v"))
      .def(
          "get_qtype",
          [](const ZXDiagram& diag, const ZXVertWrapper& v) {
            return diag.get_qtype(v);
          },
          "Returns the :py:class:`QuantumType` of the given vertex if defined, "
          "None otherwise.",
          py::arg("v"))
      .def(
          "set_vertex_ZXGen",
          [](ZXDiagram& diag, const ZXVertWrapper& v, ZXGen_ptr gen) {
            diag.set_vertex_ZXGen_ptr(v, gen);
          },
          "Updates the content of a given vertex to a particular "
          ":py:class:`ZXGen`.",
          py::arg("v"), py::arg("gen"))
      .def(
          "get_wire_qtype",
          (QuantumType(ZXDiagram::*)(const Wire&) const) & ZXDiagram::get_qtype,
          "Returns the :py:class:`QuantumType` of the given wire.",
          py::arg("w"))
      .def(
          "get_wire_type", &ZXDiagram::get_wire_type,
          "Returns the :py:class:`ZXWireType` of the given wire.", py::arg("w"))
      .def(
          "set_wire_qtype", &ZXDiagram::set_wire_qtype,
          "Updates the :py:class:`QuantumType` of the given wire.",
          py::arg("w"), py::arg("qtype"))
      .def(
          "set_wire_type", &ZXDiagram::set_wire_type,
          "Updates the :py:class:`ZXWireType` of the given wire.", py::arg("w"),
          py::arg("type"))
      .def(
          "get_wire_ends",
          [](const ZXDiagram& diag, const Wire& w) {
            return std::make_tuple(
                std::make_tuple(
                    ZXVertWrapper(diag.source(w)), diag.source_port(w)),
                std::make_tuple(
                    ZXVertWrapper(diag.target(w)), diag.target_port(w)));
          },
          "Returns a tuple ((vertex0, port0), (vertex1, port1)) describing the "
          "two ends of the wire.",
          py::arg("w"))
      .def(
          "other_end",
          [](const ZXDiagram& diag, const Wire& w, const ZXVertWrapper& v) {
            return ZXVertWrapper(diag.other_end(w, v));
          },
          "Given a wire and a vertex at one end of the wire, gives the vertex "
          "at the other end of the wire. This can be used to traverse the "
          "undirected edges of the graph.",
          py::arg("w"), py::arg("v"))
      .def(
          "check_validity", &ZXDiagram::check_validity,
          "Performs a check for the internal validity of the "
          ":py:class:`ZXDiagram` and raises an exception if it is invalid.\n"
          "- Inputs/Outputs must have degree 1 and all exist within the "
          "boundary.\n"
          "- Undirected vertices (those without ports) have no ports on "
          "incident edges.\n"
          "- Directed vertices (those with ports) have exactly one incident "
          "edge at each port.\n"
          "- :py:class:`QuantumType` of wires are compatible with the "
          ":py:class:`QuantumType` s of the ports they attach to.")
      .def(
          "symbol_substitution",
          (void(ZXDiagram::*)(const symbol_map_t&)) &
              ZXDiagram::symbol_substitution,
          "In-place substitution for symbolic expressions; iterated through "
          "each parameterised vertex and performs the substitution. This will "
          "not affect any symbols captured within boxed operations.\n\n"
          ":param symbol_map: A map from SymPy symbols to SymPy expressions.",
          py::arg("symbol_map"))
      .def(
          "symbol_substitution",
          (void(ZXDiagram::*)(
              const std::map<Sym, double, SymEngine::RCPBasicKeyLess>&)) &
              ZXDiagram::symbol_substitution,
          "In-place substitution for symbolic expressions; iterated through "
          "each parameterised vertex and performs the substitution. This will "
          "not affect any symbols captured within boxed operations.\n\n"
          ":param symbol_map: A map from SymPy symbols to floating-point "
          "values.",
          py::arg("symbol_map"))
      .def(
          "free_symbols", &ZXDiagram::free_symbols,
          "Returns the set of symbolic parameters in the diagram.")
      .def(
          "is_symbolic", &ZXDiagram::is_symbolic,
          "Returns True if the diagram contains any free symbols, False "
          "otherwise.")
      .def(
          "add_vertex",
          [](ZXDiagram& diag, ZXGen_ptr gen) {
            return ZXVertWrapper(diag.add_vertex(gen));
          },
          "Adds a new vertex to the diagram for an arbitrary "
          ":py:class:`ZXGen`.\n\n"
          ":param gen: The :py:class:`ZXGen` for the new vertex.\n"
          ":return: The handle to the new vertex.",
          py::arg("gen"))
      .def(
          "add_vertex",
          [](ZXDiagram& diag, ZXType type, QuantumType qtype) {
            return ZXVertWrapper(diag.add_vertex(type, qtype));
          },
          "Adds a new vertex to the diagram for an unparameterised, doubleable "
          "generator type.\n\n"
          ":param type: The :py:class:`ZXType` for the new vertex.\n"
          ":param qtype: The :py:class:`QuantumType` for the new vertex. "
          "Defaults to Quantum.\n"
          ":return: The handle to the new vertex.",
          py::arg("type"), py::arg("qtype") = QuantumType::Quantum)
      .def(
          "add_vertex",
          [](ZXDiagram& diag, ZXType type, const Expr& param,
             QuantumType qtype) {
            return ZXVertWrapper(diag.add_vertex(type, param, qtype));
          },
          "Adds a new vertex to the diagram for a parameterised, doubleable "
          "generator type.\n\n"
          ":param type: The :py:class:`ZXType` for the new vertex.\n"
          ":param param: The parameter for the new vertex.\n"
          ":param qtype: The :py:class:`QuantumType` for the new vertex. "
          "Defaults to Quantum.\n"
          ":return: The handle to the new vertex.",
          py::arg("type"), py::arg("param"),
          py::arg("qtype") = QuantumType::Quantum)
      .def(
          "add_zxbox",
          [](ZXDiagram& diag, const ZXDiagram& inner) {
            ZXGen_ptr box = std::make_shared<const ZXBox>(inner);
            return ZXVertWrapper(diag.add_vertex(box));
          },
          "Adds a new vertex to the diagram for a box with some inner "
          "implementation.\n\n"
          ":param inner: The :py:class:`ZXDiagram` to internalise inside the "
          "box. The current state is copied by value.\n"
          ":return: The handle to the new vertex.",
          py::arg("inner"))
      .def(
          "add_wire",
          [](ZXDiagram& diag, const ZXVertWrapper& u, const ZXVertWrapper& v,
             ZXWireType type, QuantumType qtype, std::optional<unsigned> u_port,
             std::optional<unsigned> v_port) {
            return diag.add_wire(u, v, type, qtype, u_port, v_port);
          },
          "Adds a new wire to the diagram between the given vertices.\n\n"
          ":param u: Handle to the first vertex.\n"
          ":param v: Handle to the other vertex.\n"
          ":param type: :py:class:`ZXWireType` for the wire. Defaults to "
          "Basic.\n"
          ":param qtype: :py:class:`QuantumType` for the wire. Defaults to "
          "Quantum.\n"
          ":param u_port: Port on vertex u to connect to. Defaults to None.\n"
          ":param v_port: Port on vertex v to connect to. Defaults to None.\n"
          ":return: The handle to the new wire.",
          py::arg("u"), py::arg("v"), py::arg("type") = ZXWireType::Basic,
          py::arg("qtype") = QuantumType::Quantum,
          py::arg("u_port") = std::nullopt, py::arg("v_port") = std::nullopt)
      .def(
          "remove_vertex",
          [](ZXDiagram& diag, const ZXVertWrapper& v) {
            return diag.remove_vertex(v);
          },
          "Removes the given vertex and all incident wires from the diagram. "
          "If the vertex is in the boundary, it is removed from the boundary.",
          py::arg("v"))
      .def(
          "remove_wire",
          (void(ZXDiagram::*)(const Wire&)) & ZXDiagram::remove_wire,
          "Removes the given wire from the diagram.", py::arg("w"))
      .def(
          "to_circuit", &wrapped_zx_to_circuit,
          "Extracts a unitary diagram in MBQC form as a Circuit following the "
          "routine by Backens et al. (\"There and back again: A circuit "
          "extraction tale\").\n\n"
          ":return: A pair of the generated :py:class:`Circuit`, and a map "
          "from each boundary vertex in the :py:class:`ZXDiagram` to its "
          "corresponding :py:class:`UnitID` in the :py:class:`Circuit`.")
      .def(
          "to_doubled_diagram", &ZXDiagram::to_doubled_diagram,
          "Expands any quantum vertices into pairs of classical vertices "
          "according to the doubling construction for CPM. New boundary "
          "vertices are ordered lexicographically by (b, c):\n"
          "- b boundary index in the original diagram\n"
          "- c conjugate identifier\n"
          "  + quantum boundaries are mapped to a pair with original and "
          "conjugated phases\n"
          "  + unconjugated copies are listed first\n",
          "  + classical boundaries only have the unconjugated version")
      .def(
          "to_graphviz_str",
          [](ZXDiagram& diag) { return diag.to_graphviz_str(); },
          "Returns a graphviz source string");
}

PYBIND11_MODULE(zx, m) {
  py::enum_<ZXType>(
      m, "ZXType",
      "Enum for available types of generators in :py:class:`ZXDiagram` s.")
      .value(
          "Input", ZXType::Input,
          "An input boundary vertex. Can either be Quantum or Classical. Must "
          "have degree 1. No ports.")
      .value(
          "Output", ZXType::Output,
          "An output boundary vertex. Can either be Quantum or Classical. Must "
          "have degree 1. No ports.")
      .value(
          "Open", ZXType::Open,
          "A boundary vertex that has not yet been specified as input or "
          "output. Can either be Quantum or Classical. Must have degree 1. No "
          "ports.")
      .value(
          "ZSpider", ZXType::ZSpider,
          "A Z (green) spider. Parameterised by a single phase in half-turns. "
          "Can either be Quantum or Classical - Quantum spiders can only have "
          "Quantum wires, Quantum wires on Classical spiders act as two wires. "
          "Can have arbitrary degree. No ports.")
      .value(
          "XSpider", ZXType::XSpider,
          "An X (red) spider. Parameterised by a single phase in half-turns. "
          "Can either be Quantum or Classical - Quantum spiders can only have "
          "Quantum wires, Quantum wires on Classical spiders act as two wires. "
          "Can have arbitrary degree. No ports.")
      .value(
          "Hbox", ZXType::Hbox,
          "A Hadamard box for ZH diagrams. Parameterised by a single complex "
          "value. Can either be Quantum or Classical - Quantum spiders can "
          "only have Quantum wires, Quantum wires on Classical spiders act as "
          "two wires. Can have arbitrary degree. No ports.")
      .value(
          "XY", ZXType::XY,
          "A (postselected) XY qubit in MBQC. Corresponds to a Z spider with "
          "negative phase.")
      .value(
          "XZ", ZXType::XZ,
          "A (postselected) XZ qubit in MBQC. Corresponds to a 0.5-phase "
          "(n+1)-ary Z spider connected to a phaseful 1-ary X spider.")
      .value(
          "YZ", ZXType::YZ,
          "A (postselected) YZ qubit in MBQC. Corresponds to a 0-phase "
          "(n+1)-ary Z spider connected to a phaseful 1-ary X spider.")
      .value(
          "PX", ZXType::PX,
          "A (postselected) Pauli X qubit in MBQC. Corresponds to a Z spider "
          "with phase either 0 (param=False) or 1 (param=True).")
      .value(
          "PY", ZXType::PY,
          "A (postselected) Pauli Y qubit in MBQC. Corresponds to a Z spider "
          "with phase either -0.5 (param=False) or +0.5 (param=True).")
      .value(
          "PZ", ZXType::PZ,
          "A (postselected) Pauli Z qubit in MBQC. Corresponds to a 0-phase "
          "(n+1)-ary Z spider connected to a 1-ary X spider with phase either "
          "0 (param=False) or 1 (param=True).")
      .value(
          "Triangle", ZXType::Triangle,
          "A Triangle operator, [[1, 1], [0, 1]]. Can either be Quantum or "
          "Classical, only admitting wires of the same type. Port 0 for the "
          "base of the triangle (input), port 1 for the tip (output).")
      .value(
          "ZXBox", ZXType::ZXBox,
          "A box encapsulating another :py:class:`ZXDiagram`. Inherits ports "
          "from the boundary of the internal diagram, with port numbers "
          "matching the boundary order and :py:class:`QuantumType` admitted at "
          "each port matching that of the boundary vertex.");
  py::enum_<ZXWireType>(
      m, "ZXWireType",
      "Enum for available types of wires in :py:class:`ZXDiagram` s.")
      .value("Basic", ZXWireType::Basic, "A basic identity wire.")
      .value("H", ZXWireType::H, "A Hadamard edge.");
  py::enum_<QuantumType>(
      m, "QuantumType",
      "Enum for specifying quantumness of vertices, ports, and wires in "
      ":py:class:`ZXDiagram` s for mixed quantum-classical processes.")
      .value(
          "Quantum", QuantumType::Quantum,
          "Quantum components of diagrams, represented in the framework of "
          "completely-positive maps by two parallel copies of a system related "
          "by conjugation.")
      .value(
          "Classical", QuantumType::Classical,
          "Classical components of diagrams, represented in the framework of "
          "completely-positive maps by a single self-conjugate system.");
  py::class_<ZXVertWrapper>(
      m, "ZXVert",
      "A handle to a vertex in a :py:class:`ZXDiagram`. Each instance is "
      "specific to a given :py:class:`ZXDiagram` instance and can be "
      "invalidated by rewrites. Exceptions or errors may occur if calling "
      "functions on a :py:class:`ZXVert` that is not present in the given "
      ":py:class:`ZXDiagram`.")
      .def("__repr__", &ZXVertWrapper::to_string)
      .def("__eq__", &ZXVertWrapper::operator==)
      .def("__hash__", [](const ZXVertWrapper& v) {
        return py::hash(py::str(v.to_string()));
      });
  py::class_<Wire>(
      m, "ZXWire",
      "A handle to a wire in a :py:class:`ZXDiagram`. Each instance is "
      "specific to a given :py:class:`ZXDiagram` instance and can be "
      "invalidated by rewrites. Exceptions or errors may occur if calling "
      "functions on a :py:class:`ZXWire` that is not present in the given "
      ":py:class:`ZXDiagram`.")
      .def("__eq__", [](const Wire& a, const Wire& b) { return a == b; })
      .def("__hash__", [](const Wire& w) {
        std::stringstream st;
        st << w;
        return py::hash(py::str(st.str()));
      });
  py::class_<ZXGen, std::shared_ptr<ZXGen>>(
      m, "ZXGen",
      "Encapsulates the information about the generator depicted by a given "
      "vertex in a :py:class:`ZXDiagram`.")
      .def_static(
          "create",
          [](ZXType type, QuantumType qtype) {
            return ZXGen::create_gen(type, qtype);
          },
          "Create a boundary type generator.", py::arg("type"),
          py::arg("qtype") = QuantumType::Quantum)
      .def_static(
          "create",
          [](ZXType type, const Expr& param, QuantumType qtype) {
            return ZXGen::create_gen(type, param, qtype);
          },
          "Create a boundary type generator.", py::arg("type"),
          py::arg("param"), py::arg("qtype") = QuantumType::Quantum)
      .def_property_readonly("type", &ZXGen::get_type, "The type of generator.")
      .def_property_readonly(
          "qtype", &ZXGen::get_qtype,
          "The :py:class:`QuantumType` of the generator (if applicable).")
      .def("__eq__", &ZXGen::operator==)
      .def("__repr__", [](const ZXGen& gen) { return gen.get_name(); });
  py::class_<PhasedGen, std::shared_ptr<PhasedGen>, ZXGen>(
      m, "PhasedGen",
      "Specialisation of :py:class:`ZXGen` for arbitrary-arity, symmetric "
      "generators with a single continuous parameter.")
      .def(
          py::init<ZXType, const Expr&, QuantumType>(),
          "Construct from a ZX type, parameter and quantum type.",
          py::arg("zxtype"), py::arg("param"), py::arg("qtype"))
      .def_property_readonly(
          "param", &PhasedGen::get_param, "The parameter of the generator.");
  py::class_<CliffordGen, std::shared_ptr<CliffordGen>, ZXGen>(
      m, "CliffordGen",
      "Specialisation of :py:class:`ZXGen` for arbitrary-arity, symmetric "
      "Clifford generators with a single boolean parameter.")
      .def(
          py::init<ZXType, bool, QuantumType>(),
          "Construct from a ZX type, parameter and quantum type.",
          py::arg("zxtype"), py::arg("param"), py::arg("qtype"))
      .def_property_readonly(
          "param", &CliffordGen::get_param, "The parameter of the generator.");
  py::class_<DirectedGen, std::shared_ptr<DirectedGen>, ZXGen>(
      m, "DirectedGen",
      "Specialisation of :py:class:`ZXGen` for asymmetric ZX generators which "
      "can be doubled to form a Quantum variant. Asymmetric effects handled by "
      "ports to distinguish operands.")
      .def(
          py::init<ZXType, QuantumType>(),
          "Construct from a ZX type and quantum type.", py::arg("zxtype"),
          py::arg("qtype"))
      .def_property_readonly(
          "n_ports", &DirectedGen::n_ports,
          "The number of ports on the generator.")
      .def_property_readonly(
          "signature", &DirectedGen::get_signature,
          "A list of :py:class:`QuantumType` s indicating the expected "
          ":py:class:`QuantumType` at each port.");
  ZXDiagramPybind::init_zxdiagram(m);
  py::class_<ZXBox, std::shared_ptr<ZXBox>, ZXGen>(
      m, "ZXBox",
      "Specialisation of :py:class:`ZXGen` for encapsulations of some other ZX "
      "diagrams. In general, arbitrary diagrams may be asymmetric tensors with "
      "both Quantum and Classical boundaries, so ports are used to distinguish "
      "each boundary.")
      .def(
          py::init<const ZXDiagram&>(), "Construct from a ZX diagram.",
          py::arg("zxdiag"))
      .def_property_readonly(
          "n_ports", &ZXBox::n_ports, "The number of ports on the generator.")
      .def_property_readonly(
          "signature", &ZXBox::get_signature,
          "A list of :py:class:`QuantumType` s indicating the expected "
          ":py:class:`QuantumType` at each port.")
      .def_property_readonly(
          "diagram", &ZXBox::get_diagram,
          "The internal diagram represented by the box.");
  init_rewrite(m);
  m.def(
      "circuit_to_zx", &wrapped_circuit_to_zx,
      "Construct a ZX diagram from a circuit. Return the ZX diagram and a map "
      "Between the ZX boundary vertices and the resource UIDs of the circuit.");
}

}  // namespace zx
}  // namespace tket
