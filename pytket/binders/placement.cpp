// Copyright Quantinuum
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

#include "tket/Placement/Placement.hpp"

#include <pybind11/pybind11.h>

#include "binder_json.hpp"
#include "binder_utils.hpp"
#include "typecast.hpp"

namespace py = pybind11;
using json = nlohmann::json;

namespace tket {

void place_with_map(Circuit &circ, std::map<Qubit, Node> &qmap) {
  Architecture arc;
  Placement plobj(arc);
  plobj.place_with_map(circ, qmap);
}

void place_fully_connected(
    Circuit &circ, const FullyConnected &fully_connected) {
  if (circ.n_qubits() > fully_connected.n_nodes()) {
    throw std::logic_error(
        "Circuit has more qubits than the FullyConnected graph has nodes");
  }
  std::map<Qubit, Node> qmap;
  std::vector<Qubit> qubits = circ.all_qubits();
  std::vector<Node> nodes = fully_connected.get_all_nodes_vec();
  TKET_ASSERT(nodes.size() >= qubits.size());
  for (unsigned i = 0; i < qubits.size(); i++) {
    qmap[qubits[i]] = nodes[i];
  }
  place_with_map(circ, qmap);
}

PYBIND11_MODULE(placement, m) {
  py::class_<Placement, std::shared_ptr<Placement>>(
            m, "Placement",
            "The base Placement class, contains methods for getting maps "
            "between Circuit Qubits and Architecture Nodes and for relabelling "
            "Circuit Qubits.")

            .def(py::init<Architecture &>(),
                 "The constructor for a Placement object. The Architecture object "
                 "describes the connectivity between "
                 "qubits.\n\n:param arc: An Architecture object.",
                 py::arg("arc"))
            .def("__repr__",
                 [](const Placement &) { return "<tket::Placement>"; })
            .def("place",
              [](const Placement &placement, Circuit &circ) {
                return placement.place(circ);
              },
              "Relabels Circuit Qubits to Architecture Nodes and 'unplaced'. For "
              "base Placement, all Qubits and labelled 'unplaced'. "
              "\n\n:param circuit: The Circuit being relabelled.",
              py::arg("circuit"))
            .def_static(
              "place_with_map",
              [](Circuit &circ, std::map<Qubit, Node>& qmap) {
                return Placement::place_with_map(circ, qmap);
              },
              "Relabels Circuit Qubits to Architecture Nodes using given map. "
              "\n\n:param circuit: The circuit being relabelled\n:param "
              "qmap: The map from logical to physical qubits to apply.",
              py::arg("circuit"), py::arg("qmap"))
            .def("get_placement_map", &Placement::get_placement_map,
                 "Returns a map from logical to physical qubits that is Architecture "
                 "appropriate for the given Circuit. "
                 "\n\n:param circuit: The circuit a map is designed for."
                 "\n:return: dictionary mapping " CLSOBJS(Qubit) " to "
                 CLSOBJS(Node),
                 py::arg("circuit"))
            .def("get_placement_maps", &Placement::get_all_placement_maps,
                 "Returns a list of maps from logical to physical qubits that "
                 "are Architecture appropriate for the given Circuit. Each map is "
                 "estimated to given a similar SWAP overheard after routing. "
                 "\n\n:param circuit: The circuit the maps are designed for."
                 "\n:param matches: The maximum number of maps returned by the method."
                 "\n:return: list of dictionaries mapping " CLSOBJS(Qubit) " "
                 "to " CLSOBJS(Node),
                 py::arg("circuit"), py::arg("matches")=100)
            .def(
                "to_dict", [](const Placement::Ptr &placement) {
                    return py::object(json(placement)); },
                "Return a JSON serializable dict representation of "
                "the Placement."
                "\n\n:return: dict representing the Placement.")
            .def_static(
                "from_dict", [](const py::dict &dict) {
                    return json(dict).get<Placement::Ptr>(); },
                "Construct Placement instance from JSON serializable "
                "dict representation of the Placement.");

  py::class_<LinePlacement, std::shared_ptr<LinePlacement>, Placement>(
      m, "LinePlacement",
      "The LinePlacement class, contains methods for getting maps "
      "between Circuit Qubits and Architecture Nodes and for relabelling "
      "Circuit Qubits.")
      .def(
          py::init<Architecture &, unsigned, unsigned>(),
          "The constructor for a LinePlacement object. The Architecture "
          "object describes the connectivity "
          "between qubits. In this class, a reduced qubit interaction "
          "subgraph is constructed where each node has maximum outdegree 2 "
          "and does not construct a circle (i.e. lines). "
          "To place the Circuit, a Hamiltonian Path is found in the "
          "Architecture "
          "and this subgraph of lines is assigned along it."
          "\n\n:param arc: An Architecture object."
          "\n:param maximum_line_gates: maximum number of gates in the circuit "
          "considered "
          "when constructing lines for assigning to the graph"
          "\n:param maximum_line_depth: maximum depth of circuit considered "
          "when constructing lines for assigning to the graph",
          py::arg("arc"), py::arg("maximum_line_gates") = 100,
          py::arg("maximum_line_depth") = 100)
      .def("__repr__", [](const Placement &) {
        return "<tket::LinePlacement>";
      });

  py::class_<GraphPlacement, std::shared_ptr<GraphPlacement>, Placement>(
      m, "GraphPlacement",
      "The GraphPlacement class, contains methods for getting maps "
      "between Circuit Qubits and Architecture Nodes and for relabelling "
      "Circuit Qubits.")
      .def(
          py::init<
              const Architecture &, unsigned, unsigned, unsigned, unsigned>(),
          "The constructor for a GraphPlacement object. The Architecture "
          "object describes the connectivity "
          "between qubits. To find a qubit to node assignment, this method "
          "constructs a pattern graph where vertices are Circuit qubits and "
          "edges mean a pair of qubits have an interaction in the circuit, "
          "and then tries to find a weighted subgraph monomorphsim to the "
          "architecture connectivity, or target, graph. Edges in the pattern "
          "graph are weighted by the circuit depth at which the interaction "
          "between a "
          "pair of qubit occurs. "
          "The number of edges added to the pattern graph is effected by the "
          "maximum_pattern_gates and maximum_pattern_depth arguments. "
          "If no subgraph monomorphism can be found, "
          "lower edge weights are removed from the pattern graph, are more "
          "edges "
          "are added to the target graph. Edges added to the pattern graph are "
          "weighted lower to reflect what the distance between the Nodes they  "
          "are added between was on the original target graph. "
          "\n\n:param arc: An Architecture object.\n"
          ":param maximum_matches: The total number of weighted subgraph "
          "monomorphisms that can be found before matches are "
          "returned.\n"
          ":param timeout: Total time in milliseconds before stopping "
          "search for monomorphisms.\n"
          ":param maximum_pattern_gates: The upper bound on the number of "
          "circuit gates used to construct the pattern graph for finding "
          "subgraph monomorphisms.\n"
          ":param maximum_pattern_depth: The upper bound on the circuit depth "
          "gates "
          "are added to the pattern graph to for finding subgraph "
          "monomorphisms.",
          py::arg("arc"), py::arg("maximum_matches") = 1000,
          py::arg("timeout") = 1000, py::arg("maximum_pattern_gates") = 100,
          py::arg("maximum_pattern_depth") = 100)
      .def(
          "__repr__",
          [](const Placement &) { return "<tket::GraphPlacement>"; })
      .def(
          "modify_config",
          [](GraphPlacement & /*pobj*/, py::kwargs /*kwargs*/) {
            PyErr_WarnEx(
                PyExc_DeprecationWarning,
                "GraphPlacement.modify_config no longer changes the parameters "
                "for finding solutions. Please create a new GraphPlacement "
                "object with the changed parameters.",
                1);
            return;
          },
          "Deprecated and no longer modifies parameters for finding solutions. "
          "Please create a new GraphPlacement object instead");

  py::class_<
      NoiseAwarePlacement, std::shared_ptr<NoiseAwarePlacement>, Placement>(
      m, "NoiseAwarePlacement",
      "The NoiseAwarePlacement class, contains methods for getting maps "
      "between Circuit Qubits and Architecture Nodes and for relabelling "
      "Circuit Qubits. It uses gate error rates and readout errors "
      "to find the best placement map.")
      .def(
          py::init<
              Architecture &, avg_node_errors_t, avg_link_errors_t,
              avg_readout_errors_t, unsigned, unsigned, unsigned, unsigned>(),
          "The constructor for a NoiseAwarePlacement object. The Architecture"
          " object describes the connectivity between qubits. "
          "The dictionaries passed as parameters indicate the average "
          "gate errors for single- and two-qubit gates as well as readout"
          "errors.  If no error is given for a given node or pair of nodes,"
          "the fidelity is assumed to be 1."
          "\n\n:param arc: An Architecture object\n"
          ":param node_errors: a dictionary mapping nodes in the "
          "architecture to average single-qubit gate errors\n"
          ":param link_errors: a dictionary mapping pairs of nodes in the "
          "architecture to average two-qubit gate errors\n"
          ":param readout_errors: a dictionary mapping nodes in the "
          "architecture to average measurement readout errors.\n"
          ":param maximum_matches: The total number of weighted subgraph "
          "monomorphisms that can be found before matches are returned.\n"
          ":param timeout: Total time in milliseconds before stopping search "
          "for monomorphisms.\n"
          ":param maximum_pattern_gates: The upper bound on the number of "
          "circuit gates used to construct the pattern graph for finding "
          "subgraph monomorphisms.\n"
          ":param maximum_pattern_depth: The upper bound on the circuit depth "
          "gates are added to the pattern graph to for finding subgraph "
          "monomorphisms.",
          py::arg("arc"), py::arg("node_errors") = py::dict(),
          py::arg("link_errors") = py::dict(),
          py::arg("readout_errors") = py::dict(),
          py::arg("maximum_matches") = 1000, py::arg("timeout") = 1000,
          py::arg("maximum_pattern_gates") = 100,
          py::arg("maximum_pattern_depth") = 100)
      .def(
          "__repr__",
          [](const Placement &) { return "<tket::NoiseAwarePlacement>"; })
      .def(
          "modify_config",
          [](NoiseAwarePlacement & /*pobj*/, py::kwargs /*kwargs*/) {
            PyErr_WarnEx(
                PyExc_DeprecationWarning,
                "NoiseAwarePlacement.modify_config no longer changes the "
                "parameters for finding solutions. Please create a new "
                "NoiseAwarePlacement object with the changed parameters.",
                1);
            return;
          },
          "Deprecated and no longer modifies paramters for finding solutions. "
          "Please create a new NoiseAwarePlacement object instead");

  m.def(
      "place_with_map", &place_with_map,
      "Relabels Circuit Qubits according to a map. If provided map "
      "is partial, remaining Circuit Qubits are left 'unplaced'. "
      "\n\n:param circuit: The Circuit being relabelled. \n:param qmap: "
      "The map from logical to physical qubits to apply.",
      py::arg("circuit"), py::arg("qmap"));

  m.def(
      "place_fully_connected", &place_fully_connected,
      "Relabels all Circuit Qubits to the Node objects of a FullyConnected "
      "object. "
      "\n\n:param circuit: The Circuit being relabelled. \n:param "
      "fully_connected: "
      "FullyConnected object Qubits being relabelled to match.",
      py::arg("circuit"), py::arg("fully_connected"));
}
}  // namespace tket
