// Copyright 2019-2021 Cambridge Quantum Computing
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

#include "Routing/Routing.hpp"

#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Architecture/Architecture.hpp"
#include "Circuit/Circuit.hpp"
#include "Transformations/Transform.hpp"
#include "Utils/Json.hpp"
#include "binder_json.hpp"
#include "binder_utils.hpp"
#include "typecast.hpp"

namespace py = pybind11;
using json = nlohmann::json;

namespace tket {

// definitely a better way of doing this ...
void amend_config_from_kwargs(NoiseAwarePlacement &pobj, py::kwargs kwargs) {
  PlacementConfig config_ = pobj.get_config();

  if (kwargs.contains("depth_limit"))
    config_.depth_limit = py::cast<unsigned>(kwargs["depth_limit"]);
  if (kwargs.contains("max_interaction_edges"))
    config_.max_interaction_edges =
        py::cast<unsigned>(kwargs["max_interaction_edges"]);
  if (kwargs.contains("max_matches"))
    config_.vf2_max_matches = py::cast<unsigned>(kwargs["max_matches"]);
  if (kwargs.contains("contraction_ratio"))
    config_.arc_contraction_ratio =
        py::cast<unsigned>(kwargs["contraction_ratio"]);
  if (kwargs.contains("timeout"))
    config_.timeout = py::cast<unsigned>(kwargs["timeout"]);

  pobj.set_config(config_);
}
void amend_config_from_kwargs(GraphPlacement &pobj, py::kwargs kwargs) {
  PlacementConfig config_ = pobj.get_config();

  if (kwargs.contains("depth_limit"))
    config_.depth_limit = py::cast<unsigned>(kwargs["depth_limit"]);
  if (kwargs.contains("max_interaction_edges"))
    config_.max_interaction_edges =
        py::cast<unsigned>(kwargs["max_interaction_edges"]);
  if (kwargs.contains("max_matches"))
    config_.vf2_max_matches = py::cast<unsigned>(kwargs["max_matches"]);
  if (kwargs.contains("contraction_ratio"))
    config_.arc_contraction_ratio =
        py::cast<unsigned>(kwargs["contraction_ratio"]);
  if (kwargs.contains("timeout"))
    config_.timeout = py::cast<unsigned>(kwargs["timeout"]);
  pobj.set_config(config_);
}

void place_with_map(Circuit &circ, qubit_mapping_t &qmap) {
  Architecture arc;
  Placement plobj(arc);
  plobj.place_with_map(circ, qmap);
}

std::pair<Circuit, qubit_mapping_t> route(
    const Circuit &circuit, const Architecture &arc, py::kwargs kwargs) {
  RoutingConfig config = {};
  if (kwargs.contains("swap_lookahead"))
    config.depth_limit = py::cast<unsigned>(kwargs["swap_lookahead"]);
  if (kwargs.contains("bridge_lookahead"))
    config.distrib_limit = py::cast<unsigned>(kwargs["bridge_lookahead"]);
  if (kwargs.contains("bridge_interactions"))
    config.interactions_limit =
        py::cast<unsigned>(kwargs["bridge_interactions"]);
  if (kwargs.contains("bridge_exponent"))
    config.distrib_exponent = py::cast<float>(kwargs["bridge_exponent"]);

  Routing router(circuit, arc);
  Circuit out = router.solve(config).first;
  return {out, router.return_final_map()};
}

PYBIND11_MODULE(routing, m) {
  py::class_<graphs::AbstractGraph<Node>>(
      m, "NodeGraph",
      "Abstract class for describing a device connectivity graph.");

  py::class_<Architecture, graphs::AbstractGraph<Node>, ArchitecturePtr>(
      m, "Architecture",
      "Class describing the connectivity of qubits on a general device.")
      .def(
          py::init([](const std::vector<std::pair<unsigned, unsigned>>
                          &connections) { return Architecture(connections); }),
          "The constructor for an architecture with connectivity "
          "between qubits.\n\n:param connections: A list of pairs "
          "representing qubit indices that can perform two-qubit "
          "operations",
          py::arg("connections"))
      .def(
          py::init<const std::vector<std::pair<Node, Node>> &>(),
          "The constructor for an architecture with connectivity "
          "between qubits.\n\n:param connections: A list of pairs "
          "representing Nodes that can perform two-qubit operations",
          py::arg("connections"))
      .def(
          "__repr__",
          [](const Architecture &arc) {
            return "<tket::Architecture, nodes=" +
                   std::to_string(arc.n_nodes()) + ">";
          })
      .def(
          "get_distance", &Architecture::get_distance,
          "given two nodes in Architecture, "
          "returns distance between them",
          py::arg("node_0"), py::arg("node_1"))
      .def(
          "get_adjacent_nodes", &Architecture::get_neighbour_uids,
          "given a node, returns adjacent nodes in Architecture.",
          py::arg("node"))
      .def_property_readonly(
          "nodes", &Architecture::get_all_nodes_vec,
          "Returns all nodes of architecture as Node objects.")
      .def_property_readonly(
          "coupling", &Architecture::get_all_edges_vec,
          "Returns the coupling map of the Architecture as "
          "UnitIDs. ")
      .def(
          "to_dict", [](const Architecture &arch) { return json(arch); },
          "Return a JSON serializable dict representation of "
          "the Architecture.\n"
          ":return: dict containing nodes and links.")
      .def_static(
          "from_dict", [](const json &j) { return j.get<Architecture>(); },
          "Construct Architecture instance from JSON serializable "
          "dict representation of the Architecture.")
      // as far as Python is concerned, Architectures are immutable
      .def(
          "__deepcopy__",
          [](const Architecture &arc, py::dict = py::dict()) { return arc; })
      .def(
          "__repr__",
          [](const Architecture &arc) {
            return "<tket::Architecture, nodes=" +
                   std::to_string(arc.n_nodes()) + ">";
          })
      .def(py::self == py::self);
  py::class_<
      SquareGrid, Architecture, graphs::AbstractGraph<Node>,
      std::shared_ptr<SquareGrid>>(
      m, "SquareGrid",
      "Architecture class for qubits arranged in a square lattice of "
      "given number of rows and columns. Qubits are arranged with qubits "
      "values increasing first along rows then along columns i.e. for a "
      "3 x 3 grid:\n\n 0 1 2\n\n 3 4 5\n\n 6 7 8")
      .def(
          py::init<const unsigned, const unsigned>(),
          "The constructor for a Square Grid architecture with some "
          "undirected connectivity between qubits.\n\n:param n_rows: "
          "The number of rows in the grid\n:param n_columns: The number "
          "of columns in the grid",
          py::arg("n_rows"), py::arg("n_columns"))
      .def(
          py::init<const unsigned, const unsigned, const unsigned>(),
          "The constructor for  a Square Grid architecture with some "
          "undirected connectivity between qubits.\n\n:param n_rows: "
          "The number of rows in the grid\n:param n_columns: The number "
          "of columns in the grid\n:param n_layers: The number of "
          "layers of grids",
          py::arg("n_rows"), py::arg("n_columns"), py::arg("n_layers"))
      .def(
          "squind_to_qind",
          [](const SquareGrid &self, const unsigned row, const unsigned col) {
            return self.squind_to_qind(row, col);
          },
          "Converts a (row,column) index for a square grid to a "
          "single "
          "qubit index\n\n:param row: The given row index\n:param "
          "column: The given column index\n:return: the "
          "corresponding "
          "global qubit index",
          py::arg("row"), py::arg("column"))
      .def(
          "qind_to_squind", &SquareGrid::qind_to_squind,
          "Converts a single qubit index to a (row,column) index for a "
          "square grid.\n\n:param index: The global qubit "
          "index\n:return: the corresponding grid index as a pair "
          "(row,column)",
          py::arg("index"))
      // as far as Python is concerned, Architectures are immutable
      .def(
          "__deepcopy__",
          [](const SquareGrid &arc, py::dict = py::dict()) { return arc; })
      .def("__repr__", [](const SquareGrid &arc) {
        return "<tket::SquareGrid, rows=" + std::to_string(arc.get_rows()) +
               ", columns=" + std::to_string(arc.get_columns()) +
               ", layers=" + std::to_string(arc.get_layers()) + ">";
      });
  py::class_<
      RingArch, std::shared_ptr<RingArch>, Architecture,
      graphs::AbstractGraph<Node>>(
      m, "RingArch",
      "Architecture class for number of qubits arranged in a ring.")
      .def(
          py::init<const unsigned>(),
          "The constructor for a RingArchitecture with some undirected "
          "connectivity between qubits.\n\n:param number of qubits",
          py::arg("nodes"))
      .def("__repr__", [](const RingArch &arc) {
        return "<tket::RingArch, nodes=" + std::to_string(arc.n_nodes()) + ">";
      });
  py::class_<FullyConnected, graphs::AbstractGraph<Node>>(
      m, "FullyConnected",
      "An architecture with full connectivity between qubits.")
      .def(
          py::init<unsigned>(),
          "Construct a fully-connected architecture."
          "\n\n:param n: number of qubits",
          py::arg("n"))
      .def(
          "__repr__",
          [](const FullyConnected &arc) {
            return "<tket::FullyConnected, nodes=" +
                   std::to_string(arc.n_nodes()) + ">";
          })
      .def(py::self == py::self)
      .def_property_readonly(
          "nodes", &FullyConnected::get_all_nodes_vec,
          "All nodes of the architecture as :py:class:`Node` objects.")
      .def(
          "to_dict", [](const FullyConnected &arch) { return json(arch); },
          "JSON-serializable dict representation of the architecture."
          "\n\n:return: dict containing nodes")
      .def_static(
          "from_dict", [](const json &j) { return j.get<FullyConnected>(); },
          "Construct FullyConnected instance from dict representation.");

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
            .def("place", &Placement::place,
                 "Relabels Circuit Qubits to Architecture Nodes and 'unplaced'. For "
                 "base Placement, all Qubits and labelled 'unplaced'. "
                 "\n\n:param circuit: The Circuit being relabelled.",
                 py::arg("circuit"))
            .def_static(
                    "place_with_map", &Placement::place_with_map,
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
                 "\n:return: list of dictionaries mapping " CLSOBJS(Qubit) " "
                 "to " CLSOBJS(Node),
                 py::arg("circuit"))
            .def(
                "to_dict", [](const PlacementPtr &placement) { return json(placement); },
                "Return a JSON serializable dict representation of "
                "the Placement.\n"
                ":return: dict representing the Placement.")
            .def_static(
                "from_dict", [](const json &j) { return j.get<PlacementPtr>(); },
                "Construct Placement instance from JSON serializable "
                "dict representation of the Placement.");

  py::class_<LinePlacement, std::shared_ptr<LinePlacement>, Placement>(
      m, "LinePlacement",
      "The LinePlacement class, contains methods for getting maps "
      "between Circuit Qubits and Architecture Nodes and for relabelling "
      "Circuit Qubits.")
      .def(
          py::init<Architecture &>(),
          "The constructor for a LinePlacement object. The Architecture "
          "object describes the connectivity "
          "between qubits.\n\n:param arc: An Architecture object.",
          py::arg("arc"))
      .def("__repr__", [](const Placement &) {
        return "<tket::LinePlacement>";
      });

  py::class_<GraphPlacement, std::shared_ptr<GraphPlacement>, Placement>(
      m, "GraphPlacement",
      "The GraphPlacement class, contains methods for getting maps "
      "between Circuit Qubits and Architecture Nodes and for relabelling "
      "Circuit Qubits.")
      .def(
          py::init<Architecture &>(),
          "The constructor for a GraphPlacement object. The Architecture "
          "object describes the connectivity "
          "between qubits.\n\n:param arc: An Architecture object.",
          py::arg("arc"))
      .def(
          "__repr__",
          [](const Placement &) { return "<tket::GraphPlacement>"; })
      .def(
          "modify_config",
          [](GraphPlacement &pobj, py::kwargs kwargs) {
            amend_config_from_kwargs(pobj, kwargs);
          },
          "Overides default Placement parameters to given values. Timeout is "
          "in milliseconds"
          "\n:param \\**kwargs: Parameters for placement: "
          "(int)depth_limit=5, (int)max_interaction_edges=edges in "
          "the "
          "device graph, (int)max_matches=10000, "
          "(int)contraction_ratio=10, (int)timeout=60000.");

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
              avg_readout_errors_t>(),
          "The constructor for a NoiseAwarePlacement object. The Architecture "
          "object describes the connectivity between qubits. "
          "The dictionaries passed as parameters indicate the average "
          "gate errors "
          "for single- and two-qubit gates as well as readout errors. "
          "If no error is given for a given node or pair of nodes, the "
          "fidelity is assumed to be 1."
          "\n\n:param arc: An Architecture object\n"
          ":param node_errors: a dictionary mapping nodes in the "
          "architecture to average single-qubit gate errors\n"
          ":param link_errors: a dictionary mapping pairs of nodes in the "
          "architecture to average two-qubit gate errors\n"
          ":param readout_errors: a dictionary mapping nodes in the "
          "architecture to average measurement readout errors.",
          py::arg("arc"), py::arg("node_errors") = py::dict(),
          py::arg("link_errors") = py::dict(),
          py::arg("readout_errors") = py::dict())
      .def(
          "__repr__",
          [](const Placement &) { return "<tket::NoiseAwarePlacement>"; })
      .def(
          "modify_config",
          [](NoiseAwarePlacement &pobj, py::kwargs kwargs) {
            amend_config_from_kwargs(pobj, kwargs);
          },
          "Overides default Placement parameters to given values. Timeout is "
          "in milliseconds"
          "\n:param \\**kwargs: Parameters for placement: "
          "(int)depth_limit=5, (int)max_interaction_edges=edges in "
          "the "
          "device graph, (int)max_matches=10000, "
          "(int)contraction_ratio=10, (int)timeout=60000.");

  m.def(
      "place_with_map", &place_with_map,
      "Relabels Circuit Qubits according to a map. If provided map "
      "is partial, remaining Circuit Qubits are left 'unplaced'. "
      "\n\n:param circuit: The Circuit being relabelled. \n:param qmap: "
      "The map from logical to physical qubits to apply.",
      py::arg("circuit"), py::arg("qmap"));

  m.def(
      "route",
      [](const Circuit &circuit, const Architecture &arc, py::kwargs kwargs) {
        return route(circuit, arc, kwargs).first;
      },
      "Routes the circuit subject to the connectivity of the input "
      "architecture, given configuration settings."
      "\n\n:param circuit: The circuit to be routed."
      "\n:param architecture: A representation of the qubit connectivity "
      "constraints of the device."
      "\n:param \\**kwargs: Parameters for routing: "
      "(int)swap_lookahead=50, (int)bridge_lookahead=4, "
      "(int)bridge_interactions=2, (float)bridge_exponent=0, "
      "\n:return: the routed :py:class:`Circuit`",
      py::arg("circuit"), py::arg("architecture"));
  m.def(
      "_route_return_map",
      [](const Circuit &circuit, const Architecture &arc, py::kwargs kwargs) {
        return route(circuit, arc, kwargs);
      });
}
}  // namespace tket
