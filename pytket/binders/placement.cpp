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

#include "Placement/Placement.hpp"

#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Utils/Json.hpp"
#include "binder_json.hpp"
#include "binder_utils.hpp"
#include "typecast.hpp"

namespace py = pybind11;
using json = nlohmann::json;

namespace tket {

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

void place_fully_connected(
    Circuit &circ, const FullyConnected &fully_connected) {
  if (circ.n_qubits() > fully_connected.n_nodes()) {
    throw std::logic_error(
        "Circuit has more qubits than the FullyConnected graph has nodes");
  }
  qubit_mapping_t qmap;
  unsigned index = 0;
  for (const Qubit &q : circ.all_qubits()) {
    qmap[q] = Node("fcNode", index);
    index++;
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
              [](Circuit &circ, qubit_mapping_t& qmap) {
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
      "place_fully_connected", &place_fully_connected,
      "Relabels all Circuit Qubits to the Node objects of a FullyConnected "
      "object. "
      "\n\n:param circuit: The Circuit being relabelled. \n:param "
      "fully_connected: "
      "FullyConnected object Qubits being relabelled to match.",
      py::arg("circuit"), py::arg("fully_connected"));
}
}  // namespace tket
