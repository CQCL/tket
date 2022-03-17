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

#include <pybind11/functional.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "Circuit/Circuit.hpp"
#include "Mapping/AASLabelling.hpp"
#include "Mapping/AASRoute.hpp"
#include "Mapping/BoxDecomposition.hpp"
#include "Mapping/LexiLabelling.hpp"
#include "Mapping/LexiRoute.hpp"
#include "Mapping/MappingManager.hpp"
#include "Mapping/MultiGateReorder.hpp"
#include "Mapping/RoutingMethodCircuit.hpp"
#include "binder_utils.hpp"

namespace py = pybind11;

namespace tket {
PYBIND11_MODULE(mapping, m) {
  py::class_<RoutingMethod, std::shared_ptr<RoutingMethod>>(
      m, "RoutingMethod",
      "Parent class for RoutingMethod, for inheritance purposes only, not for "
      "usage.")
      .def(py::init<>());

  py::class_<
      RoutingMethodCircuit, std::shared_ptr<RoutingMethodCircuit>,
      RoutingMethod>(
      m, "RoutingMethodCircuit",
      "The RoutingMethod class captures a method for partially mapping logical "
      "subcircuits to physical operations as permitted by some architecture. "
      "Ranked RoutingMethod objects are used by the MappingManager to route "
      "whole circuits.")
      .def(
          py::init<
              const std::function<
                  std::tuple<bool, Circuit, unit_map_t, unit_map_t>(
                      const Circuit&, const ArchitecturePtr&)>&,
              unsigned, unsigned>(),
          "Constructor for a routing method defined by partially routing "
          "subcircuits.\n\n:param route_subcircuit: A function declaration "
          "that given a Circuit and Architecture object, returns a tuple "
          "containing a bool informing MappingManager whether to substitute "
          "the returned circuit into the circuit being routed, "
          "a new modified circuit, the initial logical to physical "
          "qubit mapping of the modified circuit and the permutation of "
          "logical to physical qubit mapping given operations in the "
          "modified circuit\n:param max_size: The maximum number of gates "
          "permitted in a subcircuit\n:param max_depth: The maximum permitted "
          "depth of a subcircuit.",
          py::arg("route_subcircuit"), py::arg("max_size"),
          py::arg("max_depth"));

  py::class_<
      LexiRouteRoutingMethod, std::shared_ptr<LexiRouteRoutingMethod>,
      RoutingMethod>(
      m, "LexiRouteRoutingMethod",
      "Defines a RoutingMethod object for mapping circuits that uses the "
      "Lexicographical Comparison approach outlined in arXiv:1902.08091."
      "Only supports 1-qubit, 2-qubit and barrier gates.")
      .def(
          py::init<unsigned>(),
          "LexiRoute constructor.\n\n:param lookahead: Maximum depth of "
          "lookahead employed when picking SWAP for purpose of logical to "
          "physical mapping.",
          py::arg("lookahead") = 10);

  py::class_<
      AASRouteRoutingMethod, std::shared_ptr<AASRouteRoutingMethod>,
      RoutingMethod>(
      m, "AASRouteRoutingMethod",
      "Defines a RoutingMethod object for mapping circuits that uses the "
      "architecture aware synthesis method implemented in tket.")
      .def(
          py::init<unsigned>(),
          "AASRouteRoutingMethod constructor.\n\n:param aaslookahead: "
          "recursive interation depth of the architecture aware synthesis."
          "method.",
          py::arg("aaslookahead"));

  py::class_<
      AASLabellingMethod, std::shared_ptr<AASLabellingMethod>, RoutingMethod>(
      m, "AASLabellingMethod",
      "Defines a Labeling Method for aas for labelling all unplaced qubits in "
      "a circuit")
      .def(py::init<>(), "AASLabellingMethod constructor.");

  py::class_<
      LexiLabellingMethod, std::shared_ptr<LexiLabellingMethod>, RoutingMethod>(
      m, "LexiLabellingMethod",
      "Defines a RoutingMethod for labelling Qubits that uses the "
      "Lexicographical Comparison approach outlined in arXiv:1902.08091.")
      .def(py::init<>(), "LexiLabellingMethod constructor.");

  py::class_<
      MultiGateReorderRoutingMethod,
      std::shared_ptr<MultiGateReorderRoutingMethod>, RoutingMethod>(
      m, "MultiGateReorderRoutingMethod",
      "Defines a RoutingMethod object for commuting physically permitted "
      "multi-qubit gates to the front of the subcircuit.")
      .def(
          py::init<unsigned, unsigned>(),
          "MultiGateReorderRoutingMethod constructor.\n\n:param max_depth: "
          "Maximum number of layers of gates checked for simultaneous "
          "commutation. "
          "\n:param max_size: Maximum number of gates checked for simultaneous "
          "commutation.",
          py::arg("max_depth") = 10, py::arg("max_size") = 10);

  py::class_<
      BoxDecompositionRoutingMethod,
      std::shared_ptr<BoxDecompositionRoutingMethod>, RoutingMethod>(
      m, "BoxDecompositionRoutingMethod",
      "Defines a RoutingMethod object for decomposing boxes.")
      .def(py::init<>(), "BoxDecompositionRoutingMethod constructor.");

  py::class_<MappingManager>(
      m, "MappingManager",
      "Defined by a pytket Architecture object, maps Circuit logical qubits "
      "to physically permitted Architecture qubits. Mapping is completed by "
      "sequential routing (full or partial) of subcircuits. A custom method "
      "for routing (full or partial) of subcircuits can be defined in Python.")
      .def(
          py::init<const ArchitecturePtr&>(),
          "MappingManager constructor.\n\n:param architecture: pytket "
          "Architecture object.",
          py::arg("architecture"))
      .def(
          "route_circuit", &MappingManager::route_circuit,
          "Maps from given logical circuit to physical circuit. Modification "
          "defined by route_subcircuit, but typically this proceeds by "
          "insertion of SWAP gates that permute logical qubits on physical "
          "qubits.\n\n:param circuit: pytket circuit to be mapped"
          "\n:param routing_methods: Ranked methods to use for routing "
          "subcircuits. In given order, each method is sequentially checked "
          "for viability, with the first viable method being used."
          "\n:param label_isolated_qubits: will not label qubits without gates "
          "or only single qubit gates on them if this is set false",
          py::arg("circuit"), py::arg("routing_methods"),
          py::arg("label_isolated_qubits") = true);
}
}  // namespace tket