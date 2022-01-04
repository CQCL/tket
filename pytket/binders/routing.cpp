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
