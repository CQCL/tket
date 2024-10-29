// Copyright 2019-2024 Cambridge Quantum Computing
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

#include <optional>
#include <tklog/TketLog.hpp>

#include "binder_json.hpp"
#include "tket/Architecture/Architecture.hpp"
#include "tket/Predicates/CompilerPass.hpp"
#include "tket/Predicates/PassGenerators.hpp"
#include "tket/Predicates/PassLibrary.hpp"
#include "tket/Transformations/Transform.hpp"
#include "typecast.hpp"

namespace tket {

namespace py = pybind11;
using json = nlohmann::json;

Transform lightsabre_transform(
    const Architecture& arc, unsigned seed, unsigned optimisation_level) {
  return Transform([=](Circuit& circ) {
    py::module lightsabre_module =
        py::module::import("pytket.extras.lightsabre");
    py::object lightsabre_transform =
        lightsabre_module.attr("_gen_lightsabre_transformation");
    circ = cast<Circuit>(
        lightsabre_transform(arc, optimisation_level, seed)(circ));
    return true;
  });
}

PassPtr lightsabre_routing(
    const Architecture& arc, unsigned seed, unsigned optimisation_level) {
  // construct pass
  Transform t = lightsabre_transform(arc, seed, optimisation_level) >>
                Transforms::decompose_CX_directed(arc);

  // construct pre-conditions
  PredicatePtr twoqbpred = std::make_shared<MaxTwoQubitGatesPredicate>();
  PredicatePtr n_qubit_pred =
      std::make_shared<MaxNQubitsPredicate>(arc.n_nodes());
  PredicatePtrMap precons{
      CompilationUnit::make_type_pair(twoqbpred),
      CompilationUnit::make_type_pair(n_qubit_pred)};

  // construct post-conditions
  PredicatePtr postcon1 = std::make_shared<ConnectivityPredicate>(arc);
  std::pair<const std::type_index, PredicatePtr> pair1 =
      CompilationUnit::make_type_pair(postcon1);
  PredicatePtr postcon2 = std::make_shared<NoWireSwapsPredicate>();
  PredicatePtrMap s_postcons{pair1, CompilationUnit::make_type_pair(postcon2)};
  std::type_index gateset_ti = typeid(GateSetPredicate);

  PredicateClassGuarantees g_postcons{
      {pair1.first, Guarantee::Clear},
      {gateset_ti, Guarantee::Clear},
      {typeid(MaxTwoQubitGatesPredicate), Guarantee::Clear}};
  PostConditions pc{s_postcons, g_postcons, Guarantee::Preserve};

  // config for json
  nlohmann::json j;
  j["name"] = "LightSABREPass";
  j["architecture"] = arc;
  j["seed"] = seed;
  j["optimisation_level"] = "optimisation_level";

  PassPtr lightsabre_pass = std::make_shared<StandardPass>(precons, t, pc, j);
  PassPtr rebase_pass = gen_auto_rebase_pass(
      {OpType::CX, OpType::SX, OpType::Rz, OpType::X, OpType::TK1}, false);
  return rebase_pass >> lightsabre_pass;
}

PYBIND11_MODULE(extras, m) {
  py::module_::import("pytket._tket.passes");
  m.def(
      "LightSABRE", &lightsabre_routing,
      "Routes circuits to a given architecture using the LightSABRE method "
      "available in Qiskit.",
      py::arg("architecture"), py::arg("seed") = 0,
      py::arg("optimisation_level") = 2);
}

}  // namespace tket