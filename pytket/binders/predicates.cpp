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

#include "Predicates/Predicates.hpp"

#include "Predicates/CompilationUnit.hpp"
#include "Utils/UnitID.hpp"
#include "binder_json.hpp"
#include "typecast.hpp"

namespace py = pybind11;

namespace tket {

static std::map<UnitID, UnitID> unit_bimap_to_map(const unit_bimap_t &bimap) {
  std::map<UnitID, UnitID> res;
  for (auto iter = bimap.left.begin(); iter != bimap.left.end(); ++iter) {
    res.insert({iter->first, iter->second});
  }
  return res;
}

PYBIND11_MODULE(predicates, m) {
  /* Predicates */

  class PyPredicate : public Predicate {
   public:
    using Predicate::Predicate;

    /* Trampolines (need one for each virtual function */
    virtual bool verify(const Circuit &circ) const override {
      PYBIND11_OVERLOAD_PURE(
          bool,      /* Return type */
          Predicate, /* Parent class */
          verify,    /* Name of function in C++ (must match Python name) */
          circ       /* Argument(s) */
      );
    }
    virtual bool implies(const Predicate &other) const override {
      PYBIND11_OVERLOAD_PURE(
          bool,      /* Return type */
          Predicate, /* Parent class */
          implies,   /* Name of function in C++ (must match Python name) */
          other      /* Argument(s) */
      );
    }
    virtual PredicatePtr meet(const Predicate &other) const override {
      PYBIND11_OVERLOAD_PURE(
          PredicatePtr, /* Return type */
          Predicate,    /* Parent class */
          meet,         /* Name of function in C++ (must match Python name) */
          other         /* Argument(s) */
      );
    }
    virtual std::string to_string() const override {
      PYBIND11_OVERLOAD_PURE(
          std::string, /* Return type */
          Predicate,   /* Parent class */
          to_string    /* Name of function in C++ (must match Python name) */
      );
    }
  };

  py::class_<Predicate, PredicatePtr, PyPredicate>(
      m, "Predicate", "A predicate that may be satisfied by a circuit.")
      .def(
          "verify", &Predicate::verify,
          ":return: True if circuit satisfies predicate, else False",
          py::arg("circuit"))
      .def(
          "implies", &Predicate::implies,
          ":return: True if predicate implies another one, else False",
          py::arg("other"))
      .def("__str__", [](const Predicate &) { return "<tket::Predicate>"; })
      .def("__repr__", &Predicate::to_string)
      .def(
          "to_dict",
          [](const PredicatePtr &predicate) {
            return nlohmann::json(predicate);
          },
          "Return a JSON serializable dict representation of "
          "the Predicate.\n\n"
          ":return: dict representation of the Predicate.")
      .def_static(
          "from_dict",
          [](const nlohmann::json &j) { return j.get<PredicatePtr>(); },
          "Construct Predicate instance from JSON serializable "
          "dict representation of the Predicate.");
  py::class_<GateSetPredicate, std::shared_ptr<GateSetPredicate>, Predicate>(
      m, "GateSetPredicate",
      "Predicate asserting that the circuit contains only gates from a "
      "given set.")
      .def(
          py::init<const OpTypeSet &>(), "Construct from a set of gate types.",
          py::arg("allowed_types"))
      .def_property_readonly("gate_set", &GateSetPredicate::get_allowed_types);
  py::class_<
      NoClassicalControlPredicate, std::shared_ptr<NoClassicalControlPredicate>,
      Predicate>(
      m, "NoClassicalControlPredicate",
      "Predicate asserting that a circuit has no classical controls.")
      .def(py::init<>(), "Constructor.");
  py::class_<
      NoFastFeedforwardPredicate, std::shared_ptr<NoFastFeedforwardPredicate>,
      Predicate>(
      m, "NoFastFeedforwardPredicate",
      "Predicate asserting that a circuit has no fast feedforward.")
      .def(py::init<>(), "Constructor.");
  py::class_<
      NoClassicalBitsPredicate, std::shared_ptr<NoClassicalBitsPredicate>,
      Predicate>(
      m, "NoClassicalBitsPredicate",
      "Predicate asserting that a circuit has no classical wires.")
      .def(py::init<>(), "Constructor.");
  py::class_<
      NoWireSwapsPredicate, std::shared_ptr<NoWireSwapsPredicate>, Predicate>(
      m, "NoWireSwapsPredicate",
      "Predicate asserting that a circuit has no wire swaps.")
      .def(py::init<>(), "Constructor.");
  py::class_<
      MaxTwoQubitGatesPredicate, std::shared_ptr<MaxTwoQubitGatesPredicate>,
      Predicate>(
      m, "MaxTwoQubitGatesPredicate",
      "Predicate asserting that a circuit has no gates with more than "
      "two input wires.")
      .def(py::init<>(), "Constructor.");
  py::class_<
      ConnectivityPredicate, std::shared_ptr<ConnectivityPredicate>, Predicate>(
      m, "ConnectivityPredicate",
      "Predicate asserting that a circuit satisfies a given connectivity "
      "graph. The graph is always considered to be undirected.")
      .def(
          py::init<const Architecture &>(),
          "Construct from an :py:class:`Architecture`.",
          py::arg("architecture"));
  py::class_<
      DirectednessPredicate, std::shared_ptr<DirectednessPredicate>, Predicate>(
      m, "DirectednessPredicate",
      "Predicate asserting that a circuit satisfies a given connectivity "
      "graph. The graph is always considered to be directed.")
      .def(
          py::init<const Architecture &>(),
          "Construct from an :py:class:`Architecture`.",
          py::arg("architecture"));
  py::class_<
      CliffordCircuitPredicate, std::shared_ptr<CliffordCircuitPredicate>,
      Predicate>(
      m, "CliffordCircuitPredicate",
      "Predicate asserting that a circuit has only Clifford gates.")
      .def(py::init<>(), "Constructor.");
  py::class_<
      UserDefinedPredicate, std::shared_ptr<UserDefinedPredicate>, Predicate>(
      m, "UserDefinedPredicate", "User-defined predicate.")
      .def(
          py::init<const std::function<bool(const Circuit &)> &>(),
          "Construct from a user-defined function from "
          ":py:class:`Circuit` to `bool`.",
          py::arg("check_function"));
  py::class_<
      DefaultRegisterPredicate, std::shared_ptr<DefaultRegisterPredicate>,
      Predicate>(
      m, "DefaultRegisterPredicate",
      "Predicate asserting that a circuit only uses the default quantum "
      "and classical registers.")
      .def(py::init<>(), "Constructor.");
  py::class_<
      MaxNQubitsPredicate, std::shared_ptr<MaxNQubitsPredicate>, Predicate>(
      m, "MaxNQubitsPredicate",
      "Predicate asserting that a circuit has at most n qubits.")
      .def(py::init<unsigned>(), "Constructor.");
  py::class_<
      PlacementPredicate, std::shared_ptr<PlacementPredicate>, Predicate>(
      m, "PlacementPredicate",
      "Predicate asserting that a circuit has been acted on by some "
      "Placement object.")
      .def(
          py::init<const Architecture &>(),
          "Construct from an :py:class:`Architecture`.",
          py::arg("architecture"))
      .def(
          py::init<const node_set_t &>(), "Construct from a set of Node.",
          py::arg("nodes"));
  py::class_<
      NoBarriersPredicate, std::shared_ptr<NoBarriersPredicate>, Predicate>(
      m, "NoBarriersPredicate",
      "Predicate asserting that a circuit contains no Barrier operations.")
      .def(py::init<>(), "Constructor.");
  py::class_<
      NoMidMeasurePredicate, std::shared_ptr<NoMidMeasurePredicate>, Predicate>(
      m, "NoMidMeasurePredicate",
      "Predicate asserting that all measurements occur at the end of the "
      "circuit.")
      .def(py::init<>(), "Constructor.");
  py::class_<
      NoSymbolsPredicate, std::shared_ptr<NoSymbolsPredicate>, Predicate>(
      m, "NoSymbolsPredicate",
      "Predicate asserting that no gates in the circuit have symbolic "
      "parameters.")
      .def(py::init<>(), "Constructor.");

  /* Compilation units */

  py::class_<CompilationUnit>(
      m, "CompilationUnit",
      "This class comprises a circuit and the predicates that the "
      "circuit is required to satisfy, for example to run on a backend.")
      .def(
          py::init<const Circuit &>(),
          "Construct from a circuit, with no predicates.", py::arg("circuit"))
      .def(
          py::init<const Circuit &, const std::vector<PredicatePtr> &>(),
          "Construct from a circuit and some required predicates.",
          py::arg("circuit"), py::arg("predicates"))
      .def(
          "check_all_predicates", &CompilationUnit::check_all_predicates,
          ":return: True if all predicates are satisfied, else False")
      .def_property_readonly(
          "circuit",
          [](const CompilationUnit &cu) { return Circuit(cu.get_circ_ref()); },
          "Return a copy of the circuit.")
      .def_property_readonly(
          "initial_map",
          [](const CompilationUnit &cu) {
            return unit_bimap_to_map(cu.get_initial_map_ref());
          },
          "Returns the map from the original qubits to the "
          "corresponding qubits at the start of the current circuit.")
      .def_property_readonly(
          "final_map",
          [](const CompilationUnit &cu) {
            return unit_bimap_to_map(cu.get_final_map_ref());
          },
          "Returns the map from the original qubits to their "
          "corresponding qubits at the end of the current circuit.")
      .def(
          "__str__",
          [](const CompilationUnit &) { return "<tket::CompilationUnit>"; })
      .def("__repr__", &CompilationUnit::to_string);
}

}  // namespace tket
