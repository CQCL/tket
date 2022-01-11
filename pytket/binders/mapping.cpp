#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Mapping/LexiRoute.hpp"
#include "Mapping/MappingManager.hpp"
#include "Mapping/RoutingMethodCircuit.hpp"

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
      "The RoutingMethod class captures a method for partially mapping logical"
      "subcircuits to physical operations as permitted by some architecture. "
      "Ranked RoutingMethod objects are used by the MappingManager to route "
      "whole circuits.")
      .def(
          py::init<
              const std::function<std::tuple<Circuit, unit_map_t, unit_map_t>(
                  const Circuit&, const ArchitecturePtr&)>&,
              const std::function<bool(const Circuit&, const ArchitecturePtr&)>,
              unsigned, unsigned>(),
          "Constructor for a routing method defined by partially routing "
          "subcircuits.\n\n:param route_subcircuit: A function declaration "
          "that given a Circuit and Architecture object, returns a tuple "
          "containing a new modified circuit, the initial logical to physical "
          "qubit mapping of the modified circuit and the permutation of "
          "'logical to physical qubit mapping given operations in the "
          "modified circuit\n:param check_subcircuit: A function declaration "
          "that given a Circuit and Architecture object, returns a bool "
          "stating whether the given method can modify the "
          "given circuit\n:param max_size: The maximum number of gates "
          "permitted in a subcircuit\n:param max_depth: The maximum permitted "
          "depth of a subcircuit.",
          py::arg("route_subcircuit"), py::arg("check_subcircuit"),
          py::arg("max_size"), py::arg("max_depth"));

  py::class_<
      LexiRouteRoutingMethod, std::shared_ptr<LexiRouteRoutingMethod>,
      RoutingMethod>(
      m, "LexiRouteRoutingMethod",
      "Defines a RoutingMethod object for mapping circuits that uses the "
      "Lexicographical Comparison approach outlined in arXiv:1902.08091.")
      .def(
          py::init<unsigned>(),
          "LexiRoute constructor.\n\n:param lookahead: Maximum depth of "
          "lookahead "
          "employed when picking SWAP for purpose of logical to physical "
          "mapping.",
          py::arg("lookahead") = 10);

  py::class_<MappingManager>(
      m, "MappingManager",
      "Defined by a pytket Architecture object, maps Circuit logical Qubits "
      "to Physically permitted Architecture qubits. Mapping is completed by "
      "sequential routing (full or partial) of subcircuits. Custom method for "
      "routing (full or partial) of subcircuits can be defined in python "
      "layer.")
      .def(
          py::init<const ArchitecturePtr&>(),
          "MappingManager constructor.\n\n:param architecture: pytket "
          "Architecure object MappingManager object defined by.",
          py::arg("architecture"))
      .def(
          "route_circuit", &MappingManager::route_circuit,
          "Maps from given logical circuit to physical circuit. Modification "
          "defined by route_subcircuit, but typically this proceeds by "
          "insertion of SWAP gates that permute logical qubits on physical "
          "qubits. \n\n:param circuit: pytket circuit to be mapped"
          "\n:param routing_methods: Ranked methods to use for routing "
          "subcircuits. In given order, each method is sequentially checked "
          "for viability, with the first viable method being used.",
          py::arg("circuit"), py::arg("routing_methods"));
}
}  // namespace tket