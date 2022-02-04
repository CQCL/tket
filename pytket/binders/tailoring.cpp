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

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "Characterisation/FrameRandomisation.hpp"
#include "Clifford/CliffTableau.hpp"
#include "Converters/Converters.hpp"
#include "binder_utils.hpp"

namespace py = pybind11;

namespace tket {

QubitPauliString apply_clifford_basis_change(
    const QubitPauliString &in_pauli, const Circuit &circ) {
  CliffTableau tab = circuit_to_tableau(circ);
  QubitPauliTensor new_operator;
  for (const auto &qp : in_pauli.map) {
    switch (qp.second) {
      case Pauli::I:
        break;
      case Pauli::X:
        new_operator = new_operator * tab.get_xpauli(qp.first);
        break;
      case Pauli::Y:
        new_operator =
            new_operator * tab.get_xpauli(qp.first) * tab.get_zpauli(qp.first);
        break;
      case Pauli::Z:
        new_operator = new_operator * tab.get_zpauli(qp.first);
        break;
    }
  }
  return new_operator.string;
}

PYBIND11_MODULE(tailoring, m) {
  m.doc() = "The tailoring module provides access to noise tailoring tools.";

  py::class_<FrameRandomisation>(
      m, "FrameRandomisation",
      "The base FrameRandomisation class. FrameRandomisation finds "
      "subcircuits (cycles) of a given circuit comprised of gates with "
      "OpType only from a specified set of OpType, and wires gates into "
      "the boundary (frame) of these cycles. Input frame gates are "
      "sampled from another set of OpType, and output frame gates "
      "deduced such that the circuit unitary doesn't change, achieved by "
      "computing the action of cycle gates on frame gates.")
      .def(
          py::init([](const OpTypeSet &_cycle_types,
                      const OpTypeSet &_frame_types,
                      const std::map<OpType, std::map<py::tuple, py::tuple>>
                          &cycle_frame_actions) {
            std::map<OpType, std::map<OpTypeVector, OpTypeVector>>
                real_cycle_frame_actions;
            for (const auto &actions : cycle_frame_actions) {
              std::map<OpTypeVector, OpTypeVector> otv_map;
              for (const auto &tups : actions.second) {
                otv_map[tups.first.cast<OpTypeVector>()] =
                    tups.second.cast<OpTypeVector>();
              }
              real_cycle_frame_actions[actions.first] = otv_map;
            }
            return FrameRandomisation(
                _cycle_types, _frame_types, real_cycle_frame_actions);
          }),
          "Constructor for FrameRandomisation."
          "\n\n:param cycletypes: "
          "A set of OpType corresponding to the gates cycles found are "
          "comprised of"
          "\n:param frametypes: A set of OpType "
          "corresponding to the gates Frames are sampled from"
          "\n:param conjugates: A map from cycle OpType, to a map "
          "between Frame OptypeVector giving the required change to "
          "output frame OpType to preserve Unitary from given input "
          "frame OpType.")
      .def(
          "get_all_circuits", &FrameRandomisation::get_all_circuits,
          "For given circuit, finds all Cycles, finds all frames for "
          "each Cycle, and returns every combination of frame and cycle "
          "in a vector of Circuit.\n\n:param circuit: The circuit to "
          "find frames for.\n"
          ":return: list of " CLSOBJS(Circuit),
          py::arg("circuit"))
      .def(
          "sample_circuits", &FrameRandomisation::sample_randomisation_circuits,
          "Returns a number of instances equal to sample of frame "
          "randomisation for the given circuit. Samples individual "
          "frame gates uniformly.\n\n:param "
          "circuit: The circuit to perform frame randomisation with "
          "Pauli gates on\n:param samples: the number of frame "
          "randomised circuits to return.\n"
          ":return: list of " CLSOBJS(Circuit),
          py::arg("circuit"), py::arg("samples"))
      .def("__repr__", &FrameRandomisation::to_string);

  py::class_<PauliFrameRandomisation>(
      m, "PauliFrameRandomisation",
      "The PauliFrameRandomisation class. PauliFrameRandomisation finds "
      "subcircuits (cycles) of a given circuit comprised of gates with "
      "OpType::H, OpType::CX and OpType::S, and wires gates "
      "into "
      "the boundary (frame) of these cycles. Input frame gates are "
      "sampled from another set of OpType comprised of the Pauli gates, "
      "and output frame gates "
      "deduced such that the circuit unitary doesn't change, achieved by "
      "computing the action of cycle gates on frame gates.")
      .def(py::init<>(), "Constructor for PauliFrameRandomisation.")
      .def(
          "get_all_circuits", &PauliFrameRandomisation::get_all_circuits,
          "For given circuit, finds all Cycles, finds all frames for "
          "each Cycle, and returns every combination of frame and cycle "
          "in a vector of Circuit.\n\n:param circuit: The circuit to "
          "find frames for.\n"
          ":return: list of " CLSOBJS(Circuit),
          py::arg("circuit"))
      .def(
          "sample_circuits",
          &PauliFrameRandomisation::sample_randomisation_circuits,
          "Returns a number of instances equal to sample of frame "
          "randomisation for the given circuit. Samples individual "
          "frame gates uniformly from the Pauli gates.\n\n:param "
          "circuit: The circuit to perform frame randomisation with "
          "Pauli gates on\n:param samples: the number of frame "
          "randomised circuits to return.\n"
          ":return: list of " CLSOBJS(Circuit),
          py::arg("circuit"), py::arg("samples"))
      .def("__repr__", &PauliFrameRandomisation::to_string);

  py::class_<UniversalFrameRandomisation>(
      m, "UniversalFrameRandomisation",
      "The UniversalFrameRandomisation class. "
      "UniversalFrameRandomisation finds "
      "subcircuits (cycles) of a given circuit comprised of gates with "
      "OpType::H, OpType::CX, and OpType::Rz, and wires gates "
      "into "
      "the boundary (frame) of these cycles. Input frame gates are "
      "sampled from another set of OpType comprised of the Pauli gates, "
      "and output frame gates "
      "deduced such that the circuit unitary doesn't change, achieved by "
      "computing the action of cycle gates on frame gates. Some gates "
      "with OpType::Rz may be substituted for their dagger to achieve "
      "this.")
      .def(py::init<>(), "Constructor for UniversalFrameRandomisation.")
      .def(
          "get_all_circuits", &UniversalFrameRandomisation::get_all_circuits,
          "For given circuit, finds all Cycles, finds all frames for "
          "each Cycle, and returns every combination of frame and cycle "
          "in a vector of Circuit.\n\n:param circuit: The circuit to "
          "find frames for.\n"
          ":return: list of " CLSOBJS(Circuit),
          py::arg("circuit"))
      .def(
          "sample_circuits",
          &UniversalFrameRandomisation::sample_randomisation_circuits,
          "Returns a number of instances equal to sample of frame "
          "randomisation for the given circuit. Samples individual "
          "frame gates uniformly from the Pauli gates.\n\n:param "
          "circuit: The circuit to perform frame randomisation with "
          "Pauli gates on\n:param samples: the number of frame "
          "randomised circuits to return.\n"
          ":return: list of " CLSOBJS(Circuit),
          py::arg("circuit"), py::arg("samples"))
      .def("__repr__", &UniversalFrameRandomisation::to_string);
  m.def(
      "apply_clifford_basis_change", &apply_clifford_basis_change,
      "Given Pauli operator P and Clifford circuit C, "
      "returns C_dagger.P.C in multiplication order. This ignores any -1 "
      "phase that could be introduced.\n\n:param pauli: Pauli "
      "operator being transformed. \n:param circuit: "
      "Clifford circuit acting on Pauli operator.\n"
      ":return: :py:class:`QubitPauliString` for new operator",
      py::arg("pauli"), py::arg("circuit"));
}
}  // namespace tket
