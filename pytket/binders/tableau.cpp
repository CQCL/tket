// Copyright 2019-2023 Cambridge Quantum Computing
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

#include <sstream>

#include "tket/Converters/UnitaryTableauBox.hpp"
#include "typecast.hpp"

namespace py = pybind11;
namespace tket {

typedef py::tket_custom::SequenceVec<Qubit> py_qubit_vector_t;
PYBIND11_MODULE(tableau, m) {
  py::class_<UnitaryTableau>(
      m, "UnitaryTableau",
      "Stabilizer tableau for a unitary in the style of Aaronson&Gottesman "
      "\"Improved Simulation of Stabilizer Circuits\": rows indicate the "
      "action at the output corresponding to either an X or a Z on a single "
      "input.")
      .def(
          py::init<unsigned>(),
          "Constructs a :py:class:`UnitaryTableau` representing the identity "
          "operation over some number of qubits. Qubits will be indexed "
          "sequentially in the default register."
          "\n\n:param nqb: The number of qubits in the unitary.",
          py::arg("nqb"))
      .def(
          py::init<
              const MatrixXb&, const MatrixXb&, const VectorXb&,
              const MatrixXb&, const MatrixXb&, const VectorXb&>(),
          "Constructs a :py:class:`UnitaryTableau` from the binary tables of "
          "its components."
          "\n\n:param xx: The X component of the X rows."
          "\n:param xz: The Z component of the X rows."
          "\n:param xph: The phases of the X rows."
          "\n:param zx: The X component of the Z rows."
          "\n:param zz: The Z component of the Z rows."
          "\n:param zph: The phases of the Z rows.",
          py::arg("xx"), py::arg("xz"), py::arg("xph"), py::arg("zx"),
          py::arg("zz"), py::arg("zph"))
      .def(
          py::init<>([](const Circuit& circ) {
            return circuit_to_unitary_tableau(circ);
          }),
          "Constructs a :py:class:`UnitaryTableau` from a unitary "
          ":py:class:`Circuit`. Throws an exception if the input contains "
          "non-unitary operations."
          "\n\n:param circ: The unitary circuit to convert to a tableau.")
      .def(
          "__repr__",
          [](const UnitaryTableau& tab) {
            std::stringstream str;
            str << tab;
            return str.str();
          })
      .def(
          "get_xrow", &UnitaryTableau::get_xrow,
          "Read off an X row as a Pauli string."
          "\n\n:param qb: The qubits whose X row to read off."
          "\n:return: The Pauli string :math:`P` such that :math:`PU=UX_{qb}`.",
          py::arg("qb"))
      .def(
          "get_zrow", &UnitaryTableau::get_zrow,
          "Read off an Z row as a Pauli string."
          "\n\n:param qb: The qubits whose Z row to read off."
          "\n:return: The Pauli string :math:`P` such that :math:`PU=UZ_{qb}`.",
          py::arg("qb"))
      .def(
          "get_row_product", &UnitaryTableau::get_row_product,
          "Combine rows to yield the effect of a given Pauli string."
          "\n\n:param paulis: The Pauli string :math:`P` to consider at the "
          "input."
          "\n:return: The Pauli string :math:`Q` such that :math:`QU=UP`.",
          py::arg("paulis"))
      .def(
          "apply_gate_at_front",
          [](UnitaryTableau& self, const OpType& type,
             const py_qubit_vector_t& qbs) {
            return self.apply_gate_at_front(type, qbs);
          },
          "Update the tableau according to adding a Clifford gate before the "
          "current unitary, i.e. updates :math:`U` to :math:`UG` for a gate "
          ":math:`G`."
          "\n\n:param type: The :py:class:`OpType` of the gate to add. Must be "
          "an unparameterised Clifford gate type."
          "\n:param qbs: The qubits to apply the gate to. Length must match "
          "the arity of the given gate type.",
          py::arg("type"), py::arg("qbs"))
      .def(
          "apply_gate_at_end",
          [](UnitaryTableau& self, const OpType& type,
             const py_qubit_vector_t& qbs) {
            return self.apply_gate_at_end(type, qbs);
          },
          "Update the tableau according to adding a Clifford gate after the "
          "current unitary, i.e. updates :math:`U` to :math:`GU` for a gate "
          ":math:`G`."
          "\n\n:param type: The :py:class:`OpType` of the gate to add. Must be "
          "an unparameterised Clifford gate type."
          "\n:param qbs: The qubits to apply the gate to. Length must match "
          "the arity of the given gate type.",
          py::arg("type"), py::arg("qbs"))
      .def(
          "to_circuit", &unitary_tableau_to_circuit,
          "Synthesises a unitary :py:class:`Circuit` realising the same "
          "unitary as the tableau. Uses the method from Aaronson & Gottesman: "
          "\"Improved Simulation of Stabilizer Circuits\", Theorem 8. This is "
          "not optimised for gate count, so is not recommended for "
          "performance-sensitive usage.");
  py::class_<UnitaryTableauBox, std::shared_ptr<UnitaryTableauBox>, Op>(
      m, "UnitaryTableauBox",
      "A Clifford unitary specified by its actions on Paulis.")
      .def(
          py::init<const UnitaryTableau&>(),
          "Construct from a given tableau.\n\n"
          ":param tab: The :py:class:`UnitaryTableau` representing the desired "
          "unitary.",
          py::arg("tab"))
      .def(
          py::init<
              const MatrixXb&, const MatrixXb&, const VectorXb&,
              const MatrixXb&, const MatrixXb&, const VectorXb&>(),
          "Construct the tableau from the binary tables of its components."
          "\n\n:param xx: The X component of the X rows."
          "\n:param xz: The Z component of the X rows."
          "\n:param xph: The phases of the X rows."
          "\n:param zx: The X component of the Z rows."
          "\n:param zz: The Z component of the Z rows."
          "\n:param zph: The phases of the Z rows.",
          py::arg("xx"), py::arg("xz"), py::arg("xph"), py::arg("zx"),
          py::arg("zz"), py::arg("zph"))
      .def(
          "get_circuit",
          [](UnitaryTableauBox& ubox) { return *ubox.to_circuit(); },
          ":return: The :py:class:`Circuit` described by the box.")
      .def(
          "get_tableau", &UnitaryTableauBox::get_tableau,
          ":return: The tableau representing the unitary operation.");
}

}  // namespace tket
