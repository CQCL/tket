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

#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Utils/PauliStrings.hpp"
#include "binder_json.hpp"
#include "typecast.hpp"

namespace py = pybind11;
using json = nlohmann::json;

namespace tket {

PYBIND11_MODULE(pauli, m) {
  py::enum_<Pauli>(m, "Pauli")
      .value("I", Pauli::I)
      .value("X", Pauli::X)
      .value("Y", Pauli::Y)
      .value("Z", Pauli::Z)
      .export_values();

  py::class_<QubitPauliString>(
      m, "QubitPauliString",
      "A string of Pauli letters from the alphabet {I, X, Y, Z},"
      "implemented as a sparse list, indexed by qubit.")
      .def(py::init<>(), "Constructs an empty QubitPauliString.")
      .def(
          py::init<Qubit, Pauli>(),
          "Constructs a QubitPauliString with a single Pauli term.",
          py::arg("qubit"), py::arg("pauli"))
      .def(
          py::init<std::list<Qubit>, std::list<Pauli>>(),
          "Constructs a QubitPauliString from two matching lists of "
          "Qubits and Paulis.",
          py::arg("qubits"), py::arg("paulis"))
      .def(
          py::init<QubitPauliMap>(),
          "Construct a QubitPauliString from a QubitPauliMap.", py::arg("map"))
      .def(
          "__hash__",
          [](const QubitPauliString &qps) { return hash_value(qps); })
      .def("__repr__", &QubitPauliString::to_str)
      .def("__eq__", &QubitPauliString::operator==)
      .def("__ne__", &QubitPauliString::operator!=)
      .def("__lt__", &QubitPauliString::operator<)
      .def("__getitem__", &QubitPauliString::get)
      .def("__setitem__", &QubitPauliString::set)
      .def_property_readonly(
          "map", [](const QubitPauliString &qps) { return qps.map; },
          ":return: the QubitPauliString's underlying dict mapping "
          ":py:class:`Qubit` to :py:class:`Pauli`")
      .def(
          "to_list",
          [](const QubitPauliString &qps) {
            json j = qps;
            return j;
          },
          "A JSON-serializable representation of the QubitPauliString. "
          "\n:return: a list of :py:class:`Qubit`-to-:py:class:`Pauli` "
          "entries, "
          "represented as dicts.")
      .def_static(
          "from_list", [](const json &j) { return j.get<QubitPauliString>(); },
          "Construct a new QubitPauliString instance from a JSON serializable "
          "list "
          "representation.")
      .def(
          "compress", &QubitPauliString::compress,
          "Removes I terms to compress the sparse representation.")
      .def(
          "commutes_with", &QubitPauliString::commutes_with,
          ":return: True if the two strings commute, else False",
          py::arg("other"))
      .def(
          "to_sparse_matrix",
          (CmplxSpMat(QubitPauliString::*)(void) const) &
              QubitPauliString::to_sparse_matrix,
          "Represents the sparse string as a dense string (without "
          "padding for extra qubits) and generates the matrix for the "
          "tensor. Uses the ILO-BE convention, so ``Qubit(\"a\", 0)`` "
          "is more significant that ``Qubit(\"a\", 1)`` and "
          "``Qubit(\"b\")`` for indexing into the matrix."
          "\n\n:return: a sparse matrix corresponding to the operator")
      .def(
          "to_sparse_matrix",
          (CmplxSpMat(QubitPauliString::*)(const unsigned) const) &
              QubitPauliString::to_sparse_matrix,
          "Represents the sparse string as a dense string over "
          "`n_qubits` qubits (sequentially indexed from 0 in the "
          "default register) and generates the matrix for the tensor. "
          "Uses the ILO-BE convention, so ``Qubit(0)`` is the most "
          "significant bit for indexing into the matrix."
          "\n\n:param n_qubits: the number of qubits in the full "
          "operator"
          "\n:return: a sparse matrix corresponding to the operator",
          py::arg("n_qubits"))
      .def(
          "to_sparse_matrix",
          (CmplxSpMat(QubitPauliString::*)(const qubit_vector_t &) const) &
              QubitPauliString::to_sparse_matrix,
          "Represents the sparse string as a dense string and generates "
          "the matrix for the tensor. Orders qubits according to "
          "`qubits` (padding with identities if they are not in the "
          "sparse string), so ``qubits[0]`` is the most significant bit "
          "for indexing into the matrix."
          "\n\n:param qubits: the ordered list of qubits in the full "
          "operator"
          "\n:return: a sparse matrix corresponding to the operator",
          py::arg("qubits"))
      .def(
          "dot_state",
          (Eigen::VectorXcd(QubitPauliString::*)(const Eigen::VectorXcd &)
               const) &
              QubitPauliString::dot_state,
          "Performs the dot product of the state with the pauli string. "
          "Maps the qubits of the statevector with sequentially-indexed "
          "qubits in the default register, with ``Qubit(0)`` being the "
          "most significant qubit."
          "\n\n:param state: statevector for qubits ``Qubit(0)`` to "
          "``Qubit(n-1)``"
          "\n:return: dot product of operator with state",
          py::arg("state"))
      .def(
          "dot_state",
          (Eigen::VectorXcd(QubitPauliString::*)(
              const Eigen::VectorXcd &, const qubit_vector_t &) const) &
              QubitPauliString::dot_state,
          "Performs the dot product of the state with the pauli string. "
          "Maps the qubits of the statevector according to the ordered "
          "list `qubits`, with ``qubits[0]`` being the most significant "
          "qubit."
          "\n\n:param state: statevector"
          "\n:param qubits: order of qubits in `state` from most to "
          "least significant"
          "\n:return: dot product of operator with state",
          py::arg("state"), py::arg("qubits"))
      .def(
          "state_expectation",
          (Complex(QubitPauliString::*)(const Eigen::VectorXcd &) const) &
              QubitPauliString::state_expectation,
          "Calculates the expectation value of the state with the pauli "
          "string. Maps the qubits of the statevector with "
          "sequentially-indexed qubits in the default register, with "
          "``Qubit(0)`` being the most significant qubit."
          "\n\n:param state: statevector for qubits ``Qubit(0)`` to "
          "``Qubit(n-1)``"
          "\n:return: expectation value with respect to state",
          py::arg("state"))
      .def(
          "state_expectation",
          (Complex(QubitPauliString::*)(
              const Eigen::VectorXcd &, const qubit_vector_t &) const) &
              QubitPauliString::state_expectation,
          "Calculates the expectation value of the state with the pauli "
          "string. Maps the qubits of the statevector according to the "
          "ordered list `qubits`, with ``qubits[0]`` being the most "
          "significant qubit."
          "\n\n:param state: statevector"
          "\n:param qubits: order of qubits in `state` from most to "
          "least significant"
          "\n:return: expectation value with respect to state",
          py::arg("state"), py::arg("qubits"))

      .def(py::pickle(
          [](const QubitPauliString &qps) {
            /* Hackery to avoid pickling an opaque object */
            std::list<Qubit> qubits;
            std::list<Pauli> paulis;
            for (const std::pair<const Qubit, Pauli> &qp_pair : qps.map) {
              qubits.push_back(qp_pair.first);
              paulis.push_back(qp_pair.second);
            }
            return py::make_tuple(qubits, paulis);
          },
          [](const py::tuple &t) {
            if (t.size() != 2)
              throw std::runtime_error(
                  "Invalid state: tuple size: " + std::to_string(t.size()));
            return QubitPauliString(
                t[0].cast<std::list<Qubit>>(), t[1].cast<std::list<Pauli>>());
          }));

  m.def(
      "pauli_string_mult",
      [](const QubitPauliString &qps1, const QubitPauliString &qps2) {
        QubitPauliTensor product_tensor =
            QubitPauliTensor(qps1) * QubitPauliTensor(qps2);
        return std::pair<QubitPauliString, Complex>(
            product_tensor.string, product_tensor.coeff);
      },
      ":return: the product of two QubitPauliString objects as a pair "
      "(QubitPauliString, complex)",
      py::arg("qubitpaulistring1"), py::arg("qubitpaulistring2"));

  py::class_<PauliStabiliser>(
      m, "PauliStabiliser",
      "A string of Pauli letters from the alphabet {I, X, Y, Z} "
      "with a +/- 1 coefficient.")
      .def(py::init<>(), "Constructs an empty QubitPauliString.")
      .def(
          py::init([](const std::vector<Pauli> &string, const int &coeff) {
            if (coeff == 1) {
              return PauliStabiliser(string, true);
            }
            if (coeff == -1) {
              return PauliStabiliser(string, false);
            }
            throw std::invalid_argument("Coefficient must be -1 or 1.");
          }),
          "Constructs a PauliStabiliser with a list of Pauli terms.",
          py::arg("string"), py::arg("coeff"))
      .def_property_readonly(
          "coeff",
          [](const PauliStabiliser &stabiliser) {
            if (stabiliser.coeff) return 1;
            return -1;
          },
          "The coefficient of the stabiliser")
      .def_property_readonly(
          "string",
          [](const PauliStabiliser &stabiliser) { return stabiliser.string; },
          "The list of Pauli terms")
      .def("__eq__", &PauliStabiliser::operator==)
      .def("__ne__", &PauliStabiliser::operator!=);
}

}  // namespace tket
