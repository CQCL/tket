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

#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "binder_json.hpp"
#include "tket/Utils/PauliStrings2.hpp"
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

  py::class_<SpPauliString>(
      m, "QubitPauliString",
      "A string of Pauli letters from the alphabet {I, X, Y, Z}, "
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
          "__hash__", [](const SpPauliString &qps) { return qps.hash_value(); })
      .def("__repr__", &SpPauliString::to_str)
      .def("__eq__", &SpPauliString::operator==)
      .def("__ne__", &SpPauliString::operator!=)
      .def("__lt__", &SpPauliString::operator<)
      .def("__getitem__", &SpPauliString::get<QubitPauliMap>)
      .def("__setitem__", &SpPauliString::set<QubitPauliMap>)
      .def_property_readonly(
          "map", [](const SpPauliString &qps) { return qps.string; },
          ":return: the QubitPauliString's underlying dict mapping "
          ":py:class:`Qubit` to :py:class:`Pauli`")
      .def(
          "to_list",
          [](const SpPauliString &qps) {
            // Just return the QubitPauliMap for backwards compatibility with
            // before templated PauliTensor
            json j = qps.string;
            return j;
          },
          "A JSON-serializable representation of the QubitPauliString.\n\n"
          ":return: a list of :py:class:`Qubit`-to-:py:class:`Pauli` "
          "entries, "
          "represented as dicts.")
      .def_static(
          "from_list",
          [](const json &j) { return SpPauliString(j.get<QubitPauliMap>()); },
          "Construct a new QubitPauliString instance from a JSON serializable "
          "list "
          "representation.")
      .def(
          "compress", &SpPauliString::compress<QubitPauliMap>,
          "Removes I terms to compress the sparse representation.")
      .def(
          "commutes_with", &SpPauliString::commutes_with<no_coeff_t>,
          ":return: True if the two strings commute, else False",
          py::arg("other"))
      .def(
          "to_sparse_matrix",
          (CmplxSpMat(SpPauliString::*)(void) const) &
              SpPauliString::to_sparse_matrix,
          "Represents the sparse string as a dense string (without "
          "padding for extra qubits) and generates the matrix for the "
          "tensor. Uses the ILO-BE convention, so ``Qubit(\"a\", 0)`` "
          "is more significant that ``Qubit(\"a\", 1)`` and "
          "``Qubit(\"b\")`` for indexing into the matrix."
          "\n\n:return: a sparse matrix corresponding to the operator")
      .def(
          "to_sparse_matrix",
          (CmplxSpMat(SpPauliString::*)(const unsigned) const) &
              SpPauliString::to_sparse_matrix,
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
          (CmplxSpMat(SpPauliString::*)(const qubit_vector_t &) const) &
              SpPauliString::to_sparse_matrix,
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
          (Eigen::VectorXcd(SpPauliString::*)(const Eigen::VectorXcd &) const) &
              SpPauliString::dot_state,
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
          (Eigen::VectorXcd(SpPauliString::*)(
              const Eigen::VectorXcd &, const qubit_vector_t &) const) &
              SpPauliString::dot_state,
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
          (Complex(SpPauliString::*)(const Eigen::VectorXcd &) const) &
              SpPauliString::state_expectation,
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
          (Complex(SpPauliString::*)(
              const Eigen::VectorXcd &, const qubit_vector_t &) const) &
              SpPauliString::state_expectation,
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
          [](const SpPauliString &qps) {
            /* Hackery to avoid pickling an opaque object */
            std::list<Qubit> qubits;
            std::list<Pauli> paulis;
            for (const std::pair<const Qubit, Pauli> &qp_pair : qps.string) {
              qubits.push_back(qp_pair.first);
              paulis.push_back(qp_pair.second);
            }
            return py::make_tuple(qubits, paulis);
          },
          [](const py::tuple &t) {
            if (t.size() != 2)
              throw std::runtime_error(
                  "Invalid state: tuple size: " + std::to_string(t.size()));
            return SpPauliString(
                t[0].cast<std::list<Qubit>>(), t[1].cast<std::list<Pauli>>());
          }));

  m.def(
      "pauli_string_mult",
      [](const SpPauliString &qps1, const SpPauliString &qps2) {
        SpCxPauliTensor product_tensor =
            SpCxPauliTensor(qps1) * SpCxPauliTensor(qps2);
        return std::pair<SpPauliString, Complex>(
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
              return PauliStabiliser(string, 0);
            }
            if (coeff == -1) {
              return PauliStabiliser(string, 2);
            }
            throw std::invalid_argument("Coefficient must be -1 or 1.");
          }),
          "Constructs a PauliStabiliser with a list of Pauli terms.",
          py::arg("string"), py::arg("coeff"))
      .def_property_readonly(
          "coeff",
          [](const PauliStabiliser &stabiliser) {
            if (stabiliser.coeff % 4 == 0) return 1;
            if (stabiliser.coeff % 4 == 2) return -1;
            throw std::logic_error(
                "Cannot obtain imaginary coefficient from Pauli Stabiliser.");
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
