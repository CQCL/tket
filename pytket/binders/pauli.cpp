// Copyright Quantinuum
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
#include "deleted_hash.hpp"
#include "py_operators.hpp"
#include "tket/Utils/PauliTensor.hpp"
#include "typecast.hpp"

namespace py = pybind11;
using json = nlohmann::json;

namespace tket {

typedef py::tket_custom::SequenceVec<Qubit> py_qubit_vector_t;
PYBIND11_MODULE(pauli, m) {
  py::module::import("pytket._tket.unit_id");
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
          py::init<
              py::tket_custom::SequenceList<Qubit>,
              py::tket_custom::SequenceList<Pauli>>(),
          "Constructs a QubitPauliString from two matching lists of "
          "Qubits and Paulis.",
          py::arg("qubits"), py::arg("paulis"))
      .def(
          py::init<QubitPauliMap>(),
          "Construct a QubitPauliString from a dictionary mapping "
          ":py:class:`Qubit` to :py:class:`Pauli`.",
          py::arg("map"))
      .def(
          "__hash__", [](const SpPauliString &qps) { return qps.hash_value(); })
      .def("__repr__", &SpPauliString::to_str)
      .def("__eq__", &py_equals<SpPauliString>)
      .def("__ne__", &py_not_equals<SpPauliString>)
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
            return py::object(json(qps.string)).cast<py::list>();
          },
          "A JSON-serializable representation of the QubitPauliString.\n\n"
          ":return: a list of :py:class:`Qubit`-to-:py:class:`Pauli` "
          "entries, "
          "represented as dicts.")
      .def_static(
          "from_list",
          [](const py::list &qubit_pauli_string_list) {
            return SpPauliString(
                json(qubit_pauli_string_list).get<QubitPauliMap>());
          },
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
          [](const SpPauliString &self, const py_qubit_vector_t &qubits) {
            return self.to_sparse_matrix(qubits);
          },
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
          [](const SpPauliString &self, const Eigen::VectorXcd &state,
             const py_qubit_vector_t &qubits) {
            return self.dot_state(state, qubits);
          },
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
          [](const SpPauliString &self, const Eigen::VectorXcd &state) {
            return self.state_expectation(state);
          },
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
          [](const SpPauliString &self, const Eigen::VectorXcd &state,
             const py_qubit_vector_t &qubits) {
            return self.state_expectation(state, qubits);
          },
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
      .def(py::init<>(), "Constructs an empty PauliStabiliser.")
      .def(
          py::init([](const py::tket_custom::SequenceVec<Pauli> &string,
                      const int &coeff) {
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
            return stabiliser.is_real_negative() ? -1 : 1;
          },
          "The coefficient of the stabiliser")
      .def_property_readonly(
          "string",
          [](const PauliStabiliser &stabiliser) { return stabiliser.string; },
          "The list of Pauli terms")
      .def("__eq__", &py_equals<PauliStabiliser>)
      .def("__hash__", &deletedHash<PauliStabiliser>, deletedHashDocstring)
      .def("__ne__", &py_not_equals<PauliStabiliser>);

  py::class_<SpCxPauliTensor>(
      m, "QubitPauliTensor",
      "A tensor formed by Pauli terms, consisting of a sparse map from "
      ":py:class:`Qubit` to :py:class:`Pauli` (implemented as a "
      ":py:class:`QubitPauliString`) and a complex coefficient.")
      .def(
          py::init(
              [](const Complex &coeff) { return SpCxPauliTensor({}, coeff); }),
          "Constructs an empty QubitPauliTensor, representing the identity.",
          py::arg("coeff") = 1.)
      .def(
          py::init<Qubit, Pauli, Complex>(),
          "Constructs a QubitPauliTensor with a single Pauli term.",
          py::arg("qubit"), py::arg("pauli"), py::arg("coeff") = 1.)
      .def(
          py::init([](const py::tket_custom::SequenceList<Qubit> &qubits,
                      const py::tket_custom::SequenceList<Pauli> &paulis,
                      const Complex &coeff) {
            return SpCxPauliTensor(qubits, paulis, coeff);
          }),
          "Constructs a QubitPauliTensor from two matching lists of "
          "Qubits and Paulis.",
          py::arg("qubits"), py::arg("paulis"), py::arg("coeff") = 1.)
      .def(
          py::init<QubitPauliMap, Complex>(),
          "Construct a QubitPauliTensor from a dictionary mapping "
          ":py:class:`Qubit` to :py:class:`Pauli`.",
          py::arg("map"), py::arg("coeff") = 1.)
      .def(
          py::init([](const SpPauliString &qps, const Complex &c) {
            return SpCxPauliTensor(qps.string, c);
          }),
          "Construct a QubitPauliTensor from a QubitPauliString.",
          py::arg("string"), py::arg("coeff") = 1.)
      .def(
          "__hash__",
          [](const SpCxPauliTensor &qps) { return qps.hash_value(); })
      .def("__repr__", &SpCxPauliTensor::to_str)
      .def("__eq__", &py_equals<SpCxPauliTensor>)
      .def("__ne__", &py_not_equals<SpCxPauliTensor>)
      .def("__lt__", &SpCxPauliTensor::operator<)
      .def("__getitem__", &SpCxPauliTensor::get<QubitPauliMap>)
      .def("__setitem__", &SpCxPauliTensor::set<QubitPauliMap>)
      .def(py::self * py::self)
      .def(
          "__rmul__",
          [](const SpCxPauliTensor &qpt, const Complex &c) {
            return SpCxPauliTensor(qpt.string, qpt.coeff * c);
          },
          py::is_operator())
      .def_property(
          "string",
          [](const SpCxPauliTensor &qpt) {
            // Return as SpPauliString for backwards compatibility with before
            // templated PauliTensor
            return SpPauliString(qpt.string);
          },
          [](SpCxPauliTensor &qpt, const SpPauliString &qps) {
            qpt.string = qps.string;
          },
          "The QubitPauliTensor's underlying :py:class:`QubitPauliString`")
      .def_readwrite(
          "coeff", &SpCxPauliTensor::coeff,
          "The global coefficient of the tensor")
      .def(
          "compress", &SpCxPauliTensor::compress<QubitPauliMap>,
          "Removes I terms to compress the sparse representation.")
      .def(
          "commutes_with", &SpCxPauliTensor::commutes_with<Complex>,
          ":return: True if the two tensors commute, else False",
          py::arg("other"))
      .def(
          "to_sparse_matrix",
          [](const SpCxPauliTensor &qpt) { return qpt.to_sparse_matrix(); },
          "Represents the sparse string as a dense string (without "
          "padding for extra qubits) and generates the matrix for the "
          "tensor. Uses the ILO-BE convention, so ``Qubit(\"a\", 0)`` "
          "is more significant that ``Qubit(\"a\", 1)`` and "
          "``Qubit(\"b\")`` for indexing into the matrix."
          "\n\n:return: a sparse matrix corresponding to the tensor")
      .def(
          "to_sparse_matrix",
          [](const SpCxPauliTensor &qpt, unsigned n_qubits) {
            return qpt.to_sparse_matrix(n_qubits);
          },
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
          [](const SpCxPauliTensor &qpt, const py_qubit_vector_t &qubits) {
            return qpt.to_sparse_matrix(qubits);
          },
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
          [](const SpCxPauliTensor &qpt, const Eigen::VectorXcd &state) {
            return qpt.dot_state(state);
          },
          "Performs the dot product of the state with the pauli tensor. "
          "Maps the qubits of the statevector with sequentially-indexed "
          "qubits in the default register, with ``Qubit(0)`` being the "
          "most significant qubit."
          "\n\n:param state: statevector for qubits ``Qubit(0)`` to "
          "``Qubit(n-1)``"
          "\n:return: dot product of operator with state",
          py::arg("state"))
      .def(
          "dot_state",
          [](const SpCxPauliTensor &qpt, const Eigen::VectorXcd &state,
             const py_qubit_vector_t &qubits) {
            return qpt.dot_state(state, qubits);
          },
          "Performs the dot product of the state with the pauli tensor. "
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
          [](const SpCxPauliTensor &qpt, const Eigen::VectorXcd &state) {
            return qpt.state_expectation(state);
          },
          "Calculates the expectation value of the state with the pauli "
          "operator. Maps the qubits of the statevector with "
          "sequentially-indexed qubits in the default register, with "
          "``Qubit(0)`` being the most significant qubit."
          "\n\n:param state: statevector for qubits ``Qubit(0)`` to "
          "``Qubit(n-1)``"
          "\n:return: expectation value with respect to state",
          py::arg("state"))
      .def(
          "state_expectation",
          [](const SpCxPauliTensor &qpt, const Eigen::VectorXcd &state,
             const py_qubit_vector_t &qubits) {
            return qpt.state_expectation(state, qubits);
          },
          "Calculates the expectation value of the state with the pauli "
          "operator. Maps the qubits of the statevector according to the "
          "ordered list `qubits`, with ``qubits[0]`` being the most "
          "significant qubit."
          "\n\n:param state: statevector"
          "\n:param qubits: order of qubits in `state` from most to "
          "least significant"
          "\n:return: expectation value with respect to state",
          py::arg("state"), py::arg("qubits"))

      .def(py::pickle(
          [](const SpCxPauliTensor &qpt) {
            std::list<Qubit> qubits;
            std::list<Pauli> paulis;
            for (const std::pair<const Qubit, Pauli> &qp_pair : qpt.string) {
              qubits.push_back(qp_pair.first);
              paulis.push_back(qp_pair.second);
            }
            return py::make_tuple(qubits, paulis, qpt.coeff);
          },
          [](const py::tuple &t) {
            if (t.size() != 3)
              throw std::runtime_error(
                  "Invalid state: tuple size: " + std::to_string(t.size()));
            return SpCxPauliTensor(
                t[0].cast<std::list<Qubit>>(), t[1].cast<std::list<Pauli>>(),
                t[2].cast<Complex>());
          }));
}
}  // namespace tket
