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

#include "Circuit/Boxes.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Circuit/Circuit.hpp"
#include "Converters/PhasePoly.hpp"
#include "Utils/Json.hpp"
#include "binder_json.hpp"
#include "binder_utils.hpp"
#include "typecast.hpp"

namespace py = pybind11;
using json = nlohmann::json;

namespace tket {

void init_boxes(py::module &m) {
  py::class_<CircBox, std::shared_ptr<CircBox>, Op>(
      m, "CircBox",
      "A user-defined operation specified by a :py:class:`Circuit`.")
      .def(
          py::init<const Circuit &>(), "Construct from a :py:class:`Circuit`.",
          py::arg("circ"))
      .def(
          "get_circuit", [](CircBox &cbox) { return *cbox.to_circuit(); },
          ":return: the :py:class:`Circuit` described by the box");
  py::class_<Unitary1qBox, std::shared_ptr<Unitary1qBox>, Op>(
      m, "Unitary1qBox",
      "A user-defined one-qubit operation specified by a unitary matrix.")
      .def(
          py::init<const Eigen::Matrix2cd &>(),
          "Construct from a unitary matrix.", py::arg("m"))
      .def(
          "get_circuit", [](Unitary1qBox &ubox) { return *ubox.to_circuit(); },
          ":return: the :py:class:`Circuit` described by the box")
      .def(
          "get_matrix", &Unitary1qBox::get_matrix,
          ":return: the unitary matrix as a numpy array");
  py::class_<Unitary2qBox, std::shared_ptr<Unitary2qBox>, Op>(
      m, "Unitary2qBox",
      "A user-defined two-qubit operation specified by a unitary matrix.")
      .def(
          py::init<const Eigen::Matrix4cd &, BasisOrder>(),
          "Construct from a unitary matrix.\n\n"
          ":param m: The unitary matrix\n"
          ":param basis: Whether the provided unitary is in the ILO-BE "
          "(increasing lexicographic order of qubit ids, big-endian "
          "indexing) format, or DLO-BE (decreasing lexicographic order "
          "of ids)",
          py::arg("m"), py::arg("basis") = BasisOrder::ilo)
      .def(
          "get_circuit", [](Unitary2qBox &ubox) { return *ubox.to_circuit(); },
          ":return: the :py:class:`Circuit` described by the box")
      .def(
          "get_matrix", &Unitary2qBox::get_matrix,
          ":return: the unitary matrix (in ILO-BE format) as a numpy "
          "array");
  py::class_<Unitary3qBox, std::shared_ptr<Unitary3qBox>, Op>(
      m, "Unitary3qBox",
      "A user-defined three-qubit operation specified by a unitary matrix.")
      .def(
          py::init<const Eigen::MatrixXcd &, BasisOrder>(),
          "Construct from a unitary matrix.\n\n"
          ":param m: The unitary matrix\n"
          ":param basis: Whether the provided unitary is in the ILO-BE "
          "(increasing lexicographic order of qubit ids, big-endian "
          "indexing) format, or DLO-BE (decreasing lexicographic order "
          "of ids)",
          py::arg("m"), py::arg("basis") = BasisOrder::ilo)
      .def(
          "get_circuit", [](Unitary3qBox &ubox) { return *ubox.to_circuit(); },
          ":return: the :py:class:`Circuit` described by the box")
      .def(
          "get_matrix", &Unitary3qBox::get_matrix,
          ":return: the unitary matrix (in ILO-BE format) as a numpy array");
  py::class_<ExpBox, std::shared_ptr<ExpBox>, Op>(
      m, "ExpBox",
      "A user-defined two-qubit operation whose corresponding unitary "
      "matrix "
      "is the exponential of a user-defined hermitian matrix.")
      .def(
          py::init<const Eigen::Matrix4cd &, double, BasisOrder>(),
          "Construct :math:`e^{itA}` from a hermitian matrix :math:`A` "
          "and a parameter :math:`t`.\n\n"
          ":param A: A hermitian matrix\n"
          ":param t: Exponentiation parameter\n"
          ":param basis: Whether the provided matrix is in the ILO-BE "
          "(increasing lexicographic order of qubit ids, big-endian "
          "indexing) format, or DLO-BE (decreasing lexicographic order "
          "of ids)",
          py::arg("A"), py::arg("t"), py::arg("basis") = BasisOrder::ilo)
      .def(
          "get_circuit", [](ExpBox &ebox) { return *ebox.to_circuit(); },
          ":return: the :py:class:`Circuit` described by the box");
  py::class_<PauliExpBox, std::shared_ptr<PauliExpBox>, Op>(
      m, "PauliExpBox",
      "An operation defined as the exponential of a tensor of Pauli "
      "operations and a (possibly symbolic) phase parameter.")
      .def(
          py::init<const std::vector<Pauli> &, Expr>(),
          "Construct :math:`e^{-\\frac12 i \\pi t \\sigma_0 \\otimes "
          "\\sigma_1 \\otimes \\cdots}` from Pauli operators "
          ":math:`\\sigma_i \\in \\{I,X,Y,Z\\}` and a parameter "
          ":math:`t`.",
          py::arg("paulis"), py::arg("t"))
      .def(
          "get_circuit", [](PauliExpBox &pbox) { return *pbox.to_circuit(); },
          ":return: the :py:class:`Circuit` described by the box")
      .def(
          "get_paulis", &PauliExpBox::get_paulis,
          ":return: the corresponding list of " CLSOBJS(Pauli))
      .def(
          "get_phase", &PauliExpBox::get_phase,
          ":return: the corresponding phase parameter");
  py::class_<QControlBox, std::shared_ptr<QControlBox>, Op>(
      m, "QControlBox",
      "A user-defined controlled operation specified by an "
      ":py:class:`Op` and the number of quantum controls.")
      .def(
          py::init<Op_ptr &, unsigned>(),
          "Construct from an :py:class:`Op` and a number of quantum "
          "controls. The controls occupy the low-index ports of the "
          "resulting operation.",
          py::arg("op"), py::arg("n") = 1)
      .def(
          "get_circuit", [](QControlBox &qcbox) { return *qcbox.to_circuit(); },
          ":return: the :py:class:`Circuit` described by the box")
      .def("get_op", &QControlBox::get_op, ":return: the underlying operator")
      .def(
          "get_n_controls", &QControlBox::get_n_controls,
          ":return: the number of control qubits");

  py::class_<CompositeGateDef, composite_def_ptr_t>(
      m, "CustomGateDef",
      "A custom unitary gate definition, given as a composition of other "
      "gates")
      .def(py::init<
           const std::string &, const Circuit &, const std::vector<Sym> &>())
      .def_static(
          "define", &CompositeGateDef::define_gate,
          "Define a new custom gate as a composite of other "
          "gates\n\n:param name: Readable name for the new "
          "gate\n:param circ: The definition of the gate as a "
          "Circuit\n:param args: Symbols to be encapsulated as "
          "arguments of the custom gate",
          py::arg("name"), py::arg("circ"), py::arg("args"))
      .def_property_readonly(
          "name", &CompositeGateDef::get_name, "The readable name of the gate")
      .def_property_readonly(
          "definition", &CompositeGateDef::get_def,
          "Return definition as a circuit.")
      .def_property_readonly(
          "args", &CompositeGateDef::get_args,
          "Return symbolic arguments of gate.")
      .def_property_readonly(
          "arity", &CompositeGateDef::n_args,
          "The number of real parameters for the gate")
      .def(
          "to_dict",
          [](const CompositeGateDef &c) {
            return json(std::make_shared<CompositeGateDef>(c));
          },
          ":return: a JSON serializable dictionary representation of "
          "the CustomGateDef")
      .def_static(
          "from_dict",
          [](const json &j) { return j.get<composite_def_ptr_t>(); },
          "Construct Circuit instance from JSON serializable "
          "dictionary representation of the Circuit.");
  py::class_<CustomGate, std::shared_ptr<CustomGate>, Op>(
      m, "CustomGate",
      "A user-defined gate defined by a parametrised :py:class:`Circuit`.")
      .def(
          py::init<const composite_def_ptr_t &, const std::vector<Expr> &>(),
          "Instantiate a custom gate.", py::arg("gatedef"), py::arg("params"))
      .def_property_readonly(
          "name", [](CustomGate &cgate) { return cgate.get_name(false); },
          "The readable name of the gate.")
      .def_property_readonly(
          "params", &CustomGate::get_params, "The parameters of the gate.")
      .def_property_readonly(
          "gate", &CustomGate::get_gate, "Underlying gate object.")
      .def(
          "get_circuit",
          [](CustomGate &composite) { return *composite.to_circuit(); },
          ":return: the :py:class:`Circuit` described by the gate.");
  py::class_<PhasePolyBox, std::shared_ptr<PhasePolyBox>, Op>(
      m, "PhasePolyBox",
      "Box encapsulating any Circuit made up of CNOT and RZ as a phase "
      "polynomial + linear transformation")
      .def(
          py::init([](unsigned n_qb, const std::map<Qubit, unsigned> &q_ind,
                      const PhasePolynomial &p_p, const MatrixXb &lin_trans) {
            boost::bimap<Qubit, unsigned> bmap;
            for (const auto &pair : q_ind) {
              bmap.insert({pair.first, pair.second});
            }

            return PhasePolyBox(n_qb, bmap, p_p, lin_trans);
          }),
          "Construct from the number of qubits, the mapping from "
          "Qubit to index, the phase polynomial (map from bitstring "
          "to phase) and the linear transformation (boolean matrix)",
          py::arg("n_qubits"), py::arg("qubit_indices"),
          py::arg("phase_polynomial"), py::arg("linear_transformation"))
      .def(
          py::init([](const Circuit &circ) { return PhasePolyBox(circ); }),
          "Construct a PhasePolyBox from a given circuit containing only Rz "
          "and CX gates.",
          py::arg("circuit"))
      .def_property_readonly(
          "n_qubits", &PhasePolyBox::get_n_qubits,
          "Number of gates the polynomial acts on.")
      .def_property_readonly(
          "phase_polynomial",
          [](PhasePolyBox &ppoly) {
            const PhasePolynomial &phase_pol = ppoly.get_phase_polynomial();
            std::map<py::tuple, Expr> outmap;
            for (const auto &pair : phase_pol) {
              outmap.insert({py::tuple(py::cast(pair.first)), pair.second});
            }
            return outmap;
          },
          "Map from bitstring (basis state) to phase.")
      .def_property_readonly(
          "linear_transformation", &PhasePolyBox::get_linear_transformation,
          "Boolean matrix corresponding to linear transformation.")
      .def(
          "get_circuit", [](PhasePolyBox &ppb) { return *ppb.to_circuit(); },
          ":return: the :py:class:`Circuit` described by the box.")
      .def_property_readonly(
          "qubit_indices",
          [](PhasePolyBox &ppoly) {
            const boost::bimap<Qubit, unsigned> &bmap =
                ppoly.get_qubit_indices();
            std::map<Qubit, unsigned> outmap;
            for (const auto &pair : bmap.left) {
              outmap.insert({pair.first, pair.second});
            }
            return outmap;
          },
          "Map from Qubit to index in polynomial.");
  py::class_<ProjectorAssertionBox, std::shared_ptr<ProjectorAssertionBox>, Op>(
      m, "ProjectorAssertionBox",
      "A user-defined assertion specified by a 2x2, 4x4, or 8x8 projector "
      "matrix.")
      .def(
          py::init<const Eigen::MatrixXcd &, BasisOrder>(),
          "Construct from a projector matrix.\n\n"
          ":param m: The projector matrix\n"
          ":param basis: Whether the provided unitary is in the ILO-BE "
          "(increasing lexicographic order of qubit ids, big-endian "
          "indexing) format, or DLO-BE (decreasing lexicographic order "
          "of ids)",
          py::arg("m"), py::arg("basis") = BasisOrder::ilo)
      .def(
          "get_circuit",
          [](ProjectorAssertionBox &ubox) { return *ubox.to_circuit(); },
          ":return: the :py:class:`Circuit` described by the box")
      .def(
          "get_matrix", &ProjectorAssertionBox::get_matrix,
          ":return: the unitary matrix (in ILO-BE format) as a numpy array");
  py::class_<
      StabiliserAssertionBox, std::shared_ptr<StabiliserAssertionBox>, Op>(
      m, "StabiliserAssertionBox",
      "A user-defined assertion specified by a list of pauli stabilisers.")
      .def(
          py::init<const PauliStabiliserList>(),
          "Construct from a list of pauli stabilisers.\n\n"
          ":param m: The list of pauli stabilisers\n",
          py::arg("stabilisers"))
      .def(
          py::init([](const std::vector<std::string> &pauli_strings) {
            PauliStabiliserList stabilisers;
            for (auto &raw_string : pauli_strings) {
              std::vector<Pauli> string;
              bool coeff = true;
              for (unsigned i = 0; i < raw_string.size(); i++) {
                switch (raw_string[i]) {
                  case '-':
                    if (i == 0) {
                      coeff = false;
                    } else {
                      throw std::invalid_argument(
                          "Invalid pauli string: " + raw_string);
                    }
                    break;
                  case 'I':
                    string.push_back(Pauli::I);
                    break;
                  case 'X':
                    string.push_back(Pauli::X);
                    break;
                  case 'Y':
                    string.push_back(Pauli::Y);
                    break;
                  case 'Z':
                    string.push_back(Pauli::Z);
                    break;
                  default:
                    throw std::invalid_argument(
                        "Invalid pauli string: " + raw_string);
                }
              }
              stabilisers.push_back(PauliStabiliser(string, coeff));
            }
            return StabiliserAssertionBox(stabilisers);
          }),
          "Construct from a list of pauli stabilisers.\n\n"
          ":param m: The list of pauli stabilisers expressed as Python "
          "strings\n",
          py::arg("stabilisers"))
      .def(
          "get_circuit",
          [](StabiliserAssertionBox &ubox) { return *ubox.to_circuit(); },
          ":return: the :py:class:`Circuit` described by the box")
      .def(
          "get_stabilisers", &StabiliserAssertionBox::get_stabilisers,
          ":return: the list of pauli stabilisers");
}
}  // namespace tket
