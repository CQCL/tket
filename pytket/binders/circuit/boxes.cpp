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

#include "tket/Circuit/Boxes.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

#include "binder_json.hpp"
#include "binder_utils.hpp"
#include "tket/Circuit/Circuit.hpp"
#include "tket/Circuit/ConjugationBox.hpp"
#include "tket/Circuit/DiagonalBox.hpp"
#include "tket/Circuit/Multiplexor.hpp"
#include "tket/Circuit/PauliExpBoxes.hpp"
#include "tket/Circuit/StatePreparation.hpp"
#include "tket/Circuit/ToffoliBox.hpp"
#include "tket/Converters/PhasePoly.hpp"
#include "tket/Utils/HelperFunctions.hpp"
#include "tket/Utils/Json.hpp"
#include "typecast.hpp"

namespace py = pybind11;
using json = nlohmann::json;

namespace tket {

// The typedef PhasePolynomial leads to the python type Dict[List[bool], ...],
// which is not allowed at runtime because lists aren't hashable
typedef py::tket_custom::SequenceVec<std::pair<py::tket_custom::SequenceVec<bool>, Expr>> PyPhasePolynomialAlternate;
PhasePolynomial to_cpp_phase_poly(const PyPhasePolynomialAlternate& py_phase_poly){
    PhasePolynomial phase_poly;
    for (const auto& pair: py_phase_poly){
       phase_poly.insert_or_assign(pair.first, pair.second);
    }
    return phase_poly;
}
typedef std::map<py::tket_custom::TupleVec<bool>, Expr> PyPhasePolynomial;
PhasePolynomial to_cpp_phase_poly(const PyPhasePolynomial& py_phase_poly){
    return PhasePolynomial(std::make_move_iterator(py_phase_poly.begin()), std::make_move_iterator(py_phase_poly.end()));
}

// state_perm_t has the same hashability problem
typedef py::tket_custom::SequenceVec<std::pair<std::vector<bool>, py::tket_custom::SequenceVec<bool>>>
    py_state_perm_t;
    state_perm_t to_cpp_state_perm_t(const py_state_perm_t& py_state_perm){
        return state_perm_t(std::make_move_iterator(py_state_perm.begin()), std::make_move_iterator(py_state_perm.end()));
    }
typedef std::map<py::tket_custom::TupleVec<bool>, py::tket_custom::SequenceVec<bool>> py_state_perm_t2;
    state_perm_t to_cpp_state_perm_t(const py_state_perm_t2& py_state_perm){
        return state_perm_t(std::make_move_iterator(py_state_perm.begin()), std::make_move_iterator(py_state_perm.end()));
    }

typedef py::tket_custom::SequenceVec<std::pair<py::tket_custom::SequenceVec<bool>, Op_ptr>> py_ctrl_op_map_t_alt;
    ctrl_op_map_t to_cpp_ctrl_op_map_t(const py_ctrl_op_map_t_alt& ctrl_op_list) {
        ctrl_op_map_t ctrl_op_map;
        for (const auto &bit_op_pair: ctrl_op_list) {
            ctrl_op_map.insert_or_assign(bit_op_pair.first, bit_op_pair.second);
        }
        return ctrl_op_map;
    }
typedef std::map<py::tket_custom::TupleVec<bool>, Op_ptr> py_ctrl_op_map_t;
    ctrl_op_map_t to_cpp_ctrl_op_map_t(const py_ctrl_op_map_t& ctrl_op_map) {
        return ctrl_op_map_t(std::make_move_iterator(ctrl_op_map.begin()), std::make_move_iterator(ctrl_op_map.end()));
    }

typedef py::tket_custom::SequenceVec<std::pair<py::tket_custom::SequenceVec<bool>, py::tket_custom::SequenceVec<Op_ptr>>>
            py_ctrl_tensored_op_map_t_alt;
    ctrl_tensored_op_map_t to_cpp_ctrl_op_map_t(const py_ctrl_tensored_op_map_t_alt & ctrl_op_list) {
        ctrl_tensored_op_map_t ctrl_op_map;
        for (const auto &bit_op_pair: ctrl_op_list) {
            ctrl_op_map.insert_or_assign(bit_op_pair.first, bit_op_pair.second);
        }
        return ctrl_op_map;
    }
typedef std::map<py::tket_custom::TupleVec<bool>, py::tket_custom::SequenceVec<Op_ptr>>
            py_ctrl_tensored_op_map_t;
    ctrl_tensored_op_map_t to_cpp_ctrl_op_map_t(const py_ctrl_tensored_op_map_t & ctrl_op_map) {
        return ctrl_tensored_op_map_t(std::make_move_iterator(ctrl_op_map.begin()), std::make_move_iterator(ctrl_op_map.end()));
    }


    // Cast the std::vector keys in a map to py::tuple, since vector is not hashable
// in python
template <class T1, class T2>
std::map<py::tket_custom::TupleVec<T1>, T2> cast_keys_to_tuples(
    const std::map<std::vector<T1>, T2> &map) {
  std::map<py::tket_custom::TupleVec<T1>, T2> outmap;
  for (const auto &pair : map) {
      py::tket_custom::TupleVec<T1> tuple_vec(std::make_move_iterator(pair.first.begin()), std::make_move_iterator(pair.first.end()));
    outmap.emplace(tuple_vec, pair.second);
  }
  return outmap;
}
// Cast std::map<T1, T2> to std::vector<std::pair<T1, T2>>, useful when T1
// is a vector because vector, i.e. python list, is not hashable in python
template <class T1, class T2>
std::vector<std::pair<T1, T2>> cast_map_to_vector_of_pairs(
    const std::map<T1, T2> &map) {
  std::vector<std::pair<T1, T2>> result;
  result.reserve(map.size());
  for (const auto &pair : map) {
    result.emplace_back(pair.first, pair.second);
  }
  return result;
}

void init_boxes(py::module &m) {
  py::class_<CircBox, std::shared_ptr<CircBox>, Op>(
      m, "CircBox",
      "A user-defined operation specified by a :py:class:`Circuit`.")
      .def(
          py::init<const Circuit &>(), "Construct from a :py:class:`Circuit`.",
          py::arg("circ"))
      .def(
          "get_circuit", [](CircBox &cbox) { return *cbox.to_circuit(); },
          ":return: the :py:class:`Circuit` described by the box")
      .def(
          "symbol_substitution",
          [](CircBox &circ, const symbol_map_t &sub_map) {
            circ.symbol_substitution_in_place(sub_map);
          },
          "In-place substitution of symbolic expressions "
          "within underlying circuit; iterates "
          "through each parameterised gate/box within the circuit "
          "and performs the substitution. \n\n WARNING: "
          "This method potentially mutates the CircBox and "
          "any changes are propagated to "
          "any Circuit that the CircBox has been added to "
          "(via Circuit.add_circbox). \n\n:param symbol_map: "
          "A map from SymPy symbols to SymPy expressions",
          py::arg("symbol_map"));
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
          py::init<
              const std::vector<Pauli> &, const Expr &, const CXConfigType &>(),
          "Construct :math:`e^{-\\frac12 i \\pi t \\sigma_0 \\otimes "
          "\\sigma_1 \\otimes \\cdots}` from Pauli operators "
          ":math:`\\sigma_i \\in \\{I,X,Y,Z\\}` and a parameter "
          ":math:`t`.",
          py::arg("paulis"), py::arg("t"),
          py::arg("cx_config_type") = CXConfigType::Tree)
      .def(
          "get_circuit", [](PauliExpBox &pbox) { return *pbox.to_circuit(); },
          ":return: the :py:class:`Circuit` described by the box")
      .def(
          "get_paulis", &PauliExpBox::get_paulis,
          ":return: the corresponding list of " CLSOBJS(Pauli))
      .def(
          "get_phase", &PauliExpBox::get_phase,
          ":return: the corresponding phase parameter")
      .def(
          "get_cx_config", &PauliExpBox::get_cx_config,
          ":return: decomposition method");
  py::class_<PauliExpPairBox, std::shared_ptr<PauliExpPairBox>, Op>(
      m, "PauliExpPairBox",
      "An operation defined as a pair of exponentials of a tensor of Pauli "
      "operations and their (possibly symbolic) phase parameters.")
      .def(
          py::init<
              const std::vector<Pauli> &, const Expr &,
              const std::vector<Pauli> &, Expr, CXConfigType>(),
          "Construct a pair of Pauli exponentials of the form"
          " :math:`e^{-\\frac12 i \\pi t_j \\sigma_0 \\otimes "
          "\\sigma_1 \\otimes \\cdots}` from Pauli operator strings "
          ":math:`\\sigma_i \\in \\{I,X,Y,Z\\}` and parameters "
          ":math:`t_j, j \\in \\{0,1\\}`.",
          py::arg("paulis0"), py::arg("t0"), py::arg("paulis1"), py::arg("t1"),
          py::arg("cx_config_type") = CXConfigType::Tree)
      .def(
          "get_circuit",
          [](PauliExpPairBox &pbox) { return *pbox.to_circuit(); },
          ":return: the :py:class:`Circuit` described by the box")
      .def(
          "get_paulis_pair", &PauliExpPairBox::get_paulis_pair,
          ":return: A tuple containing the two corresponding lists of " CLSOBJS(
              Pauli))
      .def(
          "get_phase_pair", &PauliExpPairBox::get_phase_pair,
          ":return: A tuple containing the two phase parameters")
      .def(
          "get_cx_config", &PauliExpPairBox::get_cx_config,
          ":return: decomposition method");
  py::class_<
      PauliExpCommutingSetBox, std::shared_ptr<PauliExpCommutingSetBox>, Op>(
      m, "PauliExpCommutingSetBox",
      "An operation defined as a set of commuting of exponentials of a"
      "tensor of Pauli operations and their (possibly symbolic) phase "
      "parameters.")
      .def(
          py::init<
              const std::vector<std::pair<std::vector<Pauli>, Expr>>,
              const CXConfigType &>(),
          "Construct a set of necessarily commuting Pauli exponentials of the "
          "form"
          " :math:`e^{-\\frac12 i \\pi t_j \\sigma_0 \\otimes "
          "\\sigma_1 \\otimes \\cdots}` from Pauli operator strings "
          ":math:`\\sigma_i \\in \\{I,X,Y,Z\\}` and parameters "
          ":math:`t_j, j \\in \\{0, 1, \\cdots \\}`.",
          py::arg("pauli_gadgets"),
          py::arg("cx_config_type") = CXConfigType::Tree)
      .def(
          "get_circuit",
          [](PauliExpCommutingSetBox &pbox) { return *pbox.to_circuit(); },
          ":return: the :py:class:`Circuit` described by the box")
      .def(
          "get_paulis", &PauliExpCommutingSetBox::get_pauli_gadgets,
          ":return: the corresponding list of Pauli gadgets")
      .def(
          "get_cx_config", &PauliExpCommutingSetBox::get_cx_config,
          ":return: decomposition method");
  py::enum_<ToffoliBoxSynthStrat>(
      m, "ToffoliBoxSynthStrat",
      "Enum strategies for synthesising ToffoliBoxes")
      .value(
          "Matching", ToffoliBoxSynthStrat::Matching,
          "Use multiplexors to perform parallel swaps on hypercubes")
      .value(
          "Cycle", ToffoliBoxSynthStrat::Cycle,
          "Use CnX gates to perform transpositions");
  py::class_<ToffoliBox, std::shared_ptr<ToffoliBox>, Op>(
      m, "ToffoliBox",
      "An operation that constructs a circuit to implement the specified "
      "permutation of classical basis states.")
      .def(
          py::init([](const py_state_perm_t &perm, ToffoliBoxSynthStrat strat,
                      OpType optype) {
            return ToffoliBox(to_cpp_state_perm_t(perm), strat, optype);
          }),
          "Construct from a permutation of basis states\n\n"
          ":param permutation: a list of bitstring pairs\n"
          ":param strat: synthesis strategy\n"
          ":param rotation_axis: the rotation axis of the multiplexors used in "
          "the decomposition. Can be either Rx or Ry. Only applicable to the "
          "Matching strategy. Default to Ry.",
          py::arg("permutation"), py::arg("strat"),
          py::arg("rotation_axis") = OpType::Ry)
      .def(
          py::init(
              [](const py_state_perm_t &perm, const OpType &rotation_axis) {
                return ToffoliBox(to_cpp_state_perm_t(perm), ToffoliBoxSynthStrat::Matching, rotation_axis);
              }),
          "Construct from a permutation of basis states and perform synthesis "
          "using the Matching strategy\n\n"
          ":param permutation: a list of bitstring pairs\n"
          ":param rotation_axis: the rotation axis of the multiplexors used in "
          "the decomposition. Can be either Rx or Ry, default to Ry.",
          py::arg("permutation"), py::arg("rotation_axis") = OpType::Ry)
      .def(
          py::init([](unsigned n_qubits, const py_state_perm_t &perm,
                      const OpType &rotation_axis) {
            (void)n_qubits;
            PyErr_WarnEx(
                PyExc_DeprecationWarning,
                "The argument n_qubits is no longer needed. "
                "Please create ToffoliBoxes without n_qubits.",
                1);
            return ToffoliBox(
                    to_cpp_state_perm_t(perm), ToffoliBoxSynthStrat::Matching, rotation_axis);
          }),
          "Constructor for backward compatibility. Subject to deprecation.",
          py::arg("n_qubits"), py::arg("permutation"),
          py::arg("rotation_axis") = OpType::Ry)
      .def( py::init([](const py_state_perm_t2 &perm, ToffoliBoxSynthStrat strat,
                              OpType optype) {
                      return ToffoliBox(to_cpp_state_perm_t(perm), strat, optype);
                  }),
          "Construct from a permutation of basis states\n\n"
          ":param permutation: a map between bitstrings\n"
          ":param strat: synthesis strategy\n"
          ":param rotation_axis: the rotation axis of the multiplexors used in "
          "the decomposition. Can be either Rx or Ry. Only applicable to the "
          "Matching strategy. Default to Ry.",
          py::arg("permutation"), py::arg("strat"),
          py::arg("rotation_axis") = OpType::Ry)
      .def(
          py::init([](const py_state_perm_t2 &perm, const OpType &rotation_axis) {
            return ToffoliBox(
                    to_cpp_state_perm_t(perm), ToffoliBoxSynthStrat::Matching, rotation_axis);
          }),
          "Construct from a permutation of basis states and perform synthesis "
          "using the Matching strategy\n\n"
          ":param permutation: a map between bitstrings\n"
          ":param rotation_axis: the rotation axis of the multiplexors used in "
          "the decomposition. Can be either Rx or Ry, default to Ry.",
          py::arg("permutation"), py::arg("rotation_axis") = OpType::Ry)
      .def(
          py::init([](unsigned n_qubits, const py_state_perm_t2 &perm,
                      const OpType &rotation_axis) {
            (void)n_qubits;
            PyErr_WarnEx(
                PyExc_DeprecationWarning,
                "The argument n_qubits is no longer needed. "
                "Please create ToffoliBoxes without n_qubits.",
                1);
            return ToffoliBox(
                    to_cpp_state_perm_t(perm), ToffoliBoxSynthStrat::Matching, rotation_axis);
          }),
          "Constructor for backward compatibility. Subject to deprecation.",
          py::arg("n_qubits"), py::arg("permutation"),
          py::arg("rotation_axis") = OpType::Ry)
      .def(
          "get_circuit", [](ToffoliBox &tbox) { return *tbox.to_circuit(); },
          ":return: the :py:class:`Circuit` described by the box")
      .def(
          "get_permutation",
          [](ToffoliBox &box) {
            std::map<py::tuple, py::tuple> outmap;
            for (const auto &pair : box.get_permutation()) {
              outmap.insert(
                  {py::tuple(py::cast(pair.first)),
                   py::tuple(py::cast(pair.second))});
            }
            return outmap;
          },
          ":return: the permutation")
      .def(
          "get_strat", &ToffoliBox::get_strat,
          ":return: the synthesis strategy")
      .def(
          "get_rotation_axis", &ToffoliBox::get_rotation_axis,
          ":return: the rotation axis");
  py::class_<QControlBox, std::shared_ptr<QControlBox>, Op>(
      m, "QControlBox",
      "A user-defined controlled operation specified by an "
      ":py:class:`Op`, the number of quantum controls, and the control state "
      "expressed as an integer or a bit vector.")
      .def(
          py::init<Op_ptr &, unsigned, std::vector<bool> &>(),
          "Construct from an :py:class:`Op`, a number of quantum "
          "controls, and the control state expressed as a bit vector. The "
          "controls occupy the low-index ports of the "
          "resulting operation.\n\n"
          ":param op: the underlying operator\n"
          ":param n_controls: the number of control qubits. Default to 1\n"
          ":param control_state: the control state expressed as a bit vector. "
          "Default to all 1s\n",
          py::arg("op"), py::arg("n_controls") = 1,
          py::arg("control_state") = std::vector<bool>())
      .def(
          py::init([](Op_ptr &op, unsigned n_controls,
                      unsigned long long control_state) {
            return QControlBox(
                op, n_controls, dec_to_bin(control_state, n_controls));
          }),
          "Construct from an :py:class:`Op`, a number of quantum "
          "controls, and the control state expressed as an integer. The "
          "controls occupy the low-index ports of the "
          "resulting operation.\n\n"
          ":param op: the underlying operator\n"
          ":param n_controls: the number of control qubits\n"
          ":param control_state: the control state expressed as an integer. "
          "Big-endian\n",
          py::arg("op"), py::arg("n_controls"), py::arg("control_state"))
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
          ":return: the number of control qubits")
      .def(
          "get_control_state",
          [](QControlBox &qcbox) {
            return bin_to_dec(qcbox.get_control_state());
          },
          ":return: the control state as an integer (big-endian binary "
          "representation)")
      .def(
          "get_control_state_bits",
          [](QControlBox &qcbox) { return qcbox.get_control_state(); },
          ":return: the control state as a bit vector");

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
            return py::object(json(std::make_shared<CompositeGateDef>(c)))
                .

                cast<py::dict>();
          },
          ":return: a JSON serializable dictionary representation of "
          "the CustomGateDef")
      .def_static(
          "from_dict",
          [](const py::dict &composite_gate_def_dict) {
            return json(composite_gate_def_dict).get<composite_def_ptr_t>();
          },
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
                      const PyPhasePolynomial &p_p, const MatrixXb &lin_trans) {
            boost::bimap<Qubit, unsigned> bmap;
            for (const auto &pair : q_ind) {
              bmap.insert({pair.first, pair.second});
            }
            return PhasePolyBox(n_qb, bmap, to_cpp_phase_poly(p_p), lin_trans);
          }),
          "\n\nConstruct from the number of qubits, the mapping from "
          "Qubit to index, the phase polynomial (map from bitstring "
          "to phase) and the linear transformation (boolean matrix)",
          py::arg("n_qubits"), py::arg("qubit_indices"),
          py::arg("phase_polynomial"), py::arg("linear_transformation"))
      .def(
          py::init([](unsigned n_qb, const std::map<Qubit, unsigned> &q_ind,
                      const PyPhasePolynomialAlternate &p_p,
                      const MatrixXb &lin_trans) {
            boost::bimap<Qubit, unsigned> bmap;
            for (const auto &pair : q_ind) {
              bmap.insert({pair.first, pair.second});
            }
            return PhasePolyBox(n_qb, bmap, to_cpp_phase_poly(p_p), lin_trans);
          }),
          "Construct from the number of qubits, the mapping from "
          "Qubit to index, the phase polynomial (list of bitstring "
          "phase pairs) and the linear transformation (boolean matrix)\n\n"
          "If any bitstring is repeated in the phase polynomial list, the last "
          "given value for that bistring will be used",
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
          "phase_polynomial_as_list",
          [](PhasePolyBox &ppoly) {
            const PhasePolynomial &phase_pol = ppoly.get_phase_polynomial();
            return cast_map_to_vector_of_pairs(phase_pol);
          },
          "List of bitstring(basis state)-phase pairs.")
      .def_property_readonly(
          "phase_polynomial",
          [](PhasePolyBox &ppoly) {
            const PhasePolynomial &phase_pol = ppoly.get_phase_polynomial();
            return cast_keys_to_tuples(phase_pol);
          },
          "Map from bitstring (basis state) to phase.")
      .def_property_readonly(
          "phase_polynomial_as_list",
          [](PhasePolyBox &ppoly) {
            const PhasePolynomial &phase_pol = ppoly.get_phase_polynomial();
            return cast_map_to_vector_of_pairs(phase_pol);
          },
          "List of bitstring(basis state)-phase pairs.")
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
      "A user-defined assertion specified by a list of Pauli stabilisers.")
      .def(
          py::init<const PauliStabiliserList>(),
          "Construct from a list of Pauli stabilisers.\n\n"
          ":param m: The list of Pauli stabilisers\n",
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
                          "Invalid Pauli string: " + raw_string);
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
                        "Invalid Pauli string: " + raw_string);
                }
              }
              stabilisers.emplace_back(string, coeff);
            }
            return StabiliserAssertionBox(stabilisers);
          }),
          "Construct from a list of Pauli stabilisers.\n\n"
          ":param m: The list of Pauli stabilisers expressed as Python "
          "strings\n",
          py::arg("stabilisers"))
      .def(
          "get_circuit",
          [](StabiliserAssertionBox &ubox) { return *ubox.to_circuit(); },
          ":return: the :py:class:`Circuit` described by the box")
      .def(
          "get_stabilisers", &StabiliserAssertionBox::get_stabilisers,
          ":return: the list of Pauli stabilisers");


  py::class_<MultiplexorBox, std::shared_ptr<MultiplexorBox>, Op>(
      m, "MultiplexorBox",
      "A user-defined multiplexor (i.e. uniformly controlled operations) "
      "specified by a "
      "map from bitstrings to " CLSOBJS(Op)
      "or a list of bitstring-" CLSOBJS(Op) " pairs")
      .def(
              py::init([](const py_ctrl_op_map_t_alt & bitstring_op_pairs){
                  return MultiplexorBox(to_cpp_ctrl_op_map_t(bitstring_op_pairs));}),
      "Construct from a list of bitstring-" CLSOBJS(Op) "pairs\n\n"
      ":param bitstring_to_op_list: List of bitstring-" CLSOBJS(Op) "pairs\n",
      py::arg("bistring_to_op_list"))
      .def(
              py::init([](const py_ctrl_op_map_t & bitstring_op_map){
                  return MultiplexorBox(to_cpp_ctrl_op_map_t(bitstring_op_map));}),
          "Construct from a map from bitstrings to " CLSOBJS(Op) "\n\n"
          ":param op_map: Map from bitstrings to " CLSOBJS(Op) "\n",
          py::arg("op_map"))
      .def(
          "get_circuit", [](MultiplexorBox &box) { return *box.to_circuit(); },
          ":return: the :py:class:`Circuit` described by the box")
      .def(
          "get_op_map",
          [](MultiplexorBox &box) {
            return cast_keys_to_tuples(box.get_op_map());
          },
          ":return: the underlying op map")
      .def(
          "get_bitstring_op_pair_list",
          [](MultiplexorBox &box) {
          return cast_map_to_vector_of_pairs(box.get_op_map());
          },
          ":return: the underlying bistring-op pairs");
  py::class_<
      MultiplexedRotationBox, std::shared_ptr<MultiplexedRotationBox>, Op>(
      m, "MultiplexedRotationBox",
      "A user-defined multiplexed rotation gate (i.e. "
      "uniformly controlled single-axis rotations) specified by "
      "a map from bitstrings to " CLSOBJS(Op)
      "or a list of bitstring-" CLSOBJS(Op) " pairs")
      .def(
              py::init([](const py_ctrl_op_map_t_alt & bitstring_op_pairs){
          return MultiplexedRotationBox(to_cpp_ctrl_op_map_t(bitstring_op_pairs));}),
          "Construct from a list of bitstring-" CLSOBJS(Op) "pairs\n\n"
          "All " CLSOBJS(Op) "  must share the same single-qubit rotation type: "
          "Rx, Ry, or Rz.\n\n"
          ":param bitstring_to_op_list: List of bitstring-" CLSOBJS(Op) "pairs\n",
          py::arg("bistring_to_op_list"))
      .def(
              py::init([](const py_ctrl_op_map_t & bitstring_op_map){
                  return MultiplexedRotationBox(to_cpp_ctrl_op_map_t(bitstring_op_map));}),
          "Construct from a map from bitstrings to :py:class:`Op` s."
          "All " CLSOBJS(Op) "  must share the same single-qubit rotation type: "
          "Rx, Ry, or Rz.\n\n"
          ":param op_map: Map from bitstrings to " CLSOBJS(Op) "\n",
          py::arg("op_map"))
      .def(
          py::init([](const std::vector<double> &angles, const OpType &axis) {
            if (angles.size() == 0) {
              throw std::invalid_argument("Angles are empty.");
            }
            if (angles.size() & (angles.size() - 1)) {
              throw std::invalid_argument(
                  "The size of the angles is not power of 2.");
            }
            if (axis != OpType::Rx && axis != OpType::Ry &&
                axis != OpType::Rz) {
              throw std::invalid_argument(
                  "The axis must be either Rx, Ry, or Rz.");
            }
            unsigned bitstring_width = (unsigned)log2(angles.size());
            ctrl_op_map_t op_map;
            for (unsigned i = 0; i < angles.size(); i++) {
              if (std::abs(angles[i]) > EPS) {
                std::vector<bool> bits = dec_to_bin(i, bitstring_width);
                op_map.insert({bits, get_op_ptr(axis, angles[i])});
              }
            }
            return MultiplexedRotationBox(op_map);
          }),
          "Construct from a list of angles and the rotation axis.\n\n"
          ":param angles: List of rotation angles in half-turns. angles[i] is "
          "the angle activated by the binary representation of i\n"
          ":param axis: ``OpType.Rx``, ``OpType.Ry`` or ``OpType.Rz``\n",
          py::arg("angles"), py::arg("axis"))
      .def(
          "get_circuit",
          [](MultiplexedRotationBox &box) { return *box.to_circuit(); },
          ":return: the :py:class:`Circuit` described by the box")
      .def(
          "get_bitstring_op_pair_list",
          [](MultiplexedRotationBox &box) {
            return cast_map_to_vector_of_pairs(box.get_op_map());
          },
          ":return: the underlying bistring-op pairs")
      .def(
          "get_op_map",
          [](MultiplexedRotationBox &box) {
            return cast_keys_to_tuples(box.get_op_map());
          },
          ":return: the underlying op map");
  py::class_<MultiplexedU2Box, std::shared_ptr<MultiplexedU2Box>, Op>(
      m, "MultiplexedU2Box",
      "A user-defined multiplexed U2 gate (i.e. uniformly controlled U2 "
      "gate) specified by a "
      "map from bitstrings to " CLSOBJS(Op)
      "or a list of bitstring-" CLSOBJS(Op) " pairs")
      .def(
              py::init([](const py_ctrl_op_map_t_alt & bitstring_op_pairs, bool impl_diag){
          return MultiplexedU2Box(to_cpp_ctrl_op_map_t(bitstring_op_pairs), impl_diag);}),
      "Construct from a list of bitstring-" CLSOBJS(Op) "pairs\n\n"
      "Only supports single qubit unitary gate types and "
      ":py:class:`Unitary1qBox`.\n\n"
      ":param op_map: List of bitstring-" CLSOBJS(Op) "pairs\n"
      ":param impl_diag: Whether to implement the final diagonal gate, "
      "default to True.",
      py::arg("bistring_to_op_list"), py::arg("impl_diag") = true)
      .def(
              py::init([](const py_ctrl_op_map_t & bitstring_op_map, bool impl_diag){
                  return MultiplexedU2Box(to_cpp_ctrl_op_map_t(bitstring_op_map), impl_diag);}),
          "Construct from a map from bitstrings to " CLSOBJS(Op) "."
          "Only supports single qubit unitary gate types and "
          ":py:class:`Unitary1qBox`.\n\n"
          ":param op_map: Map from bitstrings to " CLSOBJS(Op) "\n"
          ":param impl_diag: Whether to implement the final diagonal gate, "
          "default to True.",
          py::arg("op_map"), py::arg("impl_diag") = true)
      .def(
          "get_circuit",
          [](MultiplexedU2Box &box) { return *box.to_circuit(); },
          ":return: the :py:class:`Circuit` described by the box")
      .def(
          "get_bitstring_op_pair_list",
          [](MultiplexedU2Box &box) {
          return cast_map_to_vector_of_pairs(box.get_op_map());
          },
          ":return: the underlying bistring-op pairs")
      .def(
          "get_op_map",
          [](MultiplexedU2Box &box) {
            return cast_keys_to_tuples(box.get_op_map());
          },
          ":return: the underlying op map")
      .def(
          "get_impl_diag", &MultiplexedU2Box::get_impl_diag,
          ":return: flag indicating whether to implement the final diagonal "
          "gate.");

  py::class_<
      MultiplexedTensoredU2Box, std::shared_ptr<MultiplexedTensoredU2Box>, Op>(
      m, "MultiplexedTensoredU2Box",
      "A user-defined multiplexed tensor product of U2 gates specified by a "
      "map from bitstrings to lists of " CLSOBJS(Op)
      "or a list of bitstring-list(" CLSOBJS(Op) ") pairs")
      .def(
              py::init([](const py_ctrl_tensored_op_map_t_alt & bitstring_op_pairs){
          return MultiplexedTensoredU2Box(to_cpp_ctrl_op_map_t(bitstring_op_pairs));}),
          "Construct from a list of bitstring-" CLSOBJS(Op) "pairs\n\n"
          "Only supports single qubit unitary gate types and "
          ":py:class:`Unitary1qBox`.\n\n"
          ":param bitstring_to_op_list: List of bitstring-List of " CLSOBJS(Op) " pairs\n",
          py::arg("bistring_to_op_list"))
      .def(
              py::init([](const py_ctrl_tensored_op_map_t & bitstring_op_map){
                  return MultiplexedTensoredU2Box(to_cpp_ctrl_op_map_t(bitstring_op_map));}),
          "Construct from a map from bitstrings to equal-sized lists of "
          CLSOBJS(Op) ". "
          "Only supports single qubit unitary gate types and "
          ":py:class:`Unitary1qBox`.\n\n"
          ":param op_map: Map from bitstrings to lists of " CLSOBJS(Op),
          py::arg("op_map"))
      .def(
          "get_circuit",
          [](MultiplexedTensoredU2Box &box) { return *box.to_circuit(); },
          ":return: the :py:class:`Circuit` described by the box")
      .def(
          "get_bitstring_op_pair_list",
          [](MultiplexedTensoredU2Box &box) {
          return cast_map_to_vector_of_pairs(box.get_op_map());
          },
          ":return: the underlying bistring-op pairs")
      .def(
          "get_op_map",
          [](MultiplexedTensoredU2Box &box) {
            return cast_keys_to_tuples(box.get_op_map());
          },
          ":return: the underlying op map");
  py::class_<StatePreparationBox, std::shared_ptr<StatePreparationBox>, Op>(
      m, "StatePreparationBox",
      "A box for preparing quantum states using multiplexed-Ry and "
      "multiplexed-Rz gates")
      .def(
          py::init<const Eigen::VectorXcd &, bool, bool>(),
          "Construct from a statevector\n\n"
          ":param statevector: normalised statevector\n"
          ":param is_inverse: whether to implement the dagger of the state "
          "preparation circuit, default to false\n"
          ":param with_initial_reset: whether to explicitly set the state to "
          "zero initially (by default the initial zero state is assumed and no "
          "explicit reset is applied)",
          py::arg("statevector"), py::arg("is_inverse") = false,
          py::arg("with_initial_reset") = false)
      .def(
          "get_circuit",
          [](StatePreparationBox &box) { return *box.to_circuit(); },
          ":return: the :py:class:`Circuit` described by the box")
      .def(
          "get_statevector", &StatePreparationBox::get_statevector,
          ":return: the statevector")
      .def(
          "is_inverse", &StatePreparationBox::is_inverse,
          ":return: flag indicating whether to implement the dagger of the "
          "state preparation circuit")
      .def(
          "with_initial_reset", &StatePreparationBox::with_initial_reset,
          ":return: flag indicating whether the qubits are explicitly "
          "set to the zero state initially");
  py::class_<DiagonalBox, std::shared_ptr<DiagonalBox>, Op>(
      m, "DiagonalBox",
      "A box for synthesising a diagonal unitary matrix into a sequence of "
      "multiplexed-Rz gates.")
      .def(
          py::init<const Eigen::VectorXcd &, bool>(),
          "Construct from the diagonal entries of the unitary operator. The "
          "size of the vector must be 2^n where n is a positive integer.\n\n"
          ":param diagonal: diagonal entries\n"
          ":param upper_triangle: indicates whether the multiplexed-Rz gates "
          "take the shape of an upper triangle or a lower triangle. Default to "
          "true.",
          py::arg("diagonal"), py::arg("upper_triangle") = true)
      .def(
          "get_circuit", [](DiagonalBox &box) { return *box.to_circuit(); },
          ":return: the :py:class:`Circuit` described by the box")
      .def(
          "get_diagonal", &DiagonalBox::get_diagonal,
          ":return: the statevector")
      .def(
          "is_upper_triangle", &DiagonalBox::is_upper_triangle,
          ":return: the upper_triangle flag");
  py::class_<ConjugationBox, std::shared_ptr<ConjugationBox>, Op>(
      m, "ConjugationBox",
      "A box to express computations that follow the compute-action-uncompute "
      "pattern.")
      .def(
          py::init<
              const Op_ptr &, const Op_ptr &, const std::optional<Op_ptr>>(),
          "Construct from operations that perform compute, action, and "
          "uncompute. All three operations need to be quantum and have the "
          "same size.\n\n"
          ":param compute: the compute operation\n"
          ":param action: the action operation\n"
          ":param uncompute: optional uncompute operation, default to "
          "compute.dagger(). If provided, the user needs to make sure that "
          "uncompute.dagger() and compute have the same unitary.",
          py::arg("compute"), py::arg("action"),
          py::arg("uncompute") = std::nullopt)
      .def(
          "get_circuit", [](ConjugationBox &box) { return *box.to_circuit(); },
          ":return: the :py:class:`Circuit` described by the box")
      .def(
          "get_compute", &ConjugationBox::get_compute,
          ":return: the compute operation")
      .def(
          "get_action", &ConjugationBox::get_action,
          ":return: the action operation")
      .def(
          "get_uncompute", &ConjugationBox::get_uncompute,
          ":return: the uncompute operation. Returns None if the default "
          "compute.dagger() is used");
}
}  // namespace tket
