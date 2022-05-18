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

#include <sstream>

#include "Circuit/Command.hpp"
#include "Gate/Gate.hpp"
#include "Gate/OpPtrFunctions.hpp"
#include "Gate/SymTable.hpp"
#include "Ops/ClassicalOps.hpp"
#include "Ops/Op.hpp"
#include "Utils/Constants.hpp"
#include "Utils/Symbols.hpp"
#include "binder_json.hpp"
#include "binder_utils.hpp"
#include "typecast.hpp"

namespace py = pybind11;
using json = nlohmann::json;

namespace tket {

void init_unitid(py::module &m);
void init_circuit(py::module &m);
void init_classical(py::module &m);
void init_boxes(py::module &m);
void init_library(py::module &m);

PYBIND11_MODULE(circuit, m) {
  init_unitid(m);
  py::class_<Op, std::shared_ptr<Op>>(
      m, "Op", "Encapsulates operation information")
      .def_static(
          "create",
          [](OpType optype) { return get_op_ptr(optype, std::vector<Expr>()); },
          "Create an :py:class:`Op` with given type")
      .def_static(
          "create",
          [](OpType optype, Expr param) { return get_op_ptr(optype, param); },
          "Create an :py:class:`Op` with given type and parameter")
      .def_static(
          "create",
          [](OpType optype, const std::vector<Expr> &params) {
            return get_op_ptr(optype, params);
          },
          "Create an :py:class:`Op` with given type and parameters")
      .def_property_readonly(
          "type", &Op::get_type, "Type of op being performed")
      .def_property_readonly(
          "params", &Op::get_params_reduced,
          "Angular parameters of the op, in half-turns (e.g. 1.0 "
          "half-turns is :math:`\\pi` radians). The parameters "
          "returned are constrained to the appropriate canonical "
          "range, which is usually the half-open interval [0,2) but "
          "for some operations (e.g. Rx, Ry and Rz) is [0,4).")
      .def_property_readonly(
          "n_qubits", &Op::n_qubits, "Number of qubits of op")
      .def_property_readonly("dagger", &Op::dagger, "Dagger of op")
      .def_property_readonly("transpose", &Op::transpose, "Transpose of op")
      .def(
          "get_name", &Op::get_name, "String representation of op",
          py::arg("latex") = false)
      .def("__eq__", &Op::operator==)
      .def("__repr__", [](const Op &op) { return op.get_name(); })
      .def("free_symbols", [](const Op &op) { return op.free_symbols(); })
      .def(
          "get_unitary",
          [](const Op *op) {
            const auto &gate = static_cast<const Gate &>(*op);
            return gate.get_unitary();
          })
      .def(
          "is_clifford_type",
          [](const Op &op) { return op.get_desc().is_clifford_gate(); })
      .def("is_gate", [](const Op &op) { return op.get_desc().is_gate(); });

  // NOTE: Sphinx does not automatically pick up the docstring for OpType
  py::enum_<OpType>(
      m, "OpType",
      "Enum for available operations compatible with "
      "tket " CLSOBJS(Circuit) ".",
      py::arithmetic())
      .value(
          "Z", OpType::Z,
          "Pauli Z: :math:`\\left[ \\begin{array}{cc} 1 & 0 \\\\ 0 & "
          "-1 \\end{array} \\right]`")
      .value(
          "X", OpType::X,
          "Pauli X: :math:`\\left[ \\begin{array}{cc} 0 & 1 \\\\ 1 & "
          "0 \\end{array} \\right]`")
      .value(
          "Y", OpType::Y,
          "Pauli Y: :math:`\\left[ \\begin{array}{cc} 0 & -i \\\\ i & "
          "0 \\end{array} \\right]`")
      .value(
          "S", OpType::S,
          ":math:`\\left[ \\begin{array}{cc} 1 & 0 \\\\ 0 & i "
          "\\end{array} \\right] = \\mathrm{U1}(\\frac12)`")
      .value(
          "Sdg", OpType::Sdg,
          ":math:`\\mathrm{S}^{\\dagger} = \\left[ \\begin{array}{cc} "
          "1 & 0 \\\\ 0 & -i \\end{array} \\right] = "
          "\\mathrm{U1}(-\\frac12)`")
      .value(
          "T", OpType::T,
          ":math:`\\left[ \\begin{array}{cc} 1 & 0 \\\\ 0 & "
          "e^{i\\pi/4} \\end{array} \\right] = "
          "\\mathrm{U1}(\\frac14)`")
      .value(
          "Tdg", OpType::Tdg,
          ":math:`\\mathrm{T}^{\\dagger} = \\left[ \\begin{array}{cc} "
          "1 & 0 \\\\ 0 & e^{-i\\pi/4} \\end{array} \\right] = "
          "\\mathrm{U1}(-\\frac14)`")
      .value(
          "V", OpType::V,
          ":math:`\\frac{1}{\\sqrt 2} \\left[ \\begin{array}{cc} 1 & "
          "-i \\\\ -i & 1 \\end{array} \\right] = "
          "\\mathrm{Rx}(\\frac12)`")
      .value(
          "Vdg", OpType::Vdg,
          ":math:`\\mathrm{V}^{\\dagger} = \\frac{1}{\\sqrt 2} "
          "\\left[ \\begin{array}{cc} 1 & i \\\\ i & 1 \\end{array} "
          "\\right] = \\mathrm{Rx}(-\\frac12)`")
      .value(
          "SX", OpType::SX,
          ":math:`\\frac{1}{2} \\left[ \\begin{array}{cc} 1 + i & "
          "1 - i \\\\ 1 - i & 1 + i \\end{array} \\right] = "
          "e^{\\frac{i\\pi}{4}}\\mathrm{Rx}(\\frac12)`")
      .value(
          "SXdg", OpType::SXdg,
          ":math:`\\mathrm{SX}^{\\dagger} = \\frac{1}{2} "
          "\\left[ \\begin{array}{cc} 1 - i & 1 + i "
          "\\\\ 1 + i & 1 - i \\end{array} "
          "\\right] = e^{\\frac{-i\\pi}{4}}\\mathrm{Rx}(-\\frac12)`")
      .value(
          "H", OpType::H,
          "Hadamard gate: :math:`\\frac{1}{\\sqrt 2} \\left[ "
          "\\begin{array}{cc} 1 & 1 \\\\ 1 & -1 \\end{array} "
          "\\right]`")
      .value(
          "Rx", OpType::Rx,
          ":math:`(\\alpha) \\mapsto e^{-\\frac12 i \\pi \\alpha "
          "\\mathrm{X}} = \\left[ \\begin{array}{cc} "
          "\\cos\\frac{\\pi\\alpha}{2} & "
          "-i\\sin\\frac{\\pi\\alpha}{2} \\\\ "
          "-i\\sin\\frac{\\pi\\alpha}{2} & "
          "\\cos\\frac{\\pi\\alpha}{2} \\end{array} \\right]`")
      .value(
          "Ry", OpType::Ry,
          ":math:`(\\alpha) \\mapsto e^{-\\frac12 i \\pi \\alpha "
          "\\mathrm{Y}} = \\left[ \\begin{array}{cc} "
          "\\cos\\frac{\\pi\\alpha}{2} & -\\sin\\frac{\\pi\\alpha}{2} "
          "\\\\ \\sin\\frac{\\pi\\alpha}{2} & "
          "\\cos\\frac{\\pi\\alpha}{2} \\end{array} \\right]`")
      .value(
          "Rz", OpType::Rz,
          ":math:`(\\alpha) \\mapsto e^{-\\frac12 i \\pi \\alpha "
          "\\mathrm{Z}} = \\left[ \\begin{array}{cc} e^{-\\frac12 i "
          "\\pi\\alpha} & 0 \\\\ 0 & e^{\\frac12 i \\pi\\alpha} "
          "\\end{array} \\right]`")
      .value(
          "U1", OpType::U1,
          ":math:`(\\lambda) \\mapsto \\mathrm{U3}(0, 0, \\lambda) = "
          "e^{\\frac12 i\\pi\\lambda} \\mathrm{Rz}(\\lambda)`. "
          "U-gates are used by IBM. See "
          "https://qiskit.org/documentation/tutorials/circuits/"
          "3_summary_of_quantum_operations.html "
          "for more information on U-gates.")
      .value(
          "U2", OpType::U2,
          ":math:`(\\phi, \\lambda) \\mapsto \\mathrm{U3}(\\frac12, "
          "\\phi, \\lambda) = e^{\\frac12 i\\pi(\\lambda+\\phi)} "
          "\\mathrm{Rz}(\\phi) \\mathrm{Ry}(\\frac12) "
          "\\mathrm{Rz}(\\lambda)`, defined by matrix multiplication")
      .value(
          "U3", OpType::U3,
          ":math:`(\\theta, \\phi, \\lambda) \\mapsto  \\left[ "
          "\\begin{array}{cc} \\cos\\frac{\\pi\\theta}{2} & "
          "-e^{i\\pi\\lambda} \\sin\\frac{\\pi\\theta}{2} \\\\ "
          "e^{i\\pi\\phi} \\sin\\frac{\\pi\\theta}{2} & "
          "e^{i\\pi(\\lambda+\\phi)} \\cos\\frac{\\pi\\theta}{2} "
          "\\end{array} \\right] = e^{\\frac12 i\\pi(\\lambda+\\phi)} "
          "\\mathrm{Rz}(\\phi) \\mathrm{Ry}(\\theta) "
          "\\mathrm{Rz}(\\lambda)`")
      .value(
          "TK1", OpType::TK1,
          ":math:`(\\alpha, \\beta, \\gamma) \\mapsto "
          "\\mathrm{Rz}(\\alpha) \\mathrm{Rx}(\\beta) "
          "\\mathrm{Rz}(\\gamma)`")
      .value(
          "TK2", OpType::TK2,
          ":math:`(\\alpha, \\beta, \\gamma) \\mapsto "
          "\\mathrm{XXPhase}(\\alpha) "
          "\\mathrm{YYPhase}(\\beta) "
          "\\mathrm{ZZPhase}(\\gamma)`")
      .value("CX", OpType::CX, "Controlled :math:`\\mathrm{X}` gate")
      .value("CY", OpType::CY, "Controlled :math:`\\mathrm{Y}` gate")
      .value("CZ", OpType::CZ, "Controlled :math:`\\mathrm{Z}` gate")
      .value("CH", OpType::CH, "Controlled :math:`\\mathrm{H}` gate")
      .value("CV", OpType::CV, "Controlled :math:`\\mathrm{V}` gate")
      .value(
          "CVdg", OpType::CVdg,
          "Controlled :math:`\\mathrm{V}^{\\dagger}` gate")
      .value("CSX", OpType::CSX, "Controlled :math:`\\mathrm{SX}` gate")
      .value(
          "CSXdg", OpType::CSXdg,
          "Controlled :math:`\\mathrm{SX}^{\\dagger}` gate")
      .value(
          "CRz", OpType::CRz,
          ":math:`(\\alpha) \\mapsto` Controlled "
          ":math:`\\mathrm{Rz}(\\alpha)` gate")
      .value(
          "CRx", OpType::CRx,
          ":math:`(\\alpha) \\mapsto` Controlled "
          ":math:`\\mathrm{Rx}(\\alpha)` gate")
      .value(
          "CRy", OpType::CRy,
          ":math:`(\\alpha) \\mapsto` Controlled "
          ":math:`\\mathrm{Ry}(\\alpha)` gate")
      .value(
          "CU1", OpType::CU1,
          ":math:`(\\lambda) \\mapsto` Controlled "
          ":math:`\\mathrm{U1}(\\lambda)` gate. Note that this is not "
          "equivalent to a :math:`\\mathrm{CRz}(\\lambda)` up to "
          "global phase, differing by an extra "
          ":math:`\\mathrm{Rz}(\\frac{\\lambda}{2})` on the control "
          "qubit.")
      .value(
          "CU3", OpType::CU3,
          ":math:`(\\theta, \\phi, \\lambda) \\mapsto` Controlled "
          ":math:`\\mathrm{U3}(\\theta, \\phi, \\lambda)` gate. "
          "Similar rules apply.")
      .value("CCX", OpType::CCX, "Toffoli gate")
      .value(
          "ECR", OpType::ECR,
          ":math:`\\frac{1}{\\sqrt 2} \\left[ "
          "\\begin{array}{cccc} 0 & 0 & 1 & i \\\\"
          "0 & 0 & i & 1 \\\\"
          "1 & -i & 0 & 0 \\\\"
          "-i & 1 & 0 & 0 \\end{array} \\right]`")
      .value("SWAP", OpType::SWAP, "Swap gate")
      .value("CSWAP", OpType::CSWAP, "Controlled swap gate")
      .value(
          "noop", OpType::noop,
          "Identity gate. These gates are not permanent and are "
          "automatically stripped by the compiler")
      .value(
          "Barrier", OpType::Barrier,
          "Meta-operation preventing compilation through it. Not "
          "automatically stripped by the compiler")
      .value(
          "Label", OpType::Label,
          "Label for control flow jumps. Does not appear within a "
          "circuit")
      .value(
          "Branch", OpType::Branch,
          "A control flow jump to a label dependent on the value of a "
          "given Bit. Does not appear within a circuit")
      .value(
          "Goto", OpType::Goto,
          "An unconditional control flow jump to a Label. Does not "
          "appear within a circuit.")
      .value(
          "Stop", OpType::Stop,
          "Halts execution immediately. Used to terminate a program. "
          "Does not appear within a circuit.")
      .value(
          "BRIDGE", OpType::BRIDGE,
          "A CX Bridge over 3 qubits. Used to apply a logical CX "
          "between the first and third qubits when they are not "
          "adjacent on the device, but both neighbour the second "
          "qubit. Acts as the identity on the second qubit")
      .value(
          "Measure", OpType::Measure,
          "Z-basis projective measurement, storing the measurement "
          "outcome in a specified bit")
      .value(
          "Reset", OpType::Reset,
          "Resets the qubit to :math:`\\left|0\\right>`")
      .value("CircBox", OpType::CircBox, "Represents an arbitrary subcircuit")
      .value(
          "PhasePolyBox", OpType::PhasePolyBox,
          "An operation representing arbitrary circuits made up of CX and Rz "
          "gates, represented as a phase polynomial together with a boolean "
          "matrix representing an additional linear transformation.")
      .value(
          "Unitary1qBox", OpType::Unitary1qBox,
          "Represents an arbitrary one-qubit unitary operation by its "
          "matrix")
      .value(
          "Unitary2qBox", OpType::Unitary2qBox,
          "Represents an arbitrary two-qubit unitary operation by its "
          "matrix")
      .value(
          "Unitary3qBox", OpType::Unitary3qBox,
          "Represents an arbitrary three-qubit unitary operation by its matrix")
      .value(
          "ExpBox", OpType::ExpBox,
          "A two-qubit operation corresponding to a unitary matrix "
          "defined as the exponential :math:`e^{itA}` of an arbitrary "
          "4x4 hermitian matrix :math:`A`.")
      .value(
          "PauliExpBox", OpType::PauliExpBox,
          "An operation defined as the exponential "
          ":math:`e^{-\\frac{i\\pi\\alpha}{2} P}` of a tensor "
          ":math:`P` of Pauli operations.")
      .value(
          "QControlBox", OpType::QControlBox,
          "An arbitrary n-controlled operation")
      .value(
          "CustomGate", OpType::CustomGate,
          ":math:`(\\alpha, \\beta, \\ldots) \\mapsto` A user-defined "
          "operation, based on a :py:class:`Circuit` :math:`C` with "
          "parameters :math:`\\alpha, \\beta, \\ldots` substituted in "
          "place of bound symbolic variables in :math:`C`, as defined "
          "by the :py:class:`CustomGateDef`.")
      .value(
          "Conditional", OpType::Conditional,
          "An operation to be applied conditionally on the value of "
          "some classical register")
      .value(
          "ISWAP", OpType::ISWAP,
          ":math:`(\\alpha) \\mapsto e^{\\frac14 i \\pi\\alpha "
          "(\\mathrm{X} \\otimes \\mathrm{X} + \\mathrm{Y} \\otimes "
          "\\mathrm{Y})} = \\left[ \\begin{array}{cccc} 1 & 0 & 0 & 0 "
          "\\\\ 0 & \\cos\\frac{\\pi\\alpha}{2} & "
          "i\\sin\\frac{\\pi\\alpha}{2} & 0 \\\\ 0 & "
          "i\\sin\\frac{\\pi\\alpha}{2} & \\cos\\frac{\\pi\\alpha}{2} "
          "& 0 \\\\ 0 & 0 & 0 & 1 \\end{array} \\right]`")
      .value(
          "PhasedISWAP", OpType::PhasedISWAP,
          ":math:`(p, t) \\mapsto \\left[ \\begin{array}{cccc} 1 & 0 "
          "& 0 & 0 \\\\ 0 & \\cos\\frac{\\pi t}{2} & "
          "i\\sin\\frac{\\pi t}{2}e^{2i\\pi p} & 0 \\\\ 0 & "
          "i\\sin\\frac{\\pi t}{2}e^{-2i\\pi p} & \\cos\\frac{\\pi "
          "t}{2} & 0 \\\\ 0 & 0 & 0 & 1 \\end{array} \\right]` "
          "(equivalent to: Rz(p)[0]; Rz(-p)[1]; ISWAP(t); "
          "Rz(-p)[0]; Rz(p)[1])")
      .value(
          "XXPhase", OpType::XXPhase,
          ":math:`(\\alpha) \\mapsto e^{-\\frac12 i \\pi\\alpha "
          "(\\mathrm{X} \\otimes \\mathrm{X})} = \\left[ "
          "\\begin{array}{cccc} \\cos\\frac{\\pi\\alpha}{2} & 0 & 0 & "
          "-i\\sin\\frac{\\pi\\alpha}{2} \\\\ 0 & "
          "\\cos\\frac{\\pi\\alpha}{2} & "
          "-i\\sin\\frac{\\pi\\alpha}{2} & 0 \\\\ 0 & "
          "-i\\sin\\frac{\\pi\\alpha}{2} & "
          "\\cos\\frac{\\pi\\alpha}{2} & 0 \\\\ "
          "-i\\sin\\frac{\\pi\\alpha}{2} & 0 & 0 & "
          "\\cos\\frac{\\pi\\alpha}{2} \\end{array} \\right]`")
      .value(
          "YYPhase", OpType::YYPhase,
          ":math:`(\\alpha) \\mapsto e^{-\\frac12 i \\pi\\alpha "
          "(\\mathrm{Y} \\otimes \\mathrm{Y})} = \\left[ "
          "\\begin{array}{cccc} \\cos\\frac{\\pi\\alpha}{2} & 0 & 0 & "
          "i\\sin\\frac{\\pi\\alpha}{2} \\\\ 0 & "
          "\\cos\\frac{\\pi\\alpha}{2} & "
          "-i\\sin\\frac{\\pi\\alpha}{2} & 0 \\\\ 0 & "
          "-i\\sin\\frac{\\pi\\alpha}{2} & "
          "\\cos\\frac{\\pi\\alpha}{2} & 0 \\\\ "
          "i\\sin\\frac{\\pi\\alpha}{2} & 0 & 0 & "
          "\\cos\\frac{\\pi\\alpha}{2} \\end{array} \\right]`")
      .value(
          "ZZPhase", OpType::ZZPhase,
          ":math:`(\\alpha) \\mapsto e^{-\\frac12 i \\pi\\alpha "
          "(\\mathrm{Z} \\otimes \\mathrm{Z})} = \\left[ "
          "\\begin{array}{cccc} e^{-\\frac12 i \\pi\\alpha} & 0 & 0 & "
          "0 \\\\ 0 & e^{\\frac12 i \\pi\\alpha} & 0 & 0 \\\\ 0 & 0 & "
          "e^{\\frac12 i \\pi\\alpha} & 0 \\\\ 0 & 0 & 0 & "
          "e^{-\\frac12 i \\pi\\alpha} \\end{array} \\right]`")
      .value(
          "XXPhase3", OpType::XXPhase3,
          "A 3-qubit gate XXPhase3(α) consists of pairwise 2-qubit XXPhase(α) "
          "interactions. "
          "Equivalent to XXPhase(α)[0, 1] XXPhase(α)[1, 2] XXPhase(α)[0, 2].")
      .value(
          "PhasedX", OpType::PhasedX,
          ":math:`(\\alpha,\\beta) \\mapsto "
          "\\mathrm{Rz}(\\beta)\\mathrm{Rx}(\\alpha)\\mathrm{Rz}(-"
          "\\beta)` (matrix-multiplication order)")
      .value(
          "NPhasedX", OpType::NPhasedX,
          ":math:`(\\alpha, \\beta) \\mapsto \\mathrm{PhasedX}(\\alpha, \\beta)"
          "^{\\otimes n}` (n-qubit gate composed of identical PhasedX in "
          "parallel.")
      .value(
          "CnRy", OpType::CnRy,
          ":math:`(\\alpha)` := n-controlled "
          ":math:`\\mathrm{Ry}(\\alpha)` gate.")
      .value("CnX", OpType::CnX, "n-controlled X gate.")
      .value(
          "ZZMax", OpType::ZZMax,
          ":math:`e^{-\\frac{i\\pi}{4}(\\mathrm{Z} \\otimes "
          "\\mathrm{Z})}`, a maximally entangling ZZPhase")
      .value(
          "ESWAP", OpType::ESWAP,
          ":math:`\\alpha \\mapsto e^{-\\frac12 i\\pi\\alpha \\cdot "
          "\\mathrm{SWAP}} = \\left[ \\begin{array}{cccc} "
          "e^{-\\frac12 i \\pi\\alpha} & 0 & 0 & 0 \\\\ 0 & "
          "\\cos\\frac{\\pi\\alpha}{2} & "
          "-i\\sin\\frac{\\pi\\alpha}{2} & 0 \\\\ 0 & "
          "-i\\sin\\frac{\\pi\\alpha}{2} & "
          "\\cos\\frac{\\pi\\alpha}{2} & 0 \\\\ 0 & 0 & 0 & "
          "e^{-\\frac12 i \\pi\\alpha} \\end{array} \\right]`")
      .value(
          "FSim", OpType::FSim,
          ":math:`(\\alpha, \\beta) \\mapsto \\left[ "
          "\\begin{array}{cccc} 1 & 0 & 0 & 0 \\\\ 0 & "
          "\\cos \\pi\\alpha & "
          "-i\\sin \\pi\\alpha & 0 \\\\ 0 & "
          "-i\\sin \\pi\\alpha & "
          "\\cos \\pi\\alpha & 0 \\\\ 0 & 0 & 0 & "
          "e^{-i\\pi\\beta} \\end{array} \\right]`")
      .value(
          "Sycamore", OpType::Sycamore,
          ":math:`\\mathrm{FSim}(\\frac12, \\frac16)`")
      .value(
          "ISWAPMax", OpType::ISWAPMax,
          ":math:`\\mathrm{ISWAP}(1) = \\left[ \\begin{array}{cccc} 1 "
          "& 0 & 0 & 0 \\\\ 0 & 0 & i & 0 \\\\ 0 & i & 0 & 0 \\\\ 0 & "
          "0 & 0 & 1 \\end{array} \\right]`")
      .value(
          "ClassicalTransform", OpType::ClassicalTransform,
          "A general classical operation where all inputs are also outputs")
      .value(
          "WASM", OpType::WASM, "Op containing a classical wasm function call")
      .value(
          "SetBits", OpType::SetBits,
          "An operation to set some bits to specified values")
      .value(
          "CopyBits", OpType::CopyBits, "An operation to copy some bit values")
      .value(
          "RangePredicate", OpType::RangePredicate,
          "A classical predicate defined by a range of values in binary "
          "encoding")
      .value(
          "ExplicitPredicate", OpType::ExplicitPredicate,
          "A classical predicate defined by a truth table")
      .value(
          "ExplicitModifier", OpType::ExplicitModifier,
          "An operation defined by a truth table that modifies one bit")
      .value(
          "MultiBit", OpType::MultiBit,
          "A classical operation applied to multiple bits simultaneously")
      .value(
          "ClassicalExpBox", OpType::ClassicalExpBox,
          "A box for holding compound classical operations on Bits.")
      .def_static(
          "from_name", [](const json &j) { return j.get<OpType>(); },
          "Construct from name");
  py::enum_<BasisOrder>(
      m, "BasisOrder",
      "Enum for readout basis and ordering.\n"
      "Readouts are viewed in increasing lexicographic order (ILO) of "
      "the bit's UnitID. This is our default convention for column "
      "indexing for ALL readout forms (shots, counts, statevector, and "
      "unitaries). e.g. :math:`\\lvert abc \\rangle` corresponds to the "
      "readout: ('c', 0) --> :math:`a`, ('c', 1) --> :math:`b`, ('d', 0) "
      "--> :math:`c`\n"
      "For statevector and unitaries, the string abc is interpreted as "
      "an index in a big-endian (BE) fashion. e.g. the statevector "
      ":math:`(a_{00}, a_{01}, a_{10}, a_{11})`\n"
      "Some backends (Qiskit, ProjectQ, etc.) use a DLO-BE (decreasing "
      "lexicographic order, big-endian) convention. This is the same as "
      "ILO-LE (little-endian) for statevectors and unitaries, but gives "
      "shot tables/readouts in a counter-intuitive manner.\n"
      "Every backend and matrix-based box has a BasisOrder option which "
      "can toggle between ILO-BE (ilo) and DLO-BE (dlo).")
      .value(
          "ilo", BasisOrder::ilo,
          "Increasing Lexicographic Order of UnitID, big-endian")
      .value(
          "dlo", BasisOrder::dlo,
          "Decreasing Lexicographic Order of UnitID, big-endian");

  py::class_<Command>(
      m, "Command",
      "A single quantum command in the circuit, defined by the Op, the "
      "qubits it acts on, and the op group name if any.")
      .def(
          py::init<const Op_ptr, unit_vector_t>(),
          "Construct from an operation and a vector of unit IDs", py::arg("op"),
          py::arg("args"))
      .def("__eq__", &Command::operator==)
      .def("__repr__", &Command::to_str)
      .def_property_readonly(
          "op", &Command::get_op_ptr, "Operation for this command.")
      .def_property_readonly(
          "args", &Command::get_args, "The qubits/bits the command acts on.")
      .def_property_readonly(
          "qubits", &Command::get_qubits, "The qubits the command acts on.")
      .def_property_readonly(
          "bits", &Command::get_bits,
          "The bits the command could write to (does not include "
          "read-only bits).")
      .def_property_readonly(
          "opgroup", &Command::get_opgroup,
          "The op group name assigned to the command (or `None` if "
          "no name is defined).")
      .def(
          "free_symbols",
          [](const Command &com) { return com.get_op_ptr()->free_symbols(); },
          ":return: set of symbolic parameters for the command");

  init_library(m);
  init_boxes(m);
  init_classical(m);
  init_circuit(m);

  m.def(
      "fresh_symbol", &SymTable::fresh_symbol,
      "Given some preferred symbol, this finds an appropriate suffix "
      "that "
      "will guarantee it has not yet been used in the current python "
      "session.\n\n:param preferred: The preferred readable symbol name "
      "as "
      "a string (default is 'a')\n\n:return: A new sympy symbol object",
      py::arg("preferred") = 'a');
}

}  // namespace tket
