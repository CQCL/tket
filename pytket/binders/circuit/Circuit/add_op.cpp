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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <optional>
#include <vector>

#include "UnitRegister.hpp"
#include "add_gate.hpp"
#include "tket/Circuit/Boxes.hpp"
#include "tket/Circuit/Circuit.hpp"
#include "tket/Circuit/ClassicalExpBox.hpp"
#include "tket/Circuit/ConjugationBox.hpp"
#include "tket/Circuit/DiagonalBox.hpp"
#include "tket/Circuit/Multiplexor.hpp"
#include "tket/Circuit/PauliExpBoxes.hpp"
#include "tket/Circuit/StatePreparation.hpp"
#include "tket/Circuit/ToffoliBox.hpp"
#include "tket/Converters/PhasePoly.hpp"
#include "tket/Gate/OpPtrFunctions.hpp"
#include "typecast.hpp"
namespace py = pybind11;

namespace tket {

const bit_vector_t no_bits;

template <typename ID>
static Circuit *add_gate_method_noparams(
    Circuit *circ, OpType type, const py::tket_custom::SequenceVec<ID> &args,
    const py::kwargs &kwargs) {
  return add_gate_method(
      circ, get_op_ptr(type, std::vector<Expr>{}, args.size()), args, kwargs);
}

template <typename ID>
static Circuit *add_gate_method_oneparam(
    Circuit *circ, OpType type, const Expr &p, const py::tket_custom::SequenceVec<ID> &args,
    const py::kwargs &kwargs) {
  return add_gate_method(circ, get_op_ptr(type, p, args.size()), args, kwargs);
}

template <typename ID>
static Circuit *add_gate_method_manyparams(
    Circuit *circ, OpType type, const py::tket_custom::SequenceVec<Expr> &ps,
    const py::tket_custom::SequenceVec<ID> &args, const py::kwargs &kwargs) {
  return add_gate_method(circ, get_op_ptr(type, ps, args.size()), args, kwargs);
}

template <typename ID>
static Circuit *add_box_method(
    Circuit *circ, Op_ptr box_ptr, const std::vector<ID> &args,
    const py::kwargs &kwargs) {
  return add_gate_method(circ, box_ptr, args, kwargs);
}

void init_circuit_add_op(py::class_<Circuit, std::shared_ptr<Circuit>> &c) {
  c.def(
       "add_gate", &add_gate_method<unsigned>,
       "Appends a single operation to the end of the circuit on some "
       "particular qubits/bits. The number of qubits/bits specified "
       "must match the arity of the gate.",
       py::arg("Op"), py::arg("args"))
      .def(
          "add_gate", &add_gate_method<UnitID>,
          "Appends a single operation to the end of the circuit on some "
          "particular qubits/bits. The number of qubits/bits specified "
          "must match the arity of the gate.",
          py::arg("Op"), py::arg("args"))
      .def(
          "add_gate", &add_gate_method_noparams<unsigned>,
          "Appends a single (non-parameterised) gate to the end of the "
          "circuit on some particular qubits from the default register "
          "('q'). The number of qubits specified must match the arity "
          "of the gate. For `OpType.Measure` operations the bit from "
          "the default register should follow the qubit."
          "\n\n>>> c.add_gate(OpType.H, [0]) # equivalent to "
          "c.H(0)\n>>> c.add_gate(OpType.CX, [0,1]) # equivalent to "
          "c.CX(0,1)"
          "\n\n:param type: The type of operation to add"
          "\n:param args: The list of indices for the qubits/bits to "
          "which the operation is applied"
          "\n:param kwargs: Additional properties for classical "
          "conditions"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("type"), py::arg("args"))
      .def(
          "add_gate", &add_gate_method_noparams<UnitID>,
          "Appends a single (non-parameterised) gate to the end of the "
          "circuit on some particular qubits from the default register "
          "('q'). The number of qubits specified must match the arity "
          "of the gate. For `OpType.Measure` operations the bit from "
          "the default register should follow the qubit."
          "\n\n>>> c.add_gate(OpType.H, [0]) # equivalent to "
          "c.H(0)\n>>> c.add_gate(OpType.CX, [0,1]) # equivalent to "
          "c.CX(0,1)"
          "\n\n:param type: The type of operation to add"
          "\n:param args: The qubits/bits to apply the gate to"
          "\n:param kwargs: Additional properties for classical "
          "conditions"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("type"), py::arg("args"))
      .def(
          "add_gate", &add_gate_method_oneparam<unsigned>,
          "Appends a single gate, parameterised by an expression, to "
          "the end of circuit on some particular qubits from the "
          "default register ('q')."
          "\n\n:param type: The type of gate to add"
          "\n:param angle: The parameter for the gate in halfturns"
          "\n:param args: The list of indices for the qubits to which "
          "the operation is applied"
          "\n:param kwargs: Additional properties for classical "
          "conditions"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("type"), py::arg("angle"), py::arg("args"))
      .def(
          "add_gate", &add_gate_method_oneparam<UnitID>,
          "Appends a single gate, parameterised by an expression, to "
          "the end of circuit on some particular qubits from the "
          "default register ('q')."
          "\n\n:param type: The type of gate to add"
          "\n:param angle: The parameter for the gate in halfturns"
          "\n:param args: The qubits/bits to apply the gate to"
          "\n:param kwargs: Additional properties for classical "
          "conditions"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("type"), py::arg("angle"), py::arg("args"))
      .def(
          "add_gate", &add_gate_method_manyparams<unsigned>,
          "Appends a single gate, parameterised with a vector of "
          "expressions corresponding to halfturns, to the end of "
          "circuit on some particular qubits from the default register "
          "('q')."
          "\n\n:param type: The type of gate to add"
          "\n:param angles: The parameters for the gate in halfturns"
          "\n:param args: The list of indices for the qubits to which "
          "the operation is applied"
          "\n:param kwargs: Additional properties for classical "
          "conditions"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("type"), py::arg("angles"), py::arg("args"))
      .def(
          "add_gate", &add_gate_method_manyparams<UnitID>,
          "Appends a single gate to the end of the circuit"
          "\n\n:param type: The type of gate to add"
          "\n:param params: The parameters for the gate in halfturns"
          "\n:param args: The qubits/bits to apply the gate to"
          "\n:param kwargs: Additional properties for classical "
          "conditions"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("type"), py::arg("angles"), py::arg("args"))
      .def(
          "add_barrier",
          [](Circuit *circ, const py::tket_custom::SequenceVec<unsigned> &qubits,
             const py::tket_custom::SequenceVec<unsigned> &bits, const std::string &data) {
            circ->add_barrier(qubits, bits, data);
            return circ;
          },
          "Append a Barrier on the given units"
          "\n\n:param data: additional data stored in the barrier"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("qubits"), py::arg("bits") = no_bits, py::arg("data") = "")
      .def(
          "add_conditional_barrier",
          [](Circuit *circ, const py::tket_custom::SequenceVec<unsigned> &barrier_qubits,
             const py::tket_custom::SequenceVec<unsigned> &barrier_bits,
             const py::tket_custom::SequenceVec<unsigned> &condition_bits, unsigned value,
             const std::string &_data) {
            circ->add_conditional_barrier(
                barrier_qubits, barrier_bits, condition_bits, value, _data);
            return circ;
          },
          "Append a Conditional Barrier on the given barrier qubits and "
          "barrier bits, conditioned on the given condition bits."
          "\n\n:param barrier_qubits: Qubit in Barrier operation."
          "\n:param barrier_bits: Bit in Barrier operation."
          "\n:param condition_bits: Bit covering classical control condition "
          "of barrier operation."
          "\n:param value: Value that classical condition must have to "
          "hold (little-endian)."
          "\n:param data: Additional data stored in Barrier operation."
          "\n:return: the new :py:class:`Circuit`",
          py::arg("barrier_qubits"), py::arg("barrier_bits"),
          py::arg("condition_bits"), py::arg("value"), py::arg("data") = "")
      .def(
          "add_circbox",
          [](Circuit *circ, const CircBox &box,
             const py::tket_custom::SequenceVec<unsigned> &args, const py::kwargs &kwargs) {
            return add_box_method(
                circ, std::make_shared<CircBox>(box), args, kwargs);
          },
          "Append a :py:class:`CircBox` to the circuit.\n\n:param "
          "circbox: The box to append\n:param args: Indices of the "
          "qubits/bits to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("circbox"), py::arg("args"))
      .def(
          "add_unitary1qbox",
          [](Circuit *circ, const Unitary1qBox &box, unsigned q0,
             const py::kwargs &kwargs) {
            return add_box_method<unsigned>(
                circ, std::make_shared<Unitary1qBox>(box), {q0}, kwargs);
          },
          "Append a :py:class:`Unitary1qBox` to the "
          "circuit.\n\n:param "
          "unitarybox: The box to append\n:param qubit_0: Index of "
          "the qubit to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("unitarybox"), py::arg("qubit_0"))
      .def(
          "add_unitary2qbox",
          [](Circuit *circ, const Unitary2qBox &box, unsigned q0, unsigned q1,
             const py::kwargs &kwargs) {
            return add_box_method<unsigned>(
                circ, std::make_shared<Unitary2qBox>(box), {q0, q1}, kwargs);
          },
          "Append a :py:class:`Unitary2qBox` to the circuit.\n\nThe "
          "matrix representation is ILO-BE.\n\n:param unitarybox: "
          "The box to append\n:param qubit_0: Index of the first "
          "target "
          "qubit\n:param qubit_1: Index of the second target qubit"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("unitarybox"), py::arg("qubit_0"), py::arg("qubit_1"))
      .def(
          "add_unitary3qbox",
          [](Circuit *circ, const Unitary3qBox &box, unsigned q0, unsigned q1,
             unsigned q2, const py::kwargs &kwargs) {
            return add_box_method<unsigned>(
                circ, std::make_shared<Unitary3qBox>(box), {q0, q1, q2},
                kwargs);
          },
          "Append a :py:class:`Unitary3qBox` to the circuit."
          "\n\n:param unitarybox: box to append"
          "\n:param qubit_0: index of target qubit 0"
          "\n:param qubit_1: index of target qubit 1"
          "\n:param qubit_2: index of target qubit 2"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("unitarybox"), py::arg("qubit_0"), py::arg("qubit_1"),
          py::arg("qubit_2"))
      .def(
          "add_expbox",
          [](Circuit *circ, const ExpBox &box, unsigned q0, unsigned q1,
             const py::kwargs &kwargs) {
            return add_box_method<unsigned>(
                circ, std::make_shared<ExpBox>(box), {q0, q1}, kwargs);
          },
          "Append an :py:class:`ExpBox` to the circuit.\n\nThe "
          "matrix representation is ILO-BE.\n\n:param expbox: The "
          "box to append\n:param qubit_0: Index of the first target "
          "qubit\n:param qubit_1: Index of the second target qubit"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("expbox"), py::arg("qubit_0"), py::arg("qubit_1"))
      .def(
          "add_pauliexpbox",
          [](Circuit *circ, const PauliExpBox &box,
             const py::tket_custom::SequenceVec<unsigned> &qubits, const py::kwargs &kwargs) {
            return add_box_method<unsigned>(
                circ, std::make_shared<PauliExpBox>(box), qubits, kwargs);
          },
          "Append a :py:class:`PauliExpBox` to the "
          "circuit.\n\n:param pauliexpbox: The box to append\n:param "
          "qubits: Indices of the qubits to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("pauliexpbox"), py::arg("qubits"))
      .def(
          "add_pauliexppairbox",
          [](Circuit *circ, const PauliExpPairBox &box,
             const py::tket_custom::SequenceVec<unsigned> &qubits, const py::kwargs &kwargs) {
            return add_box_method<unsigned>(
                circ, std::make_shared<PauliExpPairBox>(box), qubits, kwargs);
          },
          "Append a :py:class:`PauliExpPairBox` to the "
          "circuit.\n\n:param pauliexppairbox: The box to append\n:param "
          "qubits: Indices of the qubits to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("pauliexppairbox"), py::arg("qubits"))
      .def(
          "add_pauliexpcommutingsetbox",
          [](Circuit *circ, const PauliExpCommutingSetBox &box,
             const py::tket_custom::SequenceVec<unsigned> &qubits, const py::kwargs &kwargs) {
            return add_box_method<unsigned>(
                circ, std::make_shared<PauliExpCommutingSetBox>(box), qubits,
                kwargs);
          },
          "Append a :py:class:`PauliExpCommutingSetBox` to the "
          "circuit.\n\n:param pauliexpcommutingsetbox: The box to "
          "append\n:param "
          "qubits: Indices of the qubits to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("pauliexpcommutingsetbox"), py::arg("qubits"))
      .def(
          "add_toffolibox",
          [](Circuit *circ, const ToffoliBox &box,
             const py::tket_custom::SequenceVec<unsigned> &qubits, const py::kwargs &kwargs) {
            return add_box_method<unsigned>(
                circ, std::make_shared<ToffoliBox>(box), qubits, kwargs);
          },
          "Append a :py:class:`ToffoliBox` to the "
          "circuit.\n\n:param toffolibox: The box to append\n:param "
          "qubits: Indices of the qubits to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("toffolibox"), py::arg("qubits"))
      .def(
          "add_qcontrolbox",
          [](Circuit *circ, const QControlBox &box,
             const py::tket_custom::SequenceVec<unsigned> &args, const py::kwargs &kwargs) {
            return add_box_method<unsigned>(
                circ, std::make_shared<QControlBox>(box), args, kwargs);
          },
          "Append a :py:class:`QControlBox` to the circuit.\n\n"
          ":param qcontrolbox: The box to append\n"
          ":param args: Indices of the qubits to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("qcontrolbox"), py::arg("args"))
      .def(
          "add_phasepolybox",
          [](Circuit *circ, const PhasePolyBox &box,
             const py::tket_custom::SequenceVec<unsigned> &qubits, const py::kwargs &kwargs) {
            return add_box_method<unsigned>(
                circ, std::make_shared<PhasePolyBox>(box), qubits, kwargs);
          },
          "Append a :py:class:`PhasePolyBox` to the "
          "circuit.\n\n:param phasepolybox: The box to append\n:param "
          "qubits: Indices of the qubits to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("phasepolybox"), py::arg("qubits"))

      .def(
          "add_classicalexpbox_bit",
          [](Circuit *circ, const py::object &exp,
             const py::tket_custom::SequenceVec<Bit> &outputs, const py::kwargs &kwargs) {
            auto inputs = exp.attr("all_inputs")().cast<std::set<Bit>>();
              py::tket_custom::SequenceVec<Bit> o_vec, io_vec;

            // if outputs are also in inputs, add to i/o wires
            for (const Bit &out : outputs) {
              auto find = inputs.find(out);
              if (find != inputs.end()) {
                inputs.erase(find);
                io_vec.push_back(out);
              } else {
                o_vec.push_back(out);
              }
            }
            unsigned n_i = inputs.size();
            unsigned n_io = io_vec.size();
            unsigned n_o = o_vec.size();
            o_vec.insert(o_vec.begin(), io_vec.begin(), io_vec.end());
            o_vec.insert(o_vec.begin(), inputs.begin(), inputs.end());
            return add_box_method<Bit>(
                circ,
                std::make_shared<ClassicalExpBox<py::object>>(
                    n_i, n_io, n_o, exp),
                o_vec, kwargs);
          },
          "Append a :py:class:`ClassicalExpBox` over Bit to the circuit.\n\n"
          ":param classicalexpbox: The box to append\n"
          ":param args: Indices of the qubits to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("expression"), py::arg("target"))
      .def(
          "add_classicalexpbox_register",
          [](Circuit *circ, const py::object &exp,
             const py::tket_custom::SequenceVec<Bit> &outputs, const py::kwargs &kwargs) {
            auto inputs =
                exp.attr("all_inputs")().cast<std::set<BitRegister>>();
            std::set<Bit> all_bits;
            for (const BitRegister &reg : inputs) {
              for (std::size_t i = 0; i < reg.size(); i++) {
                all_bits.insert(reg[i]);
              }
            }
            py::tket_custom::SequenceVec<Bit> o_vec, io_vec;
            for (const Bit &out : outputs) {
              auto find = all_bits.find(out);
              if (find != all_bits.end()) {
                all_bits.erase(find);
                io_vec.push_back(out);
              } else {
                o_vec.push_back(out);
              }
            }
            unsigned n_i = all_bits.size();
            unsigned n_io = io_vec.size();
            unsigned n_o = o_vec.size();
            o_vec.insert(o_vec.begin(), io_vec.begin(), io_vec.end());
            o_vec.insert(o_vec.begin(), all_bits.begin(), all_bits.end());
            return add_box_method<Bit>(
                circ,
                std::make_shared<ClassicalExpBox<py::object>>(
                    n_i, n_io, n_o, exp),
                o_vec, kwargs);
          },
          "Append a :py:class:`ClassicalExpBox` over BitRegister to the "
          "circuit.\n\n"
          ":param classicalexpbox: The box to append\n"
          ":param args: Indices of the qubits to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("expression"), py::arg("target"))
      .def(
          "add_custom_gate",
          [](Circuit *circ, const composite_def_ptr_t &definition,
             const py::tket_custom::SequenceVec<Expr> &params,
             const py::tket_custom::SequenceVec<unsigned> &qubits, const py::kwargs &kwargs) {
            return add_box_method<unsigned>(
                circ, std::make_shared<CustomGate>(definition, params), qubits,
                kwargs);
          },
          "Append an instance of a :py:class:`CustomGateDef` to the "
          "circuit.\n\n:param def: The custom gate "
          "definition\n:param params: List of parameters to "
          "instantiate the gate with, in halfturns\n:param qubits: "
          "Indices of the qubits to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("definition"), py::arg("params"), py::arg("qubits"))
      .def(
          "add_barrier",
          [](Circuit *circ, const unit_vector_t &units,
             const std::string &data) {
            circ->add_barrier(units, data);
            return circ;
          },
          "Append a Barrier on the given units"
          "\n\n:param data: additional data stored in the barrier"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("units"), py::arg("data") = "")
      .def(
          "add_conditional_barrier",
          [](Circuit *circ, const unit_vector_t &barrier_args,
             const bit_vector_t &condition_bits, unsigned value,
             const std::string &_data) {
            circ->add_conditional_barrier(
                barrier_args, condition_bits, value, _data);
            return circ;
          },
          "Append a Conditional Barrier on the given barrier qubits and "
          "barrier bits, conditioned on the given condition bits."
          "\n\n:param barrier_args: Qubit and Bit in Barrier operation."
          "\n:param condition_bits: Bit covering classical control "
          " condition of barrier operation."
          "\n:param value: Value that classical condition must have to "
          "hold (little-endian)."
          "\n:param data: Additional data stored in Barrier operation."
          "\n:return: the new :py:class:`Circuit`",
          py::arg("barrier_args"), py::arg("condition_bits"), py::arg("value"),
          py::arg("data") = "")
      .def(
          "add_circbox",
          [](Circuit *circ, const CircBox &box, const unit_vector_t &args,
             const py::kwargs &kwargs) {
            return add_box_method(
                circ, std::make_shared<CircBox>(box), args, kwargs);
          },
          "Append a :py:class:`CircBox` to the circuit.\n\n:param "
          "circbox: The box to append\n:param args: The qubits/bits "
          "to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("circbox"), py::arg("args"))
      .def(
          "add_unitary1qbox",
          [](Circuit *circ, const Unitary1qBox &box, const Qubit &q0,
             const py::kwargs &kwargs) {
            return add_box_method<UnitID>(
                circ, std::make_shared<Unitary1qBox>(box), {q0}, kwargs);
          },
          "Append a :py:class:`Unitary1qBox` to the "
          "circuit.\n\n:param unitarybox: The box to append\n:param "
          "qubit_0: The qubit to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("unitarybox"), py::arg("qubit_0"))
      .def(
          "add_unitary2qbox",
          [](Circuit *circ, const Unitary2qBox &box, const Qubit &q0,
             const Qubit &q1, const py::kwargs &kwargs) {
            return add_box_method<UnitID>(
                circ, std::make_shared<Unitary2qBox>(box), {q0, q1}, kwargs);
          },
          "Append a :py:class:`Unitary2qBox` to the circuit.\n\nThe "
          "matrix representation is ILO-BE.\n\n:param unitarybox: "
          "The box to append\n:param qubit_0: The first target "
          "qubit\n:param qubit_1: The second target qubit"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("unitarybox"), py::arg("qubit_0"), py::arg("qubit_1"))
      .def(
          "add_unitary3qbox",
          [](Circuit *circ, const Unitary3qBox &box, const Qubit &q0,
             const Qubit &q1, const Qubit &q2, const py::kwargs &kwargs) {
            return add_box_method<UnitID>(
                circ, std::make_shared<Unitary3qBox>(box), {q0, q1, q2},
                kwargs);
          },
          "Append a :py:class:`Unitary3qBox` to the circuit."
          "\n\n:param unitarybox: box to append"
          "\n:param qubit_0: index of target qubit 0"
          "\n:param qubit_1: index of target qubit 1"
          "\n:param qubit_2: index of target qubit 2"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("unitarybox"), py::arg("qubit_0"), py::arg("qubit_1"),
          py::arg("qubit_2"))
      .def(
          "add_expbox",
          [](Circuit *circ, const ExpBox &box, const Qubit &q0, const Qubit &q1,
             const py::kwargs &kwargs) {
            return add_box_method<UnitID>(
                circ, std::make_shared<ExpBox>(box), {q0, q1}, kwargs);
          },
          "Append an :py:class:`ExpBox` to the circuit.\n\nThe "
          "matrix representation is ILO-BE.\n\n:param expbox: The "
          "box to append\n:param qubit_0: The first target "
          "qubit\n:param qubit_1: The second target qubit"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("expbox"), py::arg("qubit_0"), py::arg("qubit_1"))
      .def(
          "add_pauliexpbox",
          [](Circuit *circ, const PauliExpBox &box,
             const qubit_vector_t &qubits, const py::kwargs &kwargs) {
            return add_box_method<UnitID>(
                circ, std::make_shared<PauliExpBox>(box),
                {qubits.begin(), qubits.end()}, kwargs);
          },
          "Append a :py:class:`PauliExpBox` to the "
          "circuit.\n\n:param pauliexpbox: The box to append\n:param "
          "qubits: The qubits to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("pauliexpbox"), py::arg("qubits"))
      .def(
          "add_pauliexppairbox",
          [](Circuit *circ, const PauliExpPairBox &box,
             const qubit_vector_t &qubits, const py::kwargs &kwargs) {
            return add_box_method<UnitID>(
                circ, std::make_shared<PauliExpPairBox>(box),
                {qubits.begin(), qubits.end()}, kwargs);
          },
          "Append a :py:class:`PauliExpPairBox` to the "
          "circuit.\n\n:param pauliexppairbox: The box to append\n:param "
          "qubits: The qubits to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("pauliexppairbox"), py::arg("qubits"))
      .def(
          "add_pauliexpcommutingsetbox",
          [](Circuit *circ, const PauliExpCommutingSetBox &box,
             const qubit_vector_t &qubits, const py::kwargs &kwargs) {
            return add_box_method<UnitID>(
                circ, std::make_shared<PauliExpCommutingSetBox>(box),
                {qubits.begin(), qubits.end()}, kwargs);
          },
          "Append a :py:class:`PauliExpCommutingSetBox` to the "
          "circuit.\n\n:param pauliexpcommutingsetbox: The box to "
          "append\n:param "
          "qubits: The qubits to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("pauliexpcommutingsetbox"), py::arg("qubits"))
      .def(
          "add_toffolibox",
          [](Circuit *circ, const ToffoliBox &box, const qubit_vector_t &qubits,
             const py::kwargs &kwargs) {
            return add_box_method<UnitID>(
                circ, std::make_shared<ToffoliBox>(box),
                {qubits.begin(), qubits.end()}, kwargs);
          },
          "Append a :py:class:`ToffoliBox` to the "
          "circuit.\n\n:param toffolibox: The box to append\n:param "
          "qubits: Indices of the qubits to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("toffolibox"), py::arg("qubits"))
      .def(
          "add_qcontrolbox",
          [](Circuit *circ, const QControlBox &box, const unit_vector_t &args,
             const py::kwargs &kwargs) {
            return add_box_method(
                circ, std::make_shared<QControlBox>(box), args, kwargs);
          },
          "Append a :py:class:`QControlBox` to the circuit.\n\n"
          ":param qcontrolbox: The box to append\n"
          ":param args: The qubits to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("qcontrolbox"), py::arg("args"))
      .def(
          "add_phasepolybox",
          [](Circuit *circ, const PhasePolyBox &box,
             const qubit_vector_t &qubits, const py::kwargs &kwargs) {
            return add_box_method<UnitID>(
                circ, std::make_shared<PhasePolyBox>(box),
                {qubits.begin(), qubits.end()}, kwargs);
          },
          "Append a :py:class:`PhasePolyBox` to the "
          "circuit.\n\n:param phasepolybox: The box to append\n:param "
          "qubits: The qubits to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("phasepolybox"), py::arg("qubits"))
      .def(
          "add_custom_gate",
          [](Circuit *circ, const composite_def_ptr_t &definition,
             const py::tket_custom::SequenceVec<Expr> &params, const qubit_vector_t &qubits,
             const py::kwargs &kwargs) {
            return add_box_method<UnitID>(
                circ, std::make_shared<CustomGate>(definition, params),
                {qubits.begin(), qubits.end()}, kwargs);
          },
          "Append an instance of a :py:class:`CustomGateDef` to the "
          "circuit.\n\n:param def: The custom gate "
          "definition\n:param params: List of parameters to "
          "instantiate the gate with, in halfturns\n:param qubits: "
          "The qubits to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("definition"), py::arg("params"), py::arg("qubits"))
      .def(
          "add_assertion",
          [](Circuit *circ, const ProjectorAssertionBox &box,
             const py::tket_custom::SequenceVec<unsigned> &qubits,
             const std::optional<unsigned> &ancilla,
             const std::optional<std::string> &name) -> Circuit * {
            std::vector<Qubit> qubits_;
            qubits_.reserve(qubits.size());
            for (unsigned int qubit : qubits) {
              qubits_.emplace_back(qubit);
            }
            std::optional<Qubit> ancilla_;
            if (ancilla == std::nullopt) {
              ancilla_ = std::nullopt;
            } else {
              ancilla_ = Qubit(ancilla.value());
            }
            circ->add_assertion(box, qubits_, ancilla_, name);
            return circ;
          },
          "Append a :py:class:`ProjectorAssertionBox` to the circuit."
          "\n\n:param box: ProjectorAssertionBox to append"
          "\n:param qubits: indices of target qubits"
          "\n:param ancilla: index of ancilla qubit"
          "\n:param name: name used to identify this assertion"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("box"), py::arg("qubits"), py::arg("ancilla") = std::nullopt,
          py::arg("name") = std::nullopt)
      .def(
          "add_assertion",
          [](Circuit *circ, const ProjectorAssertionBox &box,
             const py::tket_custom::SequenceVec<Qubit> &qubits,
             const std::optional<Qubit> &ancilla,
             const std::optional<std::string> &name) -> Circuit * {
            circ->add_assertion(box, qubits, ancilla, name);
            return circ;
          },
          "Append a :py:class:`ProjectorAssertionBox` to the circuit."
          "\n\n:param box: ProjectorAssertionBox to append"
          "\n:param qubits: target qubits"
          "\n:param ancilla: ancilla qubit"
          "\n:param name: name used to identify this assertion"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("box"), py::arg("qubits"), py::arg("ancilla") = std::nullopt,
          py::arg("name") = std::nullopt)
      .def(
          "add_assertion",
          [](Circuit *circ, const StabiliserAssertionBox &box,
             const py::tket_custom::SequenceVec<unsigned> &qubits, const unsigned &ancilla,
             const std::optional<std::string> &name) -> Circuit * {
            std::vector<Qubit> qubits_;
            qubits_.reserve(qubits.size());
            for (unsigned int qubit : qubits) {
              qubits_.emplace_back(qubit);
            }
            Qubit ancilla_(ancilla);
            circ->add_assertion(box, qubits_, ancilla_, name);
            return circ;
          },
          "Append a :py:class:`StabiliserAssertionBox` to the circuit."
          "\n\n:param box: StabiliserAssertionBox to append"
          "\n:param qubits: indices of target qubits"
          "\n:param ancilla: index of ancilla qubit"
          "\n:param name: name used to identify this assertion"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("box"), py::arg("qubits"), py::arg("ancilla"),
          py::arg("name") = std::nullopt)
      .def(
          "add_assertion",
          [](Circuit *circ, const StabiliserAssertionBox &box,
             const py::tket_custom::SequenceVec<Qubit> &qubits, const Qubit &ancilla,
             const std::optional<std::string> &name) -> Circuit * {
            circ->add_assertion(box, qubits, ancilla, name);
            return circ;
          },
          "Append a :py:class:`StabiliserAssertionBox` to the circuit."
          "\n\n:param box: StabiliserAssertionBox to append"
          "\n:param qubits: target qubits"
          "\n:param ancilla: ancilla qubit"
          "\n:param name: name used to identify this assertion"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("box"), py::arg("qubits"), py::arg("ancilla"),
          py::arg("name") = std::nullopt)
      .def(
          "add_multiplexor",
          [](Circuit *circ, const MultiplexorBox &box,
             const unit_vector_t &args, const py::kwargs &kwargs) {
            return add_box_method(
                circ, std::make_shared<MultiplexorBox>(box), args, kwargs);
          },
          "Append a :py:class:`MultiplexorBox` to the circuit.\n\n"
          ":param box: The box to append\n"
          ":param args: The qubits to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("box"), py::arg("args"))
      .def(
          "add_multiplexor",
          [](Circuit *circ, const MultiplexorBox &box,
             const py::tket_custom::SequenceVec<unsigned> &args, const py::kwargs &kwargs) {
            return add_box_method(
                circ, std::make_shared<MultiplexorBox>(box), args, kwargs);
          },
          "Append a :py:class:`MultiplexorBox` to the circuit.\n\n"
          ":param box: The box to append\n"
          ":param args: Indices of the qubits to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("box"), py::arg("args"))
      .def(
          "add_multiplexedrotation",
          [](Circuit *circ, const MultiplexedRotationBox &box,
             const unit_vector_t &args, const py::kwargs &kwargs) {
            return add_box_method(
                circ, std::make_shared<MultiplexedRotationBox>(box), args,
                kwargs);
          },
          "Append a :py:class:`MultiplexedRotationBox` to the circuit.\n\n"
          ":param box: The box to append\n"
          ":param args: The qubits to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("box"), py::arg("args"))
      .def(
          "add_multiplexedrotation",
          [](Circuit *circ, const MultiplexedRotationBox &box,
             const py::tket_custom::SequenceVec<unsigned> &args, const py::kwargs &kwargs) {
            return add_box_method(
                circ, std::make_shared<MultiplexedRotationBox>(box), args,
                kwargs);
          },
          "Append a :py:class:`MultiplexedRotationBox` to the circuit.\n\n"
          ":param box: The box to append\n"
          ":param args: Indices of the qubits to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("box"), py::arg("args"))
      .def(
          "add_multiplexedu2",
          [](Circuit *circ, const MultiplexedU2Box &box,
             const unit_vector_t &args, const py::kwargs &kwargs) {
            return add_box_method(
                circ, std::make_shared<MultiplexedU2Box>(box), args, kwargs);
          },
          "Append a :py:class:`MultiplexedU2Box` to the circuit.\n\n"
          ":param box: The box to append\n"
          ":param args: The qubits to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("box"), py::arg("args"))
      .def(
          "add_multiplexedu2",
          [](Circuit *circ, const MultiplexedU2Box &box,
             const py::tket_custom::SequenceVec<unsigned> &args, const py::kwargs &kwargs) {
            return add_box_method(
                circ, std::make_shared<MultiplexedU2Box>(box), args, kwargs);
          },
          "Append a :py:class:`MultiplexedU2Box` to the circuit.\n\n"
          ":param box: The box to append\n"
          ":param args: Indices of the qubits to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("box"), py::arg("args"))
      .def(
          "add_multiplexed_tensored_u2",
          [](Circuit *circ, const MultiplexedTensoredU2Box &box,
             const unit_vector_t &args, const py::kwargs &kwargs) {
            return add_box_method(
                circ, std::make_shared<MultiplexedTensoredU2Box>(box), args,
                kwargs);
          },
          "Append a :py:class:`MultiplexedTensoredU2Box` to the circuit.\n\n"
          ":param box: The box to append\n"
          ":param args: The qubits to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("box"), py::arg("args"))
      .def(
          "add_multiplexed_tensored_u2",
          [](Circuit *circ, const MultiplexedTensoredU2Box &box,
             const py::tket_custom::SequenceVec<unsigned> &args, const py::kwargs &kwargs) {
            return add_box_method(
                circ, std::make_shared<MultiplexedTensoredU2Box>(box), args,
                kwargs);
          },
          "Append a :py:class:`MultiplexedTensoredU2Box` to the circuit.\n\n"
          ":param box: The box to append\n"
          ":param args: Indices of the qubits to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("box"), py::arg("args"))
      .def(
          "add_state_preparation_box",
          [](Circuit *circ, const StatePreparationBox &box,
             const unit_vector_t &args, const py::kwargs &kwargs) {
            return add_box_method(
                circ, std::make_shared<StatePreparationBox>(box), args, kwargs);
          },
          "Append a :py:class:`StatePreparationBox` to the circuit.\n\n"
          ":param box: The box to append\n"
          ":param args: The qubits to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("box"), py::arg("args"))
      .def(
          "add_state_preparation_box",
          [](Circuit *circ, const StatePreparationBox &box,
             const py::tket_custom::SequenceVec<unsigned> &args, const py::kwargs &kwargs) {
            return add_box_method(
                circ, std::make_shared<StatePreparationBox>(box), args, kwargs);
          },
          "Append a :py:class:`StatePreparationBox` to the circuit.\n\n"
          ":param box: The box to append\n"
          ":param args: Indices of the qubits to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("box"), py::arg("args"))
      .def(
          "add_diagonal_box",
          [](Circuit *circ, const DiagonalBox &box, const unit_vector_t &args,
             const py::kwargs &kwargs) {
            return add_box_method(
                circ, std::make_shared<DiagonalBox>(box), args, kwargs);
          },
          "Append a :py:class:`DiagonalBox` to the circuit.\n\n"
          ":param box: The box to append\n"
          ":param args: The qubits to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("box"), py::arg("args"))
      .def(
          "add_diagonal_box",
          [](Circuit *circ, const DiagonalBox &box,
             const py::tket_custom::SequenceVec<unsigned> &args, const py::kwargs &kwargs) {
            return add_box_method(
                circ, std::make_shared<DiagonalBox>(box), args, kwargs);
          },
          "Append a :py:class:`DiagonalBox` to the circuit.\n\n"
          ":param box: The box to append\n"
          ":param args: Indices of the qubits to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("box"), py::arg("args"))
      .def(
          "add_conjugation_box",
          [](Circuit *circ, const ConjugationBox &box,
             const unit_vector_t &args, const py::kwargs &kwargs) {
            return add_box_method(
                circ, std::make_shared<ConjugationBox>(box), args, kwargs);
          },
          "Append a :py:class:`ConjugationBox` to the circuit.\n\n"
          ":param box: The box to append\n"
          ":param args: The qubits to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("box"), py::arg("args"))
      .def(
          "add_conjugation_box",
          [](Circuit *circ, const ConjugationBox &box,
             const py::tket_custom::SequenceVec<unsigned> &args, const py::kwargs &kwargs) {
            return add_box_method(
                circ, std::make_shared<ConjugationBox>(box), args, kwargs);
          },
          "Append a :py:class:`ConjugationBox` to the circuit.\n\n"
          ":param box: The box to append\n"
          ":param args: Indices of the qubits to append the box to"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("box"), py::arg("args"))
      .def(
          "H",
          [](Circuit *circ, unsigned qb, const py::kwargs &kwargs) {
            return add_gate_method_noparams<unsigned>(
                circ, OpType::H, {qb}, kwargs);
          },
          "Appends a Hadamard gate."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit"))
      .def(
          "X",
          [](Circuit *circ, unsigned qb, const py::kwargs &kwargs) {
            return add_gate_method_noparams<unsigned>(
                circ, OpType::X, {qb}, kwargs);
          },
          "Appends an X gate.", "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit"))
      .def(
          "Y",
          [](Circuit *circ, unsigned qb, const py::kwargs &kwargs) {
            return add_gate_method_noparams<unsigned>(
                circ, OpType::Y, {qb}, kwargs);
          },
          "Appends a Y gate.", "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit"))
      .def(
          "Z",
          [](Circuit *circ, unsigned qb, const py::kwargs &kwargs) {
            return add_gate_method_noparams<unsigned>(
                circ, OpType::Z, {qb}, kwargs);
          },
          "Appends a Z gate.", "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit"))
      .def(
          "T",
          [](Circuit *circ, unsigned qb, const py::kwargs &kwargs) {
            return add_gate_method_noparams<unsigned>(
                circ, OpType::T, {qb}, kwargs);
          },
          "Appends a T gate (equivalent to U1(0.25,-))."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit"))
      .def(
          "Tdg",
          [](Circuit *circ, unsigned qb, const py::kwargs &kwargs) {
            return add_gate_method_noparams<unsigned>(
                circ, OpType::Tdg, {qb}, kwargs);
          },
          "Appends a T-dagger gate (equivalent to U1(-0.25,-))."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit"))
      .def(
          "S",
          [](Circuit *circ, unsigned qb, const py::kwargs &kwargs) {
            return add_gate_method_noparams<unsigned>(
                circ, OpType::S, {qb}, kwargs);
          },
          "Appends an S gate (equivalent to U1(0.5,-))."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit"))
      .def(
          "Sdg",
          [](Circuit *circ, unsigned qb, const py::kwargs &kwargs) {
            return add_gate_method_noparams<unsigned>(
                circ, OpType::Sdg, {qb}, kwargs);
          },
          "Appends an S-dagger gate (equivalent to U1(-0.5,-))."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit"))
      .def(
          "V",
          [](Circuit *circ, unsigned qb, const py::kwargs &kwargs) {
            return add_gate_method_noparams<unsigned>(
                circ, OpType::V, {qb}, kwargs);
          },
          "Appends a V gate (equivalent to Rx(0.5,-))."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit"))
      .def(
          "Vdg",
          [](Circuit *circ, unsigned qb, const py::kwargs &kwargs) {
            return add_gate_method_noparams<unsigned>(
                circ, OpType::Vdg, {qb}, kwargs);
          },
          "Appends a V-dagger gate (equivalent to Rx(-0.5,-))."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit"))
      .def(
          "SX",
          [](Circuit *circ, unsigned qb, const py::kwargs &kwargs) {
            return add_gate_method_noparams<unsigned>(
                circ, OpType::SX, {qb}, kwargs);
          },
          "Appends a SX gate (equivalent to Rx(0.5,-)"
          " up to a 0.25 global phase)."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit"))
      .def(
          "SXdg",
          [](Circuit *circ, unsigned qb, const py::kwargs &kwargs) {
            return add_gate_method_noparams<unsigned>(
                circ, OpType::SXdg, {qb}, kwargs);
          },
          "Appends a SXdg gate (equivalent to Rx(-0.5,-)"
          " up to a -0.25 global phase)."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit"))
      .def(
          "Measure",
          [](Circuit *circ, unsigned qb, unsigned b, const py::kwargs &kwargs) {
            return add_gate_method_noparams<unsigned>(
                circ, OpType::Measure, {qb, b}, kwargs);
          },
          "Appends a single-qubit measurement in the computational "
          "(Z) basis."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit"), py::arg("bit_index"))
      .def(
          "Rz",
          [](Circuit *circ, const Expr &angle, unsigned qb,
             const py::kwargs &kwargs) {
            return add_gate_method_oneparam<unsigned>(
                circ, OpType::Rz, angle, {qb}, kwargs);
          },
          "Appends an Rz gate with a possibly symbolic angle "
          "(specified in half-turns)."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle"), py::arg("qubit"))
      .def(
          "Rx",
          [](Circuit *circ, const Expr &angle, unsigned qb,
             const py::kwargs &kwargs) {
            return add_gate_method_oneparam<unsigned>(
                circ, OpType::Rx, angle, {qb}, kwargs);
          },
          "Appends an Rx gate with a possibly symbolic angle "
          "(specified in half-turns)."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle"), py::arg("qubit"))
      .def(
          "Ry",
          [](Circuit *circ, const Expr &angle, unsigned qb,
             const py::kwargs &kwargs) {
            return add_gate_method_oneparam<unsigned>(
                circ, OpType::Ry, angle, {qb}, kwargs);
          },
          "Appends an Ry gate with a possibly symbolic angle "
          "(specified in half-turns)."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle"), py::arg("qubit"))
      .def(
          "U1",
          [](Circuit *circ, const Expr &angle, unsigned qb,
             const py::kwargs &kwargs) {
            return add_gate_method_oneparam<unsigned>(
                circ, OpType::U1, angle, {qb}, kwargs);
          },
          "Appends a U1 gate with a possibly symbolic angle "
          "(specified in half-turns)."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle"), py::arg("qubit"))
      .def(
          "U2",
          [](Circuit *circ, const Expr &angle0, const Expr &angle1,
             const unsigned &qb, const py::kwargs &kwargs) {
            return add_gate_method_manyparams<unsigned>(
                circ, OpType::U2, {angle0, angle1}, {qb}, kwargs);
          },
          "Appends a U2 gate with possibly symbolic angles "
          "(specified in half-turns)."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle0"), py::arg("angle1"), py::arg("qubit"))
      .def(
          "U3",
          [](Circuit *circ, const Expr &angle0, const Expr &angle1,
             const Expr &angle2, const unsigned &qb, const py::kwargs &kwargs) {
            return add_gate_method_manyparams<unsigned>(
                circ, OpType::U3, {angle0, angle1, angle2}, {qb}, kwargs);
          },
          "Appends a U3 gate with possibly symbolic angles "
          "(specified in half-turns)."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle0"), py::arg("angle1"), py::arg("angle2"),
          py::arg("qubit"))
      .def(
          "TK1",
          [](Circuit *circ, const Expr &angle0, const Expr &angle1,
             const Expr &angle2, const unsigned &qb, const py::kwargs &kwargs) {
            return add_gate_method_manyparams<unsigned>(
                circ, OpType::TK1, {angle0, angle1, angle2}, {qb}, kwargs);
          },
          "Appends a TK1 gate with possibly symbolic angles "
          "(specified in half-turns)."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle0"), py::arg("angle1"), py::arg("angle2"),
          py::arg("qubit"))
      .def(
          "TK2",
          [](Circuit *circ, const Expr &angle0, const Expr &angle1,
             const Expr &angle2, const unsigned &qb0, const unsigned &qb1,
             const py::kwargs &kwargs) {
            return add_gate_method_manyparams<unsigned>(
                circ, OpType::TK2, {angle0, angle1, angle2}, {qb0, qb1},
                kwargs);
          },
          "Appends a TK2 gate with possibly symbolic angles "
          "(specified in half-turns)."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle0"), py::arg("angle1"), py::arg("angle2"),
          py::arg("qubit0"), py::arg("qubit1"))
      .def(
          "CX",
          [](Circuit *circ, unsigned ctrl, unsigned trgt,
             const py::kwargs &kwargs) {
            return add_gate_method_noparams<unsigned>(
                circ, OpType::CX, {ctrl, trgt}, kwargs);
          },
          "Appends a CX gate on the wires for the specified control "
          "and target qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("control_qubit"), py::arg("target_qubit"))
      .def(
          "CY",
          [](Circuit *circ, unsigned ctrl, unsigned trgt,
             const py::kwargs &kwargs) {
            return add_gate_method_noparams<unsigned>(
                circ, OpType::CY, {ctrl, trgt}, kwargs);
          },
          "Appends a CY gate on the wires for the specified control "
          "and target qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("control_qubit"), py::arg("target_qubit"))
      .def(
          "CZ",
          [](Circuit *circ, unsigned ctrl, unsigned trgt,
             const py::kwargs &kwargs) {
            return add_gate_method_noparams<unsigned>(
                circ, OpType::CZ, {ctrl, trgt}, kwargs);
          },
          "Appends a CZ gate on the wires for the specified control "
          "and target qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("control_qubit"), py::arg("target_qubit"))
      .def(
          "CH",
          [](Circuit *circ, unsigned ctrl, unsigned trgt,
             const py::kwargs &kwargs) {
            return add_gate_method_noparams<unsigned>(
                circ, OpType::CH, {ctrl, trgt}, kwargs);
          },
          "Appends a CH gate on the wires for the specified control "
          "and target qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("control_qubit"), py::arg("target_qubit"))
      .def(
          "CV",
          [](Circuit *circ, unsigned ctrl, unsigned trgt,
             const py::kwargs &kwargs) {
            return add_gate_method_noparams<unsigned>(
                circ, OpType::CV, {ctrl, trgt}, kwargs);
          },
          "Appends a CV gate on the wires for the specified control "
          "and target qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("control_qubit"), py::arg("target_qubit"))
      .def(
          "CVdg",
          [](Circuit *circ, unsigned ctrl, unsigned trgt,
             const py::kwargs &kwargs) {
            return add_gate_method_noparams<unsigned>(
                circ, OpType::CVdg, {ctrl, trgt}, kwargs);
          },
          "Appends a CVdg gate on the wires for the specified control "
          "and target qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("control_qubit"), py::arg("target_qubit"))
      .def(
          "CSX",
          [](Circuit *circ, unsigned ctrl, unsigned trgt,
             const py::kwargs &kwargs) {
            return add_gate_method_noparams<unsigned>(
                circ, OpType::CSX, {ctrl, trgt}, kwargs);
          },
          "Appends a CSX gate on the wires for the specified control "
          "and target qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("control_qubit"), py::arg("target_qubit"))
      .def(
          "CSXdg",
          [](Circuit *circ, unsigned ctrl, unsigned trgt,
             const py::kwargs &kwargs) {
            return add_gate_method_noparams<unsigned>(
                circ, OpType::CSXdg, {ctrl, trgt}, kwargs);
          },
          "Appends a CSXdg gate on the wires for the specified control "
          "and target qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("control_qubit"), py::arg("target_qubit"))
      .def(
          "CRz",
          [](Circuit *circ, const Expr &angle, unsigned ctrl, unsigned trgt,
             const py::kwargs &kwargs) {
            return add_gate_method_oneparam<unsigned>(
                circ, OpType::CRz, angle, {ctrl, trgt}, kwargs);
          },
          "Appends a CRz gate with a possibly symbolic angle (specified in "
          "half-turns) on the wires for the specified control and "
          "target qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle"), py::arg("control_qubit"), py::arg("target_qubit"))
      .def(
          "CRx",
          [](Circuit *circ, const Expr &angle, unsigned ctrl, unsigned trgt,
             const py::kwargs &kwargs) {
            return add_gate_method_oneparam<unsigned>(
                circ, OpType::CRx, angle, {ctrl, trgt}, kwargs);
          },
          "Appends a CRx gate with a possibly symbolic angle (specified in "
          "half-turns) on the wires for the specified control and "
          "target qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle"), py::arg("control_qubit"), py::arg("target_qubit"))
      .def(
          "CRy",
          [](Circuit *circ, const Expr &angle, unsigned ctrl, unsigned trgt,
             const py::kwargs &kwargs) {
            return add_gate_method_oneparam<unsigned>(
                circ, OpType::CRy, angle, {ctrl, trgt}, kwargs);
          },
          "Appends a CRy gate with a possibly symbolic angle (specified in "
          "half-turns) on the wires for the specified control and "
          "target qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle"), py::arg("control_qubit"), py::arg("target_qubit"))
      .def(
          "CU1",
          [](Circuit *circ, const Expr &angle, unsigned ctrl, unsigned trgt,
             const py::kwargs &kwargs) {
            return add_gate_method_oneparam<unsigned>(
                circ, OpType::CU1, angle, {ctrl, trgt}, kwargs);
          },
          "Appends a CU1 gate with a possibly symbolic angle (specified in "
          "half-turns) on the wires for the specified control and "
          "target qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle"), py::arg("control_qubit"), py::arg("target_qubit"))
      .def(
          "CU3",
          [](Circuit *circ, const Expr &angle0, const Expr &angle1,
             const Expr &angle2, unsigned ctrl, unsigned trgt,
             const py::kwargs &kwargs) {
            return add_gate_method_manyparams<unsigned>(
                circ, OpType::CU3, {angle0, angle1, angle2}, {ctrl, trgt},
                kwargs);
          },
          "Appends a CU3 gate with possibly symbolic angles (specified in "
          "half-turns) on the wires for the specified control and "
          "target qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle0"), py::arg("angle1"), py::arg("angle2"),
          py::arg("control_qubit"), py::arg("target_qubit"))
      .def(
          "ZZPhase",
          [](Circuit *circ, const Expr &angle, unsigned qb0, unsigned qb1,
             const py::kwargs &kwargs) {
            return add_gate_method_oneparam<unsigned>(
                circ, OpType::ZZPhase, angle, {qb0, qb1}, kwargs);
          },
          "Appends a ZZ gate with a possibly symbolic angle (specified in "
          "half-turns) on the wires for the specified two qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle"), py::arg("qubit0"), py::arg("qubit1"))
      .def(
          "ZZMax",
          [](Circuit *circ, unsigned qb0, unsigned qb1,
             const py::kwargs &kwargs) {
            return add_gate_method_noparams<unsigned>(
                circ, OpType::ZZMax, {qb0, qb1}, kwargs);
          },
          "Appends a ZZMax gate on the wires for the specified two qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit0"), py::arg("qubit1"))
      .def(
          "ESWAP",
          [](Circuit *circ, const Expr &angle, unsigned qb0, unsigned qb1,
             const py::kwargs &kwargs) {
            return add_gate_method_oneparam<unsigned>(
                circ, OpType::ESWAP, angle, {qb0, qb1}, kwargs);
          },
          "Appends an ESWAP gate with a possibly symbolic angle (specified in "
          "half-turns) on the wires for the specified two qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle"), py::arg("qubit0"), py::arg("qubit1"))
      .def(
          "FSim",
          [](Circuit *circ, const Expr &angle0, const Expr &angle1,
             unsigned qb0, unsigned qb1, const py::kwargs &kwargs) {
            return add_gate_method_manyparams<unsigned>(
                circ, OpType::FSim, {angle0, angle1}, {qb0, qb1}, kwargs);
          },
          "Appends an FSim gate with possibly symbolic angles (specified in "
          "half-turns) on the wires for the specified qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle0"), py::arg("angle1"), py::arg("qubit0"),
          py::arg("qubit1"))
      .def(
          "Sycamore",
          [](Circuit *circ, unsigned qubit0, unsigned qubit1,
             const py::kwargs &kwargs) {
            return add_gate_method_noparams<unsigned>(
                circ, OpType::Sycamore, {qubit0, qubit1}, kwargs);
          },
          "Appends a Sycamore gate on the wires for the specified qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit0"), py::arg("qubit1"))
      .def(
          "XXPhase",
          [](Circuit *circ, const Expr &angle, unsigned qb0, unsigned qb1,
             const py::kwargs &kwargs) {
            return add_gate_method_oneparam<unsigned>(
                circ, OpType::XXPhase, angle, {qb0, qb1}, kwargs);
          },
          "Appends a XX gate with a possibly symbolic angle (specified in "
          "half-turns) on the wires for the specified two qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle"), py::arg("qubit0"), py::arg("qubit1"))
      .def(
          "YYPhase",
          [](Circuit *circ, const Expr &angle, unsigned qb0, unsigned qb1,
             const py::kwargs &kwargs) {
            return add_gate_method_oneparam<unsigned>(
                circ, OpType::YYPhase, angle, {qb0, qb1}, kwargs);
          },
          "Appends a YY gate with a possibly symbolic angle (specified in "
          "half-turns) on the wires for the specified two qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle"), py::arg("qubit0"), py::arg("qubit1"))
      .def(
          "XXPhase3",
          [](Circuit *circ, const Expr &angle, unsigned qb0, unsigned qb1,
             unsigned qb2, const py::kwargs &kwargs) {
            return add_gate_method_oneparam<unsigned>(
                circ, OpType::XXPhase3, angle, {qb0, qb1, qb2}, kwargs);
          },
          "Appends a 3-qubit XX gate with a possibly symbolic angle (specified "
          "in "
          "half-turns) on the wires for the specified three qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle"), py::arg("qubit0"), py::arg("qubit1"),
          py::arg("qubit2"))
      .def(
          "PhasedX",
          [](Circuit *circ, const Expr &angle0, const Expr &angle1, unsigned qb,
             const py::kwargs &kwargs) {
            return add_gate_method_manyparams<unsigned>(
                circ, OpType::PhasedX, {angle0, angle1}, {qb}, kwargs);
          },
          "Appends a PhasedX gate with possibly symbolic angles (specified in "
          "half-turns) on the wires for the specified qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle0"), py::arg("angle1"), py::arg("qubit"))
      .def(
          "CCX",
          [](Circuit *circ, unsigned ctrl1, unsigned ctrl2, unsigned trgt,
             const py::kwargs &kwargs) {
            return add_gate_method_noparams<unsigned>(
                circ, OpType::CCX, {ctrl1, ctrl2, trgt}, kwargs);
          },
          "Appends a CCX gate on the wires for the specified control "
          "and target qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("control_0"), py::arg("control_1"), py::arg("target"))
      .def(
          "ECR",
          [](Circuit *circ, unsigned qb0, unsigned qb1,
             const py::kwargs &kwargs) {
            return add_gate_method_noparams<unsigned>(
                circ, OpType::ECR, {qb0, qb1}, kwargs);
          },
          "Appends an ECR gate on the wires for the specified qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit_0"), py::arg("qubit_1"))
      .def(
          "SWAP",
          [](Circuit *circ, unsigned qb0, unsigned qb1,
             const py::kwargs &kwargs) {
            return add_gate_method_noparams<unsigned>(
                circ, OpType::SWAP, {qb0, qb1}, kwargs);
          },
          "Appends a SWAP gate on the wires for the specified qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit_0"), py::arg("qubit_1"))
      .def(
          "CSWAP",
          [](Circuit *circ, unsigned ctrl, unsigned trgt0, unsigned trgt1,
             const py::kwargs &kwargs) {
            return add_gate_method_noparams<unsigned>(
                circ, OpType::CSWAP, {ctrl, trgt0, trgt1}, kwargs);
          },
          "Appends a CSWAP gate on the wires for the specified "
          "control and target qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("control"), py::arg("target_0"), py::arg("target_1"))
      .def(
          "ISWAP",
          [](Circuit *circ, const Expr &angle, unsigned qb0, unsigned qb1,
             const py::kwargs &kwargs) {
            return add_gate_method_oneparam<unsigned>(
                circ, OpType::ISWAP, angle, {qb0, qb1}, kwargs);
          },
          "Appends an ISWAP gate with a possibly symbolic angle (specified in "
          "half-turns) on the wires for the specified qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle"), py::arg("qubit0"), py::arg("qubit1"))
      .def(
          "ISWAPMax",
          [](Circuit *circ, unsigned qb0, unsigned qb1,
             const py::kwargs &kwargs) {
            return add_gate_method_noparams<unsigned>(
                circ, OpType::ISWAPMax, {qb0, qb1}, kwargs);
          },
          "Appends an ISWAPMax gate on the wires for the specified qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit0"), py::arg("qubit1"))
      .def(
          "PhasedISWAP",
          [](Circuit *circ, const Expr &angle0, const Expr &angle1,
             unsigned qb0, unsigned qb1, const py::kwargs &kwargs) {
            return add_gate_method_manyparams<unsigned>(
                circ, OpType::PhasedISWAP, {angle0, angle1}, {qb0, qb1},
                kwargs);
          },
          "Appends a PhasedISWAP gate with possibly symbolic angles (specified "
          "in "
          "half-turns) on the wires for the specified qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle0"), py::arg("angle1"), py::arg("qubit0"),
          py::arg("qubit1"))
      .def(
          "measure_all",
          [](Circuit *circ) {
            opt_reg_info_t creg_info = circ->get_reg_info("c");
            register_info_t default_info = {UnitType::Bit, 1};
            if (creg_info && creg_info.value() != default_info)
              throw CircuitInvalidity(
                  "Cannot measure all; default classical "
                  "register name is already in use");
            qubit_vector_t qbs = circ->all_qubits();
            for (unsigned i = 0; i < qbs.size(); i++) {
              Bit id(i);
              circ->add_bit(id, false);
              circ->add_measure(qbs[i], id);
            }
            return circ;
          },
          "Appends a measure gate to all qubits, storing the results "
          "in the default classical register. Bits are added to the "
          "circuit if they do not already exist."
          "\n\n:return: the new :py:class:`Circuit`")
      .def(
          "measure_register",
          [](Circuit *circ, const QubitRegister &qreg,
             const std::string &creg_name) {
            if (!circ->get_reg_info(qreg.name())) {
              throw CircuitInvalidity(
                  "The given QubitRegister is not in use, please use "
                  "add_q_register to add it to the circuit first.");
            }
            opt_reg_info_t creg_info = circ->get_reg_info(creg_name);
            if (creg_info == std::nullopt) {
              circ->add_c_register(creg_name, qreg.size());
            } else if (circ->get_reg(creg_name).size() != qreg.size()) {
              throw CircuitInvalidity(
                  "The given classical register already exists, "
                  "but its size doesn't match the given QubitRegister.");
            }
            for (unsigned i = 0; i < qreg.size(); i++) {
              circ->add_measure(qreg[i], Bit(creg_name, i));
            }
            return circ;
          },
          "Appends a measure gate to all qubits in the given register, storing "
          "the results in the given classical register with matching indices."
          "The classical register will be created if it doesn't exist."
          "\n\n:param qreg: the QubitRegister to be measured"
          "\n:param creg_name: the name of the BitRegister to store the results"
          "\n:return: the new :py:class:`Circuit`")
      .def(
          "H",
          [](Circuit *circ, const Qubit &qb, const py::kwargs &kwargs) {
            return add_gate_method_noparams<UnitID>(
                circ, OpType::H, {qb}, kwargs);
          },
          "Appends a Hadamard gate."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit"))
      .def(
          "X",
          [](Circuit *circ, const Qubit &qb, const py::kwargs &kwargs) {
            return add_gate_method_noparams<UnitID>(
                circ, OpType::X, {qb}, kwargs);
          },
          "Appends an X gate."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit"))
      .def(
          "Y",
          [](Circuit *circ, const Qubit &qb, const py::kwargs &kwargs) {
            return add_gate_method_noparams<UnitID>(
                circ, OpType::Y, {qb}, kwargs);
          },
          "Appends a Y gate."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit"))
      .def(
          "Z",
          [](Circuit *circ, const Qubit &qb, const py::kwargs &kwargs) {
            return add_gate_method_noparams<UnitID>(
                circ, OpType::Z, {qb}, kwargs);
          },
          "Appends a Z gate."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit"))
      .def(
          "T",
          [](Circuit *circ, const Qubit &qb, const py::kwargs &kwargs) {
            return add_gate_method_noparams<UnitID>(
                circ, OpType::T, {qb}, kwargs);
          },
          "Appends a T gate (equivalent to Rz(0.25,-))."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit"))
      .def(
          "Tdg",
          [](Circuit *circ, const Qubit &qb, const py::kwargs &kwargs) {
            return add_gate_method_noparams<UnitID>(
                circ, OpType::Tdg, {qb}, kwargs);
          },
          "Appends a T-dagger gate (equivalent to Rz(-0.25,-))."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit"))
      .def(
          "S",
          [](Circuit *circ, const Qubit &qb, const py::kwargs &kwargs) {
            return add_gate_method_noparams<UnitID>(
                circ, OpType::S, {qb}, kwargs);
          },
          "Appends an S gate (equivalent to Rz(0.5,-))."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit"))
      .def(
          "Sdg",
          [](Circuit *circ, const Qubit &qb, const py::kwargs &kwargs) {
            return add_gate_method_noparams<UnitID>(
                circ, OpType::Sdg, {qb}, kwargs);
          },
          "Appends an S-dagger gate (equivalent to Rz(-0.5,-))."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit"))
      .def(
          "V",
          [](Circuit *circ, const Qubit &qb, const py::kwargs &kwargs) {
            return add_gate_method_noparams<UnitID>(
                circ, OpType::V, {qb}, kwargs);
          },
          "Appends a V gate (equivalent to Rx(0.5,-))."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit"))
      .def(
          "Vdg",
          [](Circuit *circ, const Qubit &qb, const py::kwargs &kwargs) {
            return add_gate_method_noparams<UnitID>(
                circ, OpType::Vdg, {qb}, kwargs);
          },
          "Appends a V-dagger gate (equivalent to Rx(-0.5,-))."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit"))
      .def(
          "SX",
          [](Circuit *circ, const Qubit &qb, const py::kwargs &kwargs) {
            return add_gate_method_noparams<UnitID>(
                circ, OpType::SX, {qb}, kwargs);
          },
          "Appends a SX gate (equivalent to Rx(0.5,-)"
          " up to a 0.25 global phase)."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit"))
      .def(
          "SXdg",
          [](Circuit *circ, const Qubit &qb, const py::kwargs &kwargs) {
            return add_gate_method_noparams<UnitID>(
                circ, OpType::SXdg, {qb}, kwargs);
          },
          "Appends a SXdg gate (equivalent to Rx(-0.5,-)"
          " up to a -0.25 global phase)."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit"))
      .def(
          "Measure",
          [](Circuit *circ, const Qubit &qb, const Bit &b,
             const py::kwargs &kwargs) {
            return add_gate_method_noparams<UnitID>(
                circ, OpType::Measure, {qb, b}, kwargs);
          },
          "Appends a single-qubit measurement in the computational "
          "(Z) basis."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit"), py::arg("bit"))
      .def(
          "Rz",
          [](Circuit *circ, const Expr &angle, const Qubit &qb,
             const py::kwargs &kwargs) {
            return add_gate_method_oneparam<UnitID>(
                circ, OpType::Rz, angle, {qb}, kwargs);
          },
          "Appends an Rz gate with a possibly symbolic angle "
          "(specified in half-turns)."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle"), py::arg("qubit"))
      .def(
          "Rx",
          [](Circuit *circ, const Expr &angle, const Qubit &qb,
             const py::kwargs &kwargs) {
            return add_gate_method_oneparam<UnitID>(
                circ, OpType::Rx, angle, {qb}, kwargs);
          },
          "Appends an Rx gate with a possibly symbolic angle "
          "(specified in half-turns)."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle"), py::arg("qubit"))
      .def(
          "Ry",
          [](Circuit *circ, const Expr &angle, const Qubit &qb,
             const py::kwargs &kwargs) {
            return add_gate_method_oneparam<UnitID>(
                circ, OpType::Ry, angle, {qb}, kwargs);
          },
          "Appends an Ry gate with a possibly symbolic angle "
          "(specified in half-turns)."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle"), py::arg("qubit"))
      .def(
          "U1",
          [](Circuit *circ, const Expr &angle, const Qubit &qb,
             const py::kwargs &kwargs) {
            return add_gate_method_oneparam<UnitID>(
                circ, OpType::U1, angle, {qb}, kwargs);
          },
          "Appends a U1 gate with a possibly symbolic angle "
          "(specified in half-turns)."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle"), py::arg("qubit"))
      .def(
          "U2",
          [](Circuit *circ, const Expr &angle0, const Expr &angle1,
             const Qubit &qb, const py::kwargs &kwargs) {
            return add_gate_method_manyparams<UnitID>(
                circ, OpType::U2, {angle0, angle1}, {qb}, kwargs);
          },
          "Appends a U2 gate with possibly symbolic angles "
          "(specified in half-turns)."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle0"), py::arg("angle1"), py::arg("qubit"))
      .def(
          "U3",
          [](Circuit *circ, const Expr &angle0, const Expr &angle1,
             const Expr &angle2, const Qubit &qb, const py::kwargs &kwargs) {
            return add_gate_method_manyparams<UnitID>(
                circ, OpType::U3, {angle0, angle1, angle2}, {qb}, kwargs);
          },
          "Appends a U3 gate with possibly symbolic angles "
          "(specified in half-turns)."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle0"), py::arg("angle1"), py::arg("angle2"),
          py::arg("qubit"))
      .def(
          "TK1",
          [](Circuit *circ, const Expr &angle0, const Expr &angle1,
             const Expr &angle2, const Qubit &qb, const py::kwargs &kwargs) {
            return add_gate_method_manyparams<UnitID>(
                circ, OpType::TK1, {angle0, angle1, angle2}, {qb}, kwargs);
          },
          "Appends a TK1 gate with possibly symbolic angles "
          "(specified in half-turns)."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle0"), py::arg("angle1"), py::arg("angle2"),
          py::arg("qubit"))
      .def(
          "TK2",
          [](Circuit *circ, const Expr &angle0, const Expr &angle1,
             const Expr &angle2, const Qubit &qb0, const Qubit &qb1,
             const py::kwargs &kwargs) {
            return add_gate_method_manyparams<UnitID>(
                circ, OpType::TK2, {angle0, angle1, angle2}, {qb0, qb1},
                kwargs);
          },
          "Appends a TK2 gate with possibly symbolic angles "
          "(specified in half-turns)."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle0"), py::arg("angle1"), py::arg("angle2"),
          py::arg("qubit0"), py::arg("qubit1"))
      .def(
          "CX",
          [](Circuit *circ, const Qubit &ctrl, const Qubit &trgt,
             const py::kwargs &kwargs) {
            return add_gate_method_noparams<UnitID>(
                circ, OpType::CX, {ctrl, trgt}, kwargs);
          },
          "Appends a CX gate on the wires for the specified control "
          "and target qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("control_qubit"), py::arg("target_qubit"))
      .def(
          "CY",
          [](Circuit *circ, const Qubit &ctrl, const Qubit &trgt,
             const py::kwargs &kwargs) {
            return add_gate_method_noparams<UnitID>(
                circ, OpType::CY, {ctrl, trgt}, kwargs);
          },
          "Appends a CY gate on the wires for the specified control "
          "and target qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("control_qubit"), py::arg("target_qubit"))
      .def(
          "CZ",
          [](Circuit *circ, const Qubit &ctrl, const Qubit &trgt,
             const py::kwargs &kwargs) {
            return add_gate_method_noparams<UnitID>(
                circ, OpType::CZ, {ctrl, trgt}, kwargs);
          },
          "Appends a CZ gate on the wires for the specified control "
          "and target qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("control_qubit"), py::arg("target_qubit"))
      .def(
          "CH",
          [](Circuit *circ, const Qubit &ctrl, const Qubit &trgt,
             const py::kwargs &kwargs) {
            return add_gate_method_noparams<UnitID>(
                circ, OpType::CH, {ctrl, trgt}, kwargs);
          },
          "Appends a CH gate on the wires for the specified control "
          "and target qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("control_qubit"), py::arg("target_qubit"))
      .def(
          "CV",
          [](Circuit *circ, const Qubit &ctrl, const Qubit &trgt,
             const py::kwargs &kwargs) {
            return add_gate_method_noparams<UnitID>(
                circ, OpType::CV, {ctrl, trgt}, kwargs);
          },
          "Appends a CV gate on the wires for the specified control "
          "and target qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("control_qubit"), py::arg("target_qubit"))
      .def(
          "CVdg",
          [](Circuit *circ, const Qubit &ctrl, const Qubit &trgt,
             const py::kwargs &kwargs) {
            return add_gate_method_noparams<UnitID>(
                circ, OpType::CVdg, {ctrl, trgt}, kwargs);
          },
          "Appends a CVdg gate on the wires for the specified control "
          "and target qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("control_qubit"), py::arg("target_qubit"))
      .def(
          "CSX",
          [](Circuit *circ, const Qubit &ctrl, const Qubit &trgt,
             const py::kwargs &kwargs) {
            return add_gate_method_noparams<UnitID>(
                circ, OpType::CSX, {ctrl, trgt}, kwargs);
          },
          "Appends a CSX gate on the wires for the specified control "
          "and target qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("control_qubit"), py::arg("target_qubit"))
      .def(
          "CSXdg",
          [](Circuit *circ, const Qubit &ctrl, const Qubit &trgt,
             const py::kwargs &kwargs) {
            return add_gate_method_noparams<UnitID>(
                circ, OpType::CSXdg, {ctrl, trgt}, kwargs);
          },
          "Appends a CSXdg gate on the wires for the specified control "
          "and target qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("control_qubit"), py::arg("target_qubit"))
      .def(
          "CRz",
          [](Circuit *circ, const Expr &angle, const Qubit &ctrl,
             const Qubit &trgt, const py::kwargs &kwargs) {
            return add_gate_method_oneparam<UnitID>(
                circ, OpType::CRz, angle, {ctrl, trgt}, kwargs);
          },
          "Appends a CRz gate with a symbolic angle (specified in "
          "half-turns) on the wires for the specified control and "
          "target qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle"), py::arg("control_qubit"), py::arg("target_qubit"))
      .def(
          "CRx",
          [](Circuit *circ, const Expr &angle, const Qubit &ctrl,
             const Qubit &trgt, const py::kwargs &kwargs) {
            return add_gate_method_oneparam<UnitID>(
                circ, OpType::CRx, angle, {ctrl, trgt}, kwargs);
          },
          "Appends a CRx gate with a symbolic angle (specified in "
          "half-turns) on the wires for the specified control and "
          "target qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle"), py::arg("control_qubit"), py::arg("target_qubit"))
      .def(
          "CRy",
          [](Circuit *circ, const Expr &angle, const Qubit &ctrl,
             const Qubit &trgt, const py::kwargs &kwargs) {
            return add_gate_method_oneparam<UnitID>(
                circ, OpType::CRy, angle, {ctrl, trgt}, kwargs);
          },
          "Appends a CRy gate with a symbolic angle (specified in "
          "half-turns) on the wires for the specified control and "
          "target qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle"), py::arg("control_qubit"), py::arg("target_qubit"))
      .def(
          "CU1",
          [](Circuit *circ, const Expr &angle, const Qubit &ctrl,
             const Qubit &trgt, const py::kwargs &kwargs) {
            return add_gate_method_oneparam<UnitID>(
                circ, OpType::CU1, angle, {ctrl, trgt}, kwargs);
          },
          "Appends a CU1 gate with a possibly symbolic angle (specified in "
          "half-turns) on the wires for the specified control and "
          "target qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle"), py::arg("control_qubit"), py::arg("target_qubit"))
      .def(
          "CU3",
          [](Circuit *circ, const Expr &angle0, const Expr &angle1,
             const Expr &angle2, const Qubit &ctrl, const Qubit &trgt,
             const py::kwargs &kwargs) {
            return add_gate_method_manyparams<UnitID>(
                circ, OpType::CU3, {angle0, angle1, angle2}, {ctrl, trgt},
                kwargs);
          },
          "Appends a CU3 gate with possibly symbolic angles (specified in "
          "half-turns) on the wires for the specified control and "
          "target qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle0"), py::arg("angle1"), py::arg("angle2"),
          py::arg("control_qubit"), py::arg("target_qubit"))

      .def(
          "ZZPhase",
          [](Circuit *circ, const Expr &angle, const Qubit &qb0,
             const Qubit &qb1, const py::kwargs &kwargs) {
            return add_gate_method_oneparam<UnitID>(
                circ, OpType::ZZPhase, angle, {qb0, qb1}, kwargs);
          },
          "Appends a ZZ gate with a symbolic angle (specified in "
          "half-turns) on the wires for the specified two qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle"), py::arg("qubit0"), py::arg("qubit1"))
      .def(
          "ZZMax",
          [](Circuit *circ, const Qubit &qb0, const Qubit &qb1,
             const py::kwargs &kwargs) {
            return add_gate_method_noparams<UnitID>(
                circ, OpType::ZZMax, {qb0, qb1}, kwargs);
          },
          "Appends a ZZMax gate on the wires for the specified two qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit0"), py::arg("qubit1"))
      .def(
          "ESWAP",
          [](Circuit *circ, const Expr &angle, const Qubit &qb0,
             const Qubit &qb1, const py::kwargs &kwargs) {
            return add_gate_method_oneparam<UnitID>(
                circ, OpType::ESWAP, angle, {qb0, qb1}, kwargs);
          },
          "Appends an ESWAP gate with a possibly symbolic angle (specified in "
          "half-turns) on the wires for the specified two qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle"), py::arg("qubit0"), py::arg("qubit1"))
      .def(
          "FSim",
          [](Circuit *circ, const Expr &angle0, const Expr &angle1,
             const Qubit &qb0, const Qubit &qb1, const py::kwargs &kwargs) {
            return add_gate_method_manyparams<UnitID>(
                circ, OpType::FSim, {angle0, angle1}, {qb0, qb1}, kwargs);
          },
          "Appends an FSim gate with possibly symbolic angles (specified in "
          "half-turns) on the wires for the specified qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle0"), py::arg("angle1"), py::arg("qubit0"),
          py::arg("qubit1"))
      .def(
          "Sycamore",
          [](Circuit *circ, const Qubit &qubit0, const Qubit &qubit1,
             const py::kwargs &kwargs) {
            return add_gate_method_noparams<UnitID>(
                circ, OpType::Sycamore, {qubit0, qubit1}, kwargs);
          },
          "Appends a Sycamore gate on the wires for the specified qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit0"), py::arg("qubit1"))
      .def(
          "XXPhase",
          [](Circuit *circ, const Expr &angle, const Qubit &qb0,
             const Qubit &qb1, const py::kwargs &kwargs) {
            return add_gate_method_oneparam<UnitID>(
                circ, OpType::XXPhase, angle, {qb0, qb1}, kwargs);
          },
          "Appends a XX gate with a symbolic angle (specified in "
          "half-turns) on the wires for the specified two qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit0"), py::arg("qubit1"), py::arg("angle"))
      .def(
          "YYPhase",
          [](Circuit *circ, const Expr &angle, const Qubit &qb0,
             const Qubit &qb1, const py::kwargs &kwargs) {
            return add_gate_method_oneparam<UnitID>(
                circ, OpType::YYPhase, angle, {qb0, qb1}, kwargs);
          },
          "Appends a YY gate with a symbolic angle (specified in "
          "half-turns) on the wires for the specified two qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit0"), py::arg("qubit1"), py::arg("angle"))
      .def(
          "XXPhase3",
          [](Circuit *circ, const Expr &angle, const Qubit &qb0,
             const Qubit &qb1, const Qubit &qb2, const py::kwargs &kwargs) {
            return add_gate_method_oneparam<UnitID>(
                circ, OpType::XXPhase3, angle, {qb0, qb1, qb2}, kwargs);
          },
          "Appends a 3-qubit XX gate with a symbolic angle (specified in "
          "half-turns) on the wires for the specified three qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle"), py::arg("qubit0"), py::arg("qubit1"),
          py::arg("qubit2"))
      .def(
          "PhasedX",
          [](Circuit *circ, const Expr &angle0, const Expr &angle1,
             const Qubit &qubit, const py::kwargs &kwargs) {
            return add_gate_method_manyparams<UnitID>(
                circ, OpType::PhasedX, {angle0, angle1}, {qubit}, kwargs);
          },
          "Appends a PhasedX gate with possibly symbolic angles (specified in "
          "half-turns) on the wires for the specified qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle0"), py::arg("angle1"), py::arg("qubit"))
      .def(
          "CCX",
          [](Circuit *circ, const Qubit &ctrl1, const Qubit &ctrl2,
             const Qubit &trgt, const py::kwargs &kwargs) {
            return add_gate_method_noparams<UnitID>(
                circ, OpType::CCX, {ctrl1, ctrl2, trgt}, kwargs);
          },
          "Appends a CCX gate on the wires for the specified control "
          "and target qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("control_0"), py::arg("control_1"), py::arg("target"))
      .def(
          "ECR",
          [](Circuit *circ, const Qubit &qb1, const Qubit &qb2,
             const py::kwargs &kwargs) {
            return add_gate_method_noparams<UnitID>(
                circ, OpType::ECR, {qb1, qb2}, kwargs);
          },
          "Appends an ECR gate on the wires for the specified qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit_0"), py::arg("qubit_1"))
      .def(
          "SWAP",
          [](Circuit *circ, const Qubit &qb1, const Qubit &qb2,
             const py::kwargs &kwargs) {
            return add_gate_method_noparams<UnitID>(
                circ, OpType::SWAP, {qb1, qb2}, kwargs);
          },
          "Appends a SWAP gate on the wires for the specified qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit_0"), py::arg("qubit_1"))
      .def(
          "CSWAP",
          [](Circuit *circ, const Qubit &ctrl, const Qubit &trgt1,
             const Qubit &trgt2, const py::kwargs &kwargs) {
            return add_gate_method_noparams<UnitID>(
                circ, OpType::CSWAP, {ctrl, trgt1, trgt2}, kwargs);
          },
          "Appends a CSWAP gate on the wires for the specified "
          "control and target qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("control"), py::arg("target_0"), py::arg("target_1"))
      .def(
          "ISWAP",
          [](Circuit *circ, const Expr &angle, const Qubit &qb0,
             const Qubit &qb1, const py::kwargs &kwargs) {
            return add_gate_method_oneparam<UnitID>(
                circ, OpType::ISWAP, angle, {qb0, qb1}, kwargs);
          },
          "Appends an ISWAP gate with a possibly symbolic angle (specified in "
          "half-turns) on the wires for the specified qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle"), py::arg("qubit0"), py::arg("qubit1"))
      .def(
          "ISWAPMax",
          [](Circuit *circ, const Qubit &qb0, const Qubit &qb1,
             const py::kwargs &kwargs) {
            return add_gate_method_noparams<UnitID>(
                circ, OpType::ISWAPMax, {qb0, qb1}, kwargs);
          },
          "Appends an ISWAPMax gate on the wires for the specified qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("qubit0"), py::arg("qubit1"))
      .def(
          "PhasedISWAP",
          [](Circuit *circ, const Expr &angle0, const Expr &angle1,
             const Qubit &qb0, const Qubit &qb1, const py::kwargs &kwargs) {
            return add_gate_method_manyparams<UnitID>(
                circ, OpType::PhasedISWAP, {angle0, angle1}, {qb0, qb1},
                kwargs);
          },
          "Appends a PhasedISWAP gate with posisbly symbolic angles (specified "
          "in "
          "half-turns) on the wires for the specified qubits."
          "\n\n:return: the new :py:class:`Circuit`",
          py::arg("angle0"), py::arg("angle1"), py::arg("qubit0"),
          py::arg("qubit1"))
      .def(
          "Phase",
          [](Circuit *circ, const Expr &angle, const py::kwargs &kwargs) {
            return add_gate_method_oneparam<UnitID>(
                circ, OpType::Phase, angle, {}, kwargs);
          });
}

}  // namespace tket
