// Copyright 2019-2021 Cambridge Quantum Computing
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

#include "Program/Program.hpp"

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "typecast.hpp"

namespace py = pybind11;

namespace tket {

template <class ID>
static Program *add_gate_method(
    Program *prog, OpType type, const std::vector<Expr> &params,
    const std::vector<ID> &args, const py::kwargs &kwargs) {
  if (kwargs.contains("condition_bits")) {
    std::vector<ID> bits = py::cast<std::vector<ID>>(kwargs["condition_bits"]);
    unsigned n_args = args.size(), n_bits = bits.size();
    unsigned value = kwargs.contains("condition_value")
                         ? py::cast<unsigned>(kwargs["condition_value"])
                         : (1u << n_bits) - 1;
    Op_ptr con = std::make_shared<Conditional>(
        get_op_ptr(type, params, n_args), n_bits, value);
    std::vector<ID> new_args = bits;
    new_args.insert(new_args.end(), args.begin(), args.end());
    prog->add_op(con, new_args);
  } else {
    prog->add_op<ID>(type, params, args);
  }
  return prog;
}

PYBIND11_MODULE(program, m) {
  py::class_<Program>(
      m, "Program",
      "Encapsulates a control flow graph for a quantum program. Each "
      "basic block is a single quantum circuit which may include "
      "classical instructions and OpenQASM-style conditional gates. "
      "Branches are always made using a single condition bit. Allows "
      "long sequences of operations to be applied conditionally or "
      "repeatedly while some bit is true.")
      .def(py::init<>(), "Constructs an empty program.")
      .def(
          py::init<unsigned, unsigned>(),
          "Constructs a program with a given number of quantum and "
          "classical bits\n\n:param n_qubits: The number of qubits in "
          "the program\n:param c_bits: The number of classical bits in "
          "the program",
          py::arg("n_qubits"), py::arg("n_bits") = 0)
      .def("__str__", [](const Program &) { return "<tket::Program>"; })
      .def(
          "__repr__",
          [](const Program &prog) {
            std::stringstream ss;
            ss << "[";
            for (auto com : prog) {
              ss << com.to_str() << " ";
            }
            ss << "]";
            return ss.str();
          })
      .def(
          "__iter__",
          [](const Program &prog) {
            return py::make_iterator(prog.begin(), prog.end());
          },
          "Iterate through the program, a Command at a time.",
          py::keep_alive<
              0, 1>() /* Essential: keep object alive while iterator exists */)
      .def(
          "add_q_register", &Program::add_q_register,
          "Constructs a new quantum register with a given name and "
          "number of qubits.\n\n:param name: Unique readable name for "
          "the register\n:param size: Number of qubits "
          "required\n:return: a map from index to the corresponding "
          "UnitIDs",
          py::arg("name"), py::arg("size"))
      .def(
          "add_c_register", &Program::add_c_register,
          "Constructs a new classical register with a given name and "
          "number of bits.\n\n:param name: Unique readable name for the "
          "register\n:param size: Number of bits required\n:return: a "
          "map from index to the corresponding UnitIDs",
          py::arg("name"), py::arg("size"))
      .def(
          "add_qubit", &Program::add_qubit,
          "Constructs a single qubit with the given id.\n\n:param id: "
          "Unique id for the qubit\n:param reject_dups: Fail if there "
          "is already a qubit in this program with the id. Default to "
          "True",
          py::arg("id"), py::arg("reject_dups") = true)
      .def(
          "add_bit", &Program::add_bit,
          "Constructs a single bit with the given id.\n\n:param id: "
          "Unique id for the bit\n:param reject_dups: Fail if there is "
          "already a bit in this program with the id. Default to True",
          py::arg("id"), py::arg("reject_dups") = true)
      .def_property_readonly(
          "qubits", &Program::all_qubits,
          "A list of all qubit ids in the program")
      .def_property_readonly(
          "bits", &Program::all_bits,
          "A list of all classical bit ids in the program")
      .def_property_readonly(
          "bit_readout", &Program::bit_readout,
          "A map from bit to its (left-to-right) index in readouts "
          "from backends (following the increasing lexicographic "
          "order convention)")
      .def_property_readonly(
          "qubit_readout", &Program::qubit_readout,
          "A map from qubit to its (left-to-right) index in readouts "
          "from backends")
      .def(
          "get_commands",
          [](const Program &prog) {
            std::vector<Command> out;
            for (Command c : prog) out.push_back(c);
            return out;
          },
          ":return: a list of all the Commands in the program")
      .def(
          "add_gate",
          [](Program *prog, OpType type, const std::vector<unsigned> &args,
             const py::kwargs &kwargs) {
            return add_gate_method(prog, type, {}, args, kwargs);
          },
          "Appends a single (non-parameterised) gate to the end of "
          "the program on some particular qubits from the default "
          "register ('q'). The number of qubits specified must match "
          "the arity of the gate."
          "\n\n:param type: The type of operation to add"
          "\n:param args: The list of indices for the qubits/bits to "
          "which the operation is applied"
          "\n:param kwargs: Additional properties for classical "
          "conditions"
          "\n:return: the new :py:class:`Program`",
          py::arg("type"), py::arg("args"))
      .def(
          "add_gate",
          [](Program *prog, OpType type, const Expr &p,
             const std::vector<unsigned> &args, const py::kwargs &kwargs) {
            return add_gate_method(prog, type, {p}, args, kwargs);
          },
          "Appends a single gate, parameterised by an expression, to "
          "the end of the program on some particular qubits from the "
          "default register ('q')."
          "\n\n:param type: The type of gate to add"
          "\n:param angle: The parameter for the gate in halfturns"
          "\n:param args: The list of indices for the qubits/bits to "
          "which the operation is applied"
          "\n:param kwargs: Additional properties for classical "
          "conditions"
          "\n:return: the new :py:class:`Program`",
          py::arg("type"), py::arg("angle"), py::arg("args"))
      .def(
          "add_gate", &add_gate_method<unsigned>,
          "Appends a single gate, parameterised with a vector of "
          "expressions corresponding to halfturns, to the end of the "
          "program on some particular qubits from the default register "
          "('q')."
          "\n\n:param type: The type of gate to add"
          "\n:param angles: The parameters for the gate in halfturns"
          "\n:param args: The list of indices for the qubits/bits to "
          "which the operation is applied"
          "\n:param kwargs: Additional properties for classical "
          "conditions"
          "\n:return: the new :py:class:`Program`",
          py::arg("type"), py::arg("angles"), py::arg("args"))
      .def(
          "add_gate", &add_gate_method<UnitID>,
          "Appends a single gate to the end of the program"
          "\n\n:param type: The type of gate to add"
          "\n:param params: The parameters for the gate in halfturns"
          "\n:param args: The qubits/bits to apply the gate to"
          "\n:param kwargs: Additional properties for classical "
          "conditions"
          "\n:return: the new :py:class:`Program`",
          py::arg("type"), py::arg("params"), py::arg("args"))
      .def(
          "append_circuit",
          [](Program *prog, const Circuit &circ) {
            prog->add_block(circ);
            return prog;
          },
          "Appends a circuit to the end of the program"
          "\n\n:param circuit: The circuit to add"
          "\n:return: the new :py:class:`Program`",
          py::arg("circuit"))
      .def(
          "append", &Program::append,
          "In-place sequential composition of programs, appending a "
          "copy of the argument onto the end of `self`."
          "\n\n:param prog: The program to be appended to the end of "
          "`self`",
          py::arg("prog"))
      .def(
          "append_if", &Program::append_if,
          "In-place sequential composition of programs, performing "
          "`body` after `self` if the `condition_bit` is found to be 1."
          "\n\n:param condition_bit: A single bit condition."
          "\n\n:param body: The program to be applied after `self` if "
          "`condition_bit` is 1.",
          py::arg("condition_bit"), py::arg("body"))
      .def(
          "append_if_else", &Program::append_if_else,
          "In-place sequential composition of programs, performing "
          "`if_body` after `self` if the `condition_bit` is found to be "
          "1, and `else_body` if it is 0."
          "\n\n:param condition_bit: A single bit condition."
          "\n\n:param if_body: The program to be applied after `self` "
          "if `condition_bit` is 1.",
          "\n\n:param else_body: The program to be applied after `self` "
          "if `condition_bit` is 0.",
          py::arg("condition_bit"), py::arg("if_body"), py::arg("else_body"))
      .def(
          "append_while", &Program::append_while,
          "In-place sequential composition of programs, performing "
          "`body` after `self` repeatedly whilst the `condition_bit` is "
          "found to be 1."
          "\n\n:param condition_bit: A single bit condition."
          "\n\n:param body: The program to be applied after `self` "
          "repeatedly whilst `condition_bit` is 1.",
          py::arg("condition_bit"), py::arg("body"));
}

}  // namespace tket
