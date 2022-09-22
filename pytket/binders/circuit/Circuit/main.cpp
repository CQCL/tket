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
#include <pybind11/functional.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <optional>
#include <sstream>

#include "Circuit/Boxes.hpp"
#include "Circuit/Circuit.hpp"
#include "Circuit/Command.hpp"
#include "Gate/OpPtrFunctions.hpp"
#include "Gate/SymTable.hpp"
#include "Mapping/Verification.hpp"
#include "Ops/Op.hpp"
#include "Simulation/CircuitSimulator.hpp"
#include "UnitRegister.hpp"
#include "Utils/Json.hpp"
#include "binder_json.hpp"
#include "binder_utils.hpp"
#include "boost/graph/iteration_macros.hpp"
#include "typecast.hpp"

namespace py = pybind11;
using json = nlohmann::json;

namespace tket {

const bit_vector_t no_bits;

void init_circuit_add_op(py::class_<Circuit, std::shared_ptr<Circuit>> &c);
void init_circuit_add_classical_op(
    py::class_<Circuit, std::shared_ptr<Circuit>> &c);

void init_circuit(py::module &m) {
  py::class_<Circuit, std::shared_ptr<Circuit>> circuit_cls(
      m, "Circuit", py::dynamic_attr(),
      "Encapsulates a quantum circuit using a DAG representation.\n\n>>> "
      "from pytket import Circuit\n>>> c = Circuit(4,2) # Create a circuit "
      "with 4 qubits and 2 classical bits"
      "\n>>> c.H(0) # Apply a gate to qubit 0\n>>> "
      "c.Rx(0.5,1) # Angles of rotation are expressed in half-turns "
      "(i.e. 0.5 means PI/2)\n>>> c.Measure(1,0) # Measure qubit 1, saving "
      "result in bit 0");
  init_circuit_add_op(circuit_cls);
  init_circuit_add_classical_op(circuit_cls);
  circuit_cls
      .def(py::init<>(), "Constructs a circuit with a completely empty DAG.")
      .def(
          py::init<const std::string &>(),
          "Constructs a named circuit with a completely empty DAG."
          "\n\n:param name: name for the circuit",
          py::arg("name"))
      .def(
          py::init<unsigned, std::optional<std::string>>(),
          "Constructs a circuit with a given number of qubits/blank "
          "wires.\n\n>>> c = Circuit()\n>>> c.add_blank_wires(3)\n\nis "
          "equivalent to\n\n>>> c = Circuit(3)\n\n:param n_qubits: The "
          "number of qubits in the circuit\n:param name: Optional name "
          "for the circuit.",
          py::arg("n_qubits"), py::arg("name") = std::nullopt)
      .def(
          py::init<unsigned, unsigned, std::optional<std::string>>(),
          "Constructs a circuit with a given number of quantum and "
          "classical bits\n\n:param n_qubits: The number of qubits in "
          "the circuit\n:param n_bits: The number of classical bits in "
          "the circuit\n:param name: Optional name for the circuit.",
          py::arg("n_qubits"), py::arg("n_bits"),
          py::arg("name") = std::nullopt)
      .def("__eq__", &Circuit::operator==)
      .def(
          "__str__",
          [](const Circuit &circ) {
            return "<tket::Circuit, qubits=" + std::to_string(circ.n_qubits()) +
                   ", gates=" + std::to_string(circ.n_gates()) + ">";
          })
      .def(
          "__repr__",
          [](const Circuit &circ) {
            std::stringstream ss;
            ss << "[";
            for (auto com : circ.get_commands()) {
              ss << com.to_str() << " ";
            }
            ss << "]";
            return ss.str();
          })
      .def(
          "__iter__",
          [](const Circuit &circ) {
            return py::make_iterator(circ.begin(), circ.end());
          },
          "Iterate through the circuit, a Command at a time.",
          py::keep_alive<
              0, 1>() /* Essential: keep object alive while iterator exists */)
      .def(
          "get_commands",
          [](const Circuit &circ) {
            std::vector<Command> out;
            for (Command c : circ) out.push_back(c);
            return out;
          },
          ":return: a list of all the Commands in the circuit")
      .def(
          "get_unitary",
          [](const Circuit &circ) { return tket_sim::get_unitary(circ); },
          ":return: The numerical unitary matrix of the circuit, using ILO-BE "
          "convention.")
      .def(
          "get_unitary_times_other",
          [](const Circuit &circ, Eigen::MatrixXcd matr) {
            tket_sim::apply_unitary(circ, matr);
            return matr;
          },
          "Calculate UM, where U is the numerical unitary matrix of "
          "the circuit, with ILO-BE convention, and M is another matrix. "
          "This is more efficient than calculating U separately, if M has "
          "fewer columns than U."
          "\n\n:param matr: The matrix to be multiplied."
          "\n:return: The product of the circuit unitary and the given matrix.",
          py::arg("matr"))
      .def(
          "get_statevector",
          [](const Circuit &circ) { return tket_sim::get_statevector(circ); },
          "Calculate the unitary matrix of the circuit, using ILO-BE "
          "convention, applied to the column vector (1,0,0...), "
          "which is thus another column vector. Due to "
          "pybind11 and numpy peculiarities, to treat the "
          "result as a genuine column vector and perform further "
          "matrix multiplication, you need to call "
          ".reshape(rows,1) to get a 2D matrix with "
          "the correct dimensions."
          "\n\n:return: The calculated vector.")
      .def(
          "add_q_register",
          [](Circuit &circ, const std::string &name, const std::size_t &size) {
            circ.add_q_register(name, size);
            return QubitRegister(name, size);
          },
          "Constructs a new quantum register with a given name and "
          "number of qubits.\n\n:param name: Unique readable name for "
          "the register\n:param size: Number of qubits "
          "required\n:return: a map from index to the corresponding "
          "UnitIDs",
          py::arg("name"), py::arg("size"))
      .def(
          "add_q_register",
          [](Circuit &circ, const QubitRegister &reg) {
            const std::string &name = reg.name();
            const std::size_t &size = reg.size();
            register_t existing = circ.get_reg(name);

            if (existing.size() > 0) {
              if (existing.size() != size) {
                throw CircuitInvalidity(
                    "Existing register with name \"" + name +
                    "\" already exists, and does not match "
                    "requested size.");
              }
              return reg;
            }
            circ.add_q_register(name, size);
            return reg;
          },
          "Adds QubitRegister to Circuit"
          "\n\n:param register: QubitRegister ",
          py::arg("register"))
      .def(
          "add_c_register",
          [](Circuit &circ, const std::string &name, const std::size_t &size) {
            circ.add_c_register(name, size);
            return BitRegister(name, size);
          },
          "Constructs a new classical register with a given name and "
          "number of bits.\n\n:param name: Unique readable name for the "
          "register\n:param size: Number of bits required\n:return: a "
          "map from index to the corresponding UnitIDs",
          py::arg("name"), py::arg("size"))
      .def(
          "add_c_register",
          [](Circuit &circ, const BitRegister &reg) {
            const std::string &name = reg.name();
            const std::size_t &size = reg.size();
            register_t existing = circ.get_reg(name);

            if (existing.size() > 0) {
              if (existing.size() != size) {
                throw CircuitInvalidity(
                    "Existing register with name \"" + name +
                    "\" already exists, and does not match "
                    "requested size.");
              }
              return reg;
            }
            circ.add_c_register(name, size);
            return reg;
          },
          "Adds BitRegister to Circuit"
          "\n\n:param register: BitRegister ",
          py::arg("register"))
      .def(
          "get_c_register",
          [](Circuit &circ, const std::string &name) {
            register_t reg = circ.get_reg(name);
            if (reg.size() == 0 ||
                reg.begin()->second.type() != UnitType::Bit) {
              throw CircuitInvalidity(
                  "Cannot find classical register with name \"" + name + "\".");
            }
            return BitRegister(name, reg.size());
          },
          "Get the classical register with the given name.\n\n:param name: "
          "name for the register\n:return: the retrieved "
          ":py:class:`BitRegister`",
          py::arg("name"))
      .def_property_readonly(
          "c_registers",
          [](Circuit &circ) {
            bit_vector_t all_bits = circ.all_bits();
            std::map<std::string, unsigned> bits_map;
            std::vector<BitRegister> b_regs;
            for (Bit bit : all_bits) {
              auto it = bits_map.find(bit.reg_name());
              if (it == bits_map.end()) {
                bits_map.insert({bit.reg_name(), 1});
              } else {
                it->second++;
              }
            }
            for (auto const &it : bits_map) {
              b_regs.push_back(BitRegister(it.first, it.second));
            }
            return b_regs;
          },
          "Get all classical registers.\n\n:return: List of "
          ":py:class:`BitRegister`")
      .def(
          "get_q_register",
          [](Circuit &circ, const std::string &name) {
            register_t reg = circ.get_reg(name);
            if (reg.size() == 0 ||
                reg.begin()->second.type() != UnitType::Qubit) {
              throw CircuitInvalidity(
                  "Cannot find quantum register with name \"" + name + "\".");
            }
            return QubitRegister(name, reg.size());
          },
          "Get the quantum register with the given name.\n\n:param name: "
          "name for the register\n:return: the retrieved "
          ":py:class:`QubitRegister`",
          py::arg("name"))
      .def_property_readonly(
          "q_registers",
          [](Circuit &circ) {
            qubit_vector_t all_qbs = circ.all_qubits();
            std::map<std::string, unsigned> qbs_map;
            std::vector<QubitRegister> q_regs;
            for (Qubit qb : all_qbs) {
              auto it = qbs_map.find(qb.reg_name());
              if (it == qbs_map.end()) {
                qbs_map.insert({qb.reg_name(), 1});
              } else {
                it->second++;
              }
            }
            for (auto const &it : qbs_map) {
              q_regs.push_back(QubitRegister(it.first, it.second));
            }
            return q_regs;
          },
          "Get all quantum registers.\n\n:return: List of "
          ":py:class:`QubitRegister`")
      .def(
          "add_qubit", &Circuit::add_qubit,
          "Constructs a single qubit with the given id.\n\n:param id: "
          "Unique id for the qubit\n:param reject_dups: Fail if there "
          "is already a qubit in this circuit with the id. Default to "
          "True",
          py::arg("id"), py::arg("reject_dups") = true)
      .def(
          "add_bit", &Circuit::add_bit,
          "Constructs a single bit with the given id.\n\n:param id: "
          "Unique id for the bit\n:param reject_dups: Fail if there is "
          "already a bit in this circuit with the id. Default to True",
          py::arg("id"), py::arg("reject_dups") = true)
      .def_property_readonly(
          "qubits", &Circuit::all_qubits,
          "A list of all qubit ids in the circuit")
      .def_property_readonly(
          "bits", &Circuit::all_bits,
          "A list of all classical bit ids in the circuit")
      .def_property_readonly(
          "bit_readout", &Circuit::bit_readout,
          "A map from bit to its (left-to-right) index in readouts "
          "from backends (following the increasing lexicographic "
          "order convention)")
      .def_property_readonly(
          "qubit_readout", &Circuit::qubit_readout,
          "A map from qubit to its (left-to-right) index in readouts "
          "from backends. A qubit will feature in this map if it is "
          "measured and neither it nor the bit containing the "
          "measurement result is subsequently acted on")
      .def_property_readonly(
          "qubit_to_bit_map", &Circuit::qubit_to_bit_map,
          "A map from qubit to the bit it is measured to. "
          "A qubit will feature in this map if it is "
          "measured and neither it nor the bit containing the "
          "measurement result is subsequently acted on")
      .def_property_readonly(
          "opgroups", &Circuit::get_opgroups,
          "A set of all opgroup names in the circuit")
      .def(
          "flatten_registers", &Circuit::flatten_registers,
          "Combines all qubits into a single register namespace with "
          "the default name, and likewise for bits")

      // Circuit composition:
      .def(
          "add_circuit",
          [](Circuit &circ, const Circuit &circ2, const std::vector<Qubit> &qbs,
             const std::vector<Bit> &bits) {
            unit_map_t umap;
            unsigned i = 0;
            for (const Qubit &q : qbs) {
              umap.insert({Qubit(i), q});
              ++i;
            }
            i = 0;
            for (const Bit &b : bits) {
              umap.insert({Bit(i), b});
              ++i;
            }
            circ.append_with_map(circ2, umap);
            return &circ;
          },
          "In-place sequential composition of circuits, appending a "
          "copy of the argument onto the end of the circuit. "
          "Connects qubits and bits with the same behaviour as "
          ":py:meth:`add_gate`."
          "\n\n:param circuit: The circuit to be appended to the end "
          "of `self`"
          "\n:param qubits: List mapping the (default register) "
          "qubits of `circuit` to the qubits of `self`"
          "\n:param bits: List mapping the (default register) bits "
          "of `circuit` to the bits of `self`"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("circuit"), py::arg("qubits"), py::arg("bits") = no_bits)
      .def(
          "add_circuit",
          [](Circuit &circ, const Circuit &circ2,
             const std::vector<unsigned> &qbs,
             const std::vector<unsigned> &bits) {
            circ.append_qubits(circ2, qbs, bits);
            return &circ;
          },
          "In-place sequential composition of circuits, appending a "
          "copy of the argument onto the end of the circuit. "
          "Connects qubits and bits with the same behaviour as "
          ":py:meth:`add_gate`."
          "\n\n:param circuit: The circuit to be appended to the end "
          "of `self`"
          "\n:param qubits: List mapping the (default register) "
          "qubits of `circuit` to the (default register) qubits of "
          "`self`"
          "\n:param bits: List mapping the (default register) bits "
          "of `circuit` to the (default register) bits of `self`"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("circuit"), py::arg("qubits"), py::arg("bits") = no_bits)
      .def(
          "append", (void(Circuit::*)(const Circuit &)) & Circuit::append,
          "In-place sequential composition of circuits, appending a "
          "copy of the argument onto the end of the circuit. Inputs and "
          "Outputs are unified if they share the same id, defaulting to "
          "parallel composition if there is no match."
          "\n\n:param circuit: The circuit to be appended to the end of "
          "`self`",
          py::arg("circuit"))
      .def(
          "add_phase",
          [](Circuit &circ, Expr a) {
            circ.add_phase(a);
            return &circ;
          },
          "Add a global phase to the circuit.\n\n:param a: Phase to "
          "add, in halfturns\n\n:return: circuit with added phase",
          py::arg("a"))
      .def(
          "_n_vertices", &Circuit::n_vertices,
          ":return: the number of vertices in the DAG, i.e. the sum of "
          "the number of operations, inputs, and outputs")
      .def_property_readonly(
          "is_simple", [](const Circuit &circ) { return circ.is_simple(); },
          "Checks that the circuit has only 1 quantum and 1 classic "
          "register using the default names ('q' and 'c'). This "
          "means it is suitable to refer to qubits simply by their "
          "integer indices.")
      .def_property_readonly(
          "n_gates", &Circuit::n_gates,
          ":return: the number of gates in the Circuit")
      .def_property_readonly(
          "n_qubits", &Circuit::n_qubits,
          ":return: the number of qubits in the circuit")
      .def_property_readonly(
          "n_bits", &Circuit::n_bits,
          ":return: the number of classiclal bits in the circuit")
      .def_property_readonly(
          "phase", &Circuit::get_phase,
          ":return: the global phase applied to the circuit, in "
          "halfturns (not meaningful for circuits with classical "
          "interactions)")
      .def_property("name", &Circuit::get_name, &Circuit::set_name)
      .def(
          "remove_blank_wires", &Circuit::remove_blank_wires,
          "Removes any Input-Output pairs in the DAG with no "
          "intervening operations, i.e. removes untouched qubits/bits "
          "from the circuit. This may occur when optimisations "
          "recognise that the operations on a qubit reduce to the "
          "identity, or when routing adds wires to \"fill out\" the "
          "architecture.")
      .def(
          "add_blank_wires", &Circuit::add_blank_wires,
          "Adds a number of new qubits to the circuit. These will be "
          "added to the default register ('q') if possible, filling out "
          "the unused indices from 0.\n\n:param number: Number of "
          "qubits to add",
          py::arg("number"))
      .def(
          "rename_units", &Circuit::rename_units<UnitID, UnitID>,
          "Rename qubits and bits simultaneously according to the map "
          "of ids provided\n\n:param map: Dictionary from current ids "
          "to new ids",
          py::arg("map"))
      .def(
          "depth", [](const Circuit &circ) { return circ.depth(); },
          "Returns the number of interior vertices on the longest path through "
          "the DAG, excluding vertices representing barrier operations."
          "\n\n>>> c = Circuit(3)"
          "\n>>> c.depth()"
          "\n0"
          "\n>>> c.CX(0,1)"
          "\n>>> c.CX(1,2)"
          "\n>>> c.CX(2,0)"
          "\n>>> c.depth()"
          "\n3"
          "\n\n:return: the circuit depth")
      .def(
          "n_gates_of_type",
          [](const Circuit &circ, const OpType &_type) {
            return circ.count_gates(_type);
          },
          "Returns the number of vertices in the dag of a given "
          "operation type.\n\n>>> c.CX(0,1)\n>>> c.H(0)\n>>> "
          "c.CX(0,1)\n>>> c.n_gates_of_type(OpType.CX)\n2\n\n:param "
          "type: The operation type to search for\n:return: the "
          "number of operations matching `type`",
          py::arg("type"))
      .def(
          "depth_by_type", &Circuit::depth_by_type,
          "Returns the number of vertices in the longest path through the "
          "sub-DAG consisting of vertices representing operations of the given "
          "type."
          "\n\n>>> c = Circuit(3)"
          "\n>>> c.CX(0,1)"
          "\n>>> c.Z(1)"
          "\n>>> c.CX(1,2)"
          "\n>>> c.depth_by_type(OpType.CX)"
          "\n2"
          "\n\n:param type: the operation type of interest"
          "\n:return: the circuit depth with respect to operations matching "
          "`type`",
          py::arg("type"))
      .def(
          "depth_by_type", &Circuit::depth_by_types,
          "Returns the number of vertices in the longest path through the "
          "sub-DAG consisting of vertices representing operations of the given "
          "types."
          "\n\n>>> c = Circuit(3)"
          "\n>>> c.CZ(0,1)"
          "\n>>> c.Z(1)"
          "\n>>> c.CX(1,2)"
          "\n>>> c.depth_by_type({OpType.CZ, OpType.CX})"
          "\n2"
          "\n\n:param types: the set of operation types of interest"
          "\n:return: the circuit depth with respect to operations matching an "
          "element of `types`",
          py::arg("types"))
      .def(
          "_to_graphviz_file", &Circuit::to_graphviz_file,
          "Saves a visualisation of a circuit's DAG to a \".dot\" file",
          py::arg("filename"))
      .def(
          "to_dict", [](const Circuit &c) { return json(c); },
          ":return: a JSON serializable dictionary representation of "
          "the Circuit")
      .def_static(
          "from_dict", [](const json &j) { return j.get<Circuit>(); },
          "Construct Circuit instance from JSON serializable "
          "dictionary representation of the Circuit.")
      .def(py::pickle(
          [](py::object self) {  // __getstate__
            return py::make_tuple(self.attr("to_dict")());
          },
          [](const py::tuple &t) {  // __setstate__
            const json j = t[0].cast<json>();
            return j.get<Circuit>();
          }))
      .def(
          "to_latex_file", &Circuit::to_latex_file,
          "Produces a latex file with a visualisation of the circuit "
          "using the Quantikz package.\n\n:param filename: Name of file "
          "to write output to (must end in \".tex\")",
          py::arg("filename"))
      .def(
          py::self >> py::self,
          "Creates a new Circuit, corresponding to the sequential "
          "composition of the given Circuits. Any qubits/bits with the "
          "same ids will be unified. Any ids without a match will be "
          "added in parallel.")
      .def(
          py::self * py::self,
          "Creates a new Circuit, corresponding to the parallel "
          "composition of the given Circuits. This will fail if the "
          "circuits share qubits/bits with the same ids.")
      .def(
          "dagger", &Circuit::dagger,
          "Given a pure circuit (i.e. without any measurements or "
          "conditional gates), produces a new circuit for the "
          "inverse/adjoint operation.\n\n:return: a new "
          ":py:class:`Circuit` corresponding to the inverse operation")
      .def(
          "transpose", &Circuit::transpose,
          "Given a pure circuit (i.e. without any measurements or "
          "conditional gates), produces a new circuit for the "
          "transpose operation.\n\n:return: a new "
          ":py:class:`Circuit` corresponding to the transpose operation")
      .def(
          "copy", [](Circuit &circ) { return Circuit(circ); },
          ":return: an identical copy of the circuit")
      .def(
          "symbol_substitution",
          (void(Circuit::*)(const symbol_map_t &)) &
              Circuit::symbol_substitution,
          "In-place substitution for symbolic expressions; iterates "
          "through each parameterised gate and performs the "
          "substitution. This will not affect any symbols captured "
          "within boxed operations.\n\n:param symbol_map: A map from "
          "SymPy symbols to SymPy expressions",
          py::arg("symbol_map"))
      .def(
          "symbol_substitution",
          (void(Circuit::*)(
              const std::map<Sym, double, SymEngine::RCPBasicKeyLess> &)) &
              Circuit::symbol_substitution,
          "In-place substitution for symbolic expressions; iterates "
          "through each parameterised gate and performs the "
          "substitution. This will not affect any symbols captured "
          "within boxed operations.\n\n:param symbol_map: A map from "
          "SymPy symbols to floating-point values",
          py::arg("symbol_map"))
      .def(
          "free_symbols", &Circuit::free_symbols,
          ":return: set of symbolic parameters in the circuit")
      .def(
          "is_symbolic", &Circuit::is_symbolic,
          ":return: True if the circuit "
          "contains any free symbols, False otherwise.")
      .def(
          "substitute_named",
          [](Circuit &circ, Op_ptr op, const std::string opgroup) {
            return circ.substitute_named(op, opgroup);
          },
          "Substitute all ops with the given name for the given op."
          "The replacement operations retain the same name.\n\n"
          ":param op: the replacement operation\n"
          ":param opgroup: the name of the operations group to replace\n"
          ":return: whether any replacements were made",
          py::arg("op"), py::arg("opgroup"))
      .def(
          "substitute_named",
          [](Circuit &circ, const Circuit &repl, const std::string opgroup) {
            return circ.substitute_named(repl, opgroup);
          },
          "Substitute all ops with the given name for the given circuit."
          "Named operations in the replacement circuit must not match "
          "any named operations in the circuit being modified.\n\n"
          ":param repl: the replacement circuit\n"
          ":param opgroup: the name of the operations group to replace\n"
          ":return: whether any replacements were made",
          py::arg("repl"), py::arg("opgroup"))
      .def(
          "substitute_named",
          [](Circuit &circ, const CircBox &box, const std::string opgroup) {
            return circ.substitute_named(box, opgroup);
          },
          "Substitute all ops with the given name for the given box."
          "The replacement boxes retain the same name.\n\n"
          ":param box: the replacement CircBox\n"
          ":param opgroup: the name of the operations group to replace\n"
          ":return: whether any replacements were made",
          py::arg("box"), py::arg("opgroup"))
      .def(
          "substitute_named",
          [](Circuit &circ, const Unitary1qBox &box,
             const std::string opgroup) {
            return circ.substitute_named(box, opgroup);
          },
          "Substitute all ops with the given name for the given box."
          "The replacement boxes retain the same name.\n\n"
          ":param box: the replacement Unitary1qBox\n"
          ":param opgroup: the name of the operations group to replace\n"
          ":return: whether any replacements were made",
          py::arg("box"), py::arg("opgroup"))
      .def(
          "substitute_named",
          [](Circuit &circ, const Unitary2qBox &box,
             const std::string opgroup) {
            return circ.substitute_named(box, opgroup);
          },
          "Substitute all ops with the given name for the given box."
          "The replacement boxes retain the same name.\n\n"
          ":param box: the replacement Unitary2qBox\n"
          ":param opgroup: the name of the operations group to replace\n"
          ":return: whether any replacements were made",
          py::arg("box"), py::arg("opgroup"))
      .def(
          "substitute_named",
          [](Circuit &circ, const Unitary3qBox &box,
             const std::string opgroup) {
            return circ.substitute_named(box, opgroup);
          },
          "Substitute all ops with the given name for the given box."
          "The replacement boxes retain the same name.\n\n"
          ":param box: the replacement Unitary3qBox\n"
          ":param opgroup: the name of the operations group to replace\n"
          ":return: whether any replacements were made",
          py::arg("box"), py::arg("opgroup"))
      .def(
          "substitute_named",
          [](Circuit &circ, const ExpBox &box, const std::string opgroup) {
            return circ.substitute_named(box, opgroup);
          },
          "Substitute all ops with the given name for the given box."
          "The replacement boxes retain the same name.\n\n"
          ":param box: the replacement ExpBox\n"
          ":param opgroup: the name of the operations group to replace\n"
          ":return: whether any replacements were made",
          py::arg("box"), py::arg("opgroup"))
      .def(
          "substitute_named",
          [](Circuit &circ, const PauliExpBox &box, const std::string opgroup) {
            return circ.substitute_named(box, opgroup);
          },
          "Substitute all ops with the given name for the given box."
          "The replacement boxes retain the same name.\n\n"
          ":param box: the replacement PauliExpBox\n"
          ":param opgroup: the name of the operations group to replace\n"
          ":return: whether any replacements were made",
          py::arg("box"), py::arg("opgroup"))
      .def(
          "substitute_named",
          [](Circuit &circ, const QControlBox &box, const std::string opgroup) {
            return circ.substitute_named(box, opgroup);
          },
          "Substitute all ops with the given name for the given box."
          "The replacement boxes retain the same name.\n\n"
          ":param box: the replacement QControlBox\n"
          ":param opgroup: the name of the operations group to replace\n"
          ":return: whether any replacements were made",
          py::arg("box"), py::arg("opgroup"))
      .def(
          "substitute_named",
          [](Circuit &circ, const CustomGate &box, const std::string opgroup) {
            return circ.substitute_named(box, opgroup);
          },
          "Substitute all ops with the given name for the given box."
          "The replacement boxes retain the same name.\n\n"
          ":param box: the replacement CustomGate\n"
          ":param opgroup: the name of the operations group to replace\n"
          ":return: whether any replacements were made",
          py::arg("box"), py::arg("opgroup"))

      // Methods for contextual optimization
      .def(
          "qubit_create", &Circuit::qubit_create,
          "Make a quantum input a Create operation (initialized to 0")
      .def(
          "qubit_discard", &Circuit::qubit_discard,
          "Make a quantum output a Discard operation")
      .def(
          "qubit_create_all", &Circuit::qubit_create_all,
          "Make all quantum inputs Create operations (initialized to 0)")
      .def(
          "qubit_discard_all", &Circuit::qubit_discard_all,
          "Make all quantum outputs Discard operations")
      .def(
          "qubit_is_created", &Circuit::is_created,
          "Query whether a qubit has its initial state set to zero")
      .def(
          "qubit_is_discarded", &Circuit::is_discarded,
          "Query whether a qubit has its final state discarded")
      .def("_classical_eval", &Circuit::classical_eval)

      .def(
          "valid_connectivity", &respects_connectivity_constraints,
          "Confirms whether all two qubit gates in given circuit are "
          "along some edge of the architecture."
          "\n\n:param arch: The architecture capturing the desired "
          "connectivity"
          "\n:param directed: If true, also checks that CX or ECR gates are in "
          "the same direction as the edges of the architecture"
          "\n:param allow_bridge: Accept BRIDGEs as valid, assuming the "
          "middle qubit neighbours the others"
          "\n\n:return: True or False",
          py::arg("arch"), py::arg("directed"), py::arg("allow_bridge") = false)
      .def(
          "implicit_qubit_permutation", &Circuit::implicit_qubit_permutation,
          ":return: dictionary mapping input qubit to output qubit on "
          "the same path")
      .def(
          "ops_of_type",
          [](const Circuit &circ, OpType optype) {
            VertexSet vset = circ.get_gates_of_type(optype);
            std::list<Op_ptr> ops;
            for (const auto &v : vset) {
              ops.push_back(circ.dag[v].op);
            }
            return ops;
          },
          "Get all operations in the circuit of a given type."
          "\n\nThe order is not guaranteed."
          "\n\n:param optype: operation type"
          "\n\n:return: list of :py:class:`Op`",
          py::arg("optype"))
      .def(
          "commands_of_type", &Circuit::get_commands_of_type,
          "Get all commands in a circuit of a given type."
          "\n\nThe order is consistent with the causal order of the "
          "operations in the circuit."
          "\n\n:param optype: operation type"
          "\n\n:return: list of :py:class:`Command`",
          py::arg("optype"))
      .def_property_readonly(
          "_dag_data",
          [](Circuit &circ) {
            IndexMap im = circ.index_map();

            // subset of quantum input nodes
            std::set<unsigned> q_inputs;
            for (const Vertex &v : circ.q_inputs()) {
              q_inputs.insert(im[v]);
            }

            // subset of classical input nodes
            std::set<unsigned> c_inputs;
            for (const Vertex &v : circ.c_inputs()) {
              c_inputs.insert(im[v]);
            }

            // subset of quantum output nodes
            std::set<unsigned> q_outputs;
            for (const Vertex &v : circ.q_outputs()) {
              q_outputs.insert(im[v]);
            }

            // subset of classical output nodes
            std::set<unsigned> c_outputs;
            for (const Vertex &v : circ.c_outputs()) {
              c_outputs.insert(im[v]);
            }

            // maps from input and output nodes to unit names
            std::map<unsigned, std::string> input_names, output_names;
            for (const BoundaryElement &b : circ.boundary) {
              std::string bname = b.id_.repr();
              input_names[im[b.in_]] = bname;
              output_names[im[b.out_]] = bname;
            }

            // map from node to string description
            std::map<unsigned, std::string> node_data;
            BGL_FORALL_VERTICES(v, circ.dag, DAG) {
              node_data[im[v]] = circ.get_Op_ptr_from_Vertex(v)->get_name();
            }

            // set of tuples (source node, target node, source port, target
            // port, edge type)
            std::set<
                std::tuple<unsigned, unsigned, unsigned, unsigned, unsigned>>
                edge_data;
            BGL_FORALL_EDGES(e, circ.dag, DAG) {
              Vertex v_so = circ.source(e);
              Vertex v_ta = circ.target(e);
              unsigned v_s = im[v_so];
              unsigned v_t = im[v_ta];
              // EdgeType converted to unsigned because of some weird
              // behaviour with pybind11 conversions being
              // overwritten. TODO Do this properly.
              EdgeType etype = circ.dag[e].type;
              unsigned edge_type = (etype == EdgeType::Quantum)   ? 0
                                   : (etype == EdgeType::Boolean) ? 1
                                                                  : 2;
              edge_data.insert(
                  {v_s, v_t, circ.get_source_port(e), circ.get_target_port(e),
                   edge_type});
            }

            return std::make_tuple(
                q_inputs, c_inputs, q_outputs, c_outputs, input_names,
                output_names, node_data, edge_data);
          },
          "DAG data for circuit");
}

}  // namespace tket
