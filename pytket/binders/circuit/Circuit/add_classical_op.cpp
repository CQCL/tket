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

#define STR(x) #x

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <bitset>

#include "Circuit/Circuit.hpp"
#include "Gate/OpPtrFunctions.hpp"
#include "Ops/ClassicalOps.hpp"
#include "Ops/Op.hpp"
#include "UnitRegister.hpp"
#include "add_gate.hpp"
#include "binder_utils.hpp"

namespace py = pybind11;

namespace tket {

static void apply_classical_op_to_registers(
    Circuit &circ, std::shared_ptr<const ClassicalEvalOp> op,
    const std::vector<BitRegister> &registers, const py::kwargs &kwargs) {
  unsigned n_op_args = registers.size();
  const unsigned n_bits = std::min_element(
                              registers.begin(), registers.end(),
                              [](const BitRegister &i, const BitRegister &j) {
                                return i.size() < j.size();
                              })
                              ->size();
  std::vector<Bit> args(n_bits * n_op_args);
  for (unsigned i = 0; i < n_bits; i++) {
    for (unsigned j = 0; j < n_op_args; j++) {
      args[n_op_args * i + j] = registers[j][i];
    }
  }
  std::shared_ptr<MultiBitOp> mbop = std::make_shared<MultiBitOp>(op, n_bits);
  add_gate_method<Bit>(&circ, mbop, args, kwargs);
}

void init_circuit_add_classical_op(
    py::class_<Circuit, std::shared_ptr<Circuit>> &c) {
  c.def(
       "add_c_transform",
       [](Circuit &circ, const std::vector<_tket_uint_t> &values,
          const std::vector<unsigned> &args, const std::string &name,
          const py::kwargs &kwargs) {
         unsigned n_args = args.size();
         std::shared_ptr<ClassicalTransformOp> op =
             std::make_shared<ClassicalTransformOp>(n_args, values, name);
         return add_gate_method<unsigned>(&circ, op, args, kwargs);
       },
       "Appends a purely classical transformation, defined by a table of "
       "values, to "
       "the end of the circuit."
       "\n\n:param values: table of values: bit :math:`j` (in little-endian "
       "order) of the "
       "term indexed by :math:`sum_i a_i 2^i` is output :math:`j` of the "
       "transform "
       "applied to inputs :math:`(a_i)`."
       "\n:param args: bits to which the transform is applied"
       "\n:param name: operation name"
       "\n:param kwargs: additional arguments passed to `add_gate_method` ."
       " Allowed parameters are `opgroup`,  `condition` , `condition_bits`,"
       " `condition_value`"
       "\n:return: the new :py:class:`Circuit`",
       py::arg("values"), py::arg("args"),
       py::arg("name") = "ClassicalTransform")
      .def(
          "add_c_transform",
          [](Circuit &circ, const std::vector<_tket_uint_t> &values,
             const std::vector<Bit> &args, const std::string &name,
             const py::kwargs &kwargs) {
            unsigned n_args = args.size();
            std::shared_ptr<ClassicalTransformOp> op =
                std::make_shared<ClassicalTransformOp>(n_args, values, name);
            return add_gate_method<Bit>(&circ, op, args, kwargs);
          },
          "See :py:meth:`add_c_transform`.", py::arg("values"), py::arg("args"),
          py::arg("name") = "ClassicalTransform")
      .def(
          "_add_wasm",
          [](Circuit &circ, const std::string &funcname,
             const std::string &wasm_uid,
             const std::vector<unsigned> &i32list_i,
             const std::vector<unsigned> &i32list_o,
             const std::vector<unsigned> &args,
             const py::kwargs &kwargs) -> Circuit * {
            unsigned n_args = args.size();
            std::shared_ptr<WASMOp> op = std::make_shared<WASMOp>(
                n_args, i32list_i, i32list_o, funcname, wasm_uid);

            return add_gate_method<unsigned>(&circ, op, args, kwargs);
          },
          "Add a classical function call from a wasm file to the circuit. "
          "\n\n:param funcname: name of the function that is called"
          "\n:param wasm_uid: unit id to identify the wasm file"
          "\n:param i32list_i: list of the number of bits in the input "
          "variables"
          "\n:param i32list_o: list of the number of bits in the output "
          "variables"
          "\n:param args: vector of circuit bits the wasm op should be added to"
          "\n:param kwargs: additional arguments passed to `add_gate_method` ."
          " Allowed parameters are `opgroup`,  `condition` , `condition_bits`,"
          " `condition_value`"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("funcname"), py::arg("wasm_uid"), py::arg("i32list_i"),
          py::arg("i32list_o"), py::arg("args"))
      .def(
          "_add_wasm",
          [](Circuit &circ, const std::string &funcname,
             const std::string &wasm_uid,
             const std::vector<unsigned> &i32list_i,
             const std::vector<unsigned> &i32list_o,
             const std::vector<Bit> &args,
             const py::kwargs &kwargs) -> Circuit * {
            unsigned n_args = args.size();
            std::shared_ptr<WASMOp> op = std::make_shared<WASMOp>(
                n_args, i32list_i, i32list_o, funcname, wasm_uid);

            return add_gate_method<Bit>(&circ, op, args, kwargs);
          },
          "Add a classical function call from a wasm file to the circuit. "
          "\n\n:param funcname: name of the function that is called"
          "\n:param wasm_uid: unit id to identify the wasm file"
          "\n:param i32list_i: list of the number of bits in the input "
          "variables"
          "\n:param i32list_o: list of the number of bits in the output "
          "variables"
          "\n:param args: vector of circuit bits the wasm op should be added to"
          "\n:param kwargs: additional arguments passed to `add_gate_method` ."
          " Allowed parameters are `opgroup`,  `condition` , `condition_bits`,"
          " `condition_value`"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("funcname"), py::arg("wasm_uid"), py::arg("i32list_i"),
          py::arg("i32list_o"), py::arg("args"))
      .def(
          "_add_wasm",
          [](Circuit &circ, const std::string &funcname,
             const std::string &wasm_uid,
             const std::vector<BitRegister> &list_reg_in,
             const std::vector<BitRegister> &list_reg_out,
             const py::kwargs &kwargs) -> Circuit * {
            unsigned n_args = 0;

            for (auto r : list_reg_in) {
              n_args += r.size();
            }

            for (auto r : list_reg_out) {
              n_args += r.size();
            }

            std::vector<Bit> args(n_args);
            std::vector<unsigned> i32list_i(list_reg_in.size());
            std::vector<unsigned> i32list_o(list_reg_out.size());

            unsigned i = 0;
            unsigned j = 0;
            for (auto r : list_reg_in) {
              i32list_i[i] = r.size();
              for (unsigned k = 0; k < r.size(); ++k) {
                args[j] = r[k];
                ++j;
              }
              ++i;
            }

            i = 0;
            for (auto r : list_reg_out) {
              i32list_o[i] = r.size();
              for (unsigned k = 0; k < r.size(); ++k) {
                args[j] = r[k];
                ++j;
              }
              ++i;
            }

            std::shared_ptr<WASMOp> op = std::make_shared<WASMOp>(
                n_args, i32list_i, i32list_o, funcname, wasm_uid);

            return add_gate_method<Bit>(&circ, op, args, kwargs);
          },
          "Add a classical function call from a wasm file to the circuit. "
          "\n\n:param funcname: name of the function that is called"
          "\n:param wasm_uid: unit id to identify the wasm file"
          "\n:param list_reg_in: list of the classical registers in the "
          "circuit used as inputs"
          "\n:param list_reg_out: list of the classical registers in the "
          "circuit used as outputs"
          "\n:param kwargs: additional arguments passed to `add_gate_method` ."
          " Allowed parameters are `opgroup`,  `condition` , `condition_bits`,"
          " `condition_value`"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("funcname"), py::arg("wasm_uid"), py::arg("list_reg_in"),
          py::arg("list_reg_out"))
      .def(
          "add_c_setbits",
          [](Circuit &circ, const std::vector<bool> &values,
             const std::vector<unsigned> args, const py::kwargs &kwargs) {
            std::shared_ptr<SetBitsOp> op = std::make_shared<SetBitsOp>(values);
            return add_gate_method<unsigned>(&circ, op, args, kwargs);
          },
          "Appends an operation to set some bit values."
          "\n\n:param values: values to set"
          "\n:param args: bits to set"
          "\n:param kwargs: additional arguments passed to `add_gate_method` ."
          " Allowed parameters are `opgroup`,  `condition` , `condition_bits`,"
          " `condition_value`"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("values"), py::arg("args"))
      .def(
          "add_c_setbits",
          [](Circuit &circ, const std::vector<bool> &values,
             const std::vector<Bit> args, const py::kwargs &kwargs) {
            std::shared_ptr<SetBitsOp> op = std::make_shared<SetBitsOp>(values);
            return add_gate_method<Bit>(&circ, op, args, kwargs);
          },
          "See :py:meth:`add_c_setbits`.", py::arg("values"), py::arg("args"))
      .def(
          "add_c_setreg",
          [](Circuit &circ, const _tket_uint_t value, const BitRegister &reg,
             const py::kwargs &kwargs) {
            auto bs = std::bitset<_TKET_REG_WIDTH>(value);
            std::vector<Bit> args(reg.size());
            std::vector<bool> vals(reg.size());
            for (unsigned i = 0; i < reg.size(); i++) {
              args[i] = reg[i];
              vals[i] = bs[i];
            }

            std::shared_ptr<SetBitsOp> op = std::make_shared<SetBitsOp>(vals);
            return add_gate_method<Bit>(&circ, op, args, kwargs);
          },
          "Set a classical register to an unsigned integer value. The "
          "little-endian bitwise representation of the integer is truncated to "
          "the register size, up to " STR(_TKET_REG_WIDTH) " bit width.",
          py::arg("value"), py::arg("arg"))
      .def(
          "add_c_copybits",
          [](Circuit &circ, const std::vector<unsigned> &args_in,
             const std::vector<unsigned> &args_out, const py::kwargs &kwargs) {
            unsigned n_args_in = args_in.size();
            std::shared_ptr<CopyBitsOp> op =
                std::make_shared<CopyBitsOp>(n_args_in);
            std::vector<unsigned> args = args_in;
            args.insert(args.end(), args_out.begin(), args_out.end());
            return add_gate_method<unsigned>(&circ, op, args, kwargs);
          },
          "Appends a classical copy operation"
          "\n\n:param args_in: source bits"
          "\n:param args_out: destination bits"
          "\n:param kwargs: additional arguments passed to `add_gate_method` ."
          " Allowed parameters are `opgroup`,  `condition` , `condition_bits`,"
          " `condition_value`"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("args_in"), py::arg("args_out"))
      .def(
          "add_c_copybits",
          [](Circuit &circ, const std::vector<Bit> &args_in,
             const std::vector<Bit> &args_out, const py::kwargs &kwargs) {
            unsigned n_args_in = args_in.size();
            std::shared_ptr<CopyBitsOp> op =
                std::make_shared<CopyBitsOp>(n_args_in);
            std::vector<Bit> args = args_in;
            args.insert(args.end(), args_out.begin(), args_out.end());
            return add_gate_method<Bit>(&circ, op, args, kwargs);
          },
          "See :py:meth:`add_c_copybits`.", py::arg("args_in"),
          py::arg("args_out"))
      .def(
          "add_c_copyreg",
          [](Circuit &circ, const BitRegister &input_reg,
             const BitRegister &output_reg, const py::kwargs &kwargs) {
            const unsigned width =
                std::min(input_reg.size(), output_reg.size());

            std::shared_ptr<CopyBitsOp> op =
                std::make_shared<CopyBitsOp>(width);
            std::vector<Bit> args(width * 2);
            for (unsigned i = 0; i < width; i++) {
              args[i] = input_reg[i];
              args[i + width] = output_reg[i];
            }
            return add_gate_method<Bit>(&circ, op, args, kwargs);
          },
          "Copy a classical register to another. Copying is truncated to the "
          "size of the smaller of the two registers.",
          py::arg("input_reg"), py::arg("output_reg"))
      .def(
          "add_c_predicate",
          [](Circuit &circ, const std::vector<bool> &values,
             const std::vector<unsigned> &args_in, unsigned arg_out,
             const std::string &name, const py::kwargs &kwargs) {
            unsigned n_args_in = args_in.size();
            std::shared_ptr<ExplicitPredicateOp> op =
                std::make_shared<ExplicitPredicateOp>(n_args_in, values, name);
            std::vector<unsigned> args = args_in;
            args.push_back(arg_out);
            return add_gate_method<unsigned>(&circ, op, args, kwargs);
          },
          "Appends a classical predicate, defined by a truth table, to the end "
          "of the "
          "circuit."
          "\n\n:param values: table of values: the term indexed by "
          ":math:`sum_i a_i 2^i` "
          "is the value of the predicate for inputs :math:`(a_i)`."
          "\n:param args_in: input bits for the predicate"
          "\n:param arg_out: output bit, distinct from all inputs",
          "\n:param name: operation name"
          "\n:param kwargs: additional arguments passed to `add_gate_method` ."
          " Allowed parameters are `opgroup`,  `condition` , `condition_bits`,"
          " `condition_value`"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("values"), py::arg("args_in"), py::arg("arg_out"),
          py::arg("name") = "ExplicitPredicate")
      .def(
          "add_c_predicate",
          [](Circuit &circ, const std::vector<bool> &values,
             const std::vector<Bit> &args_in, Bit arg_out,
             const std::string &name, const py::kwargs &kwargs) {
            unsigned n_args_in = args_in.size();
            std::shared_ptr<ExplicitPredicateOp> op =
                std::make_shared<ExplicitPredicateOp>(n_args_in, values, name);
            std::vector<Bit> args = args_in;
            args.push_back(arg_out);
            return add_gate_method<Bit>(&circ, op, args, kwargs);
          },
          "See :py:meth:`add_c_predicate`.", py::arg("values"),
          py::arg("args_in"), py::arg("arg_out"),
          py::arg("name") = "ExplicitPredicate")
      .def(
          "add_c_modifier",
          [](Circuit &circ, const std::vector<bool> &values,
             const std::vector<unsigned> &args_in, unsigned arg_inout,
             const std::string &name, const py::kwargs &kwargs) {
            unsigned n_args_in = args_in.size();
            std::shared_ptr<ExplicitModifierOp> op =
                std::make_shared<ExplicitModifierOp>(n_args_in, values, name);
            std::vector<unsigned> args = args_in;
            args.push_back(arg_inout);
            return add_gate_method<unsigned>(&circ, op, args, kwargs);
          },
          "Appends a classical modifying operation, defined by a truth table, "
          "to the "
          "end of the circuit."
          "\n\n:param values: table of values: the term indexed by "
          ":math:`sum_i a_i 2^i` "
          "is the value of the predicate for inputs :math:`(a_i)`, where the "
          "modified "
          "bit is the last of the :math:`a_i`."
          "\n:param args_in: input bits, excluding the modified bit"
          "\n:param arg_out: modified bit",
          "\n:param name: operation name"
          "\n:param kwargs: additional arguments passed to `add_gate_method` ."
          " Allowed parameters are `opgroup`,  `condition` , `condition_bits`,"
          " `condition_value`"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("values"), py::arg("args_in"), py::arg("arg_inout"),
          py::arg("name") = "ExplicitModifier")
      .def(
          "add_c_modifier",
          [](Circuit &circ, const std::vector<bool> &values,
             const std::vector<Bit> &args_in, Bit arg_inout,
             const std::string &name, const py::kwargs &kwargs) {
            unsigned n_args_in = args_in.size();
            std::shared_ptr<ExplicitModifierOp> op =
                std::make_shared<ExplicitModifierOp>(n_args_in, values, name);
            std::vector<Bit> args = args_in;
            args.push_back(arg_inout);
            return add_gate_method<Bit>(&circ, op, args, kwargs);
          },
          "See :py:meth:`add_c_modifier`.", py::arg("values"),
          py::arg("args_in"), py::arg("arg_inout"),
          py::arg("name") = "ExplicitModifier")
      .def(
          "add_c_and",
          [](Circuit &circ, unsigned arg0_in, unsigned arg1_in,
             unsigned arg_out, const py::kwargs &kwargs) {
            Op_ptr op;
            std::vector<unsigned> args;
            if (arg0_in == arg_out) {
              op = AndWithOp();
              args = {arg1_in, arg_out};
            } else if (arg1_in == arg_out) {
              op = AndWithOp();
              args = {arg0_in, arg_out};
            } else {
              op = AndOp();
              args = {arg0_in, arg1_in, arg_out};
            }
            return add_gate_method<unsigned>(&circ, op, args, kwargs);
          },
          "Appends a binary AND operation to the end of the circuit."
          "\n\n:param arg0_in: first input bit"
          "\n:param arg1_in: second input bit"
          "\n:param arg_out: output bit"
          "\n:param kwargs: additional arguments passed to `add_gate_method` ."
          " Allowed parameters are `opgroup`,  `condition` , `condition_bits`,"
          " `condition_value`"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("arg0_in"), py::arg("arg1_in"), py::arg("arg_out"))
      .def(
          "add_c_and",
          [](Circuit &circ, Bit arg0_in, Bit arg1_in, Bit arg_out,
             const py::kwargs &kwargs) {
            Op_ptr op;
            std::vector<Bit> args;
            if (arg0_in == arg_out) {
              op = AndWithOp();
              args = {arg1_in, arg_out};
            } else if (arg1_in == arg_out) {
              op = AndWithOp();
              args = {arg0_in, arg_out};
            } else {
              op = AndOp();
              args = {arg0_in, arg1_in, arg_out};
            }
            return add_gate_method<Bit>(&circ, op, args, kwargs);
          },
          "See :py:meth:`add_c_and`.", py::arg("arg0_in"), py::arg("arg1_in"),
          py::arg("arg_out"))
      .def(
          "add_c_or",
          [](Circuit &circ, unsigned arg0_in, unsigned arg1_in,
             unsigned arg_out, const py::kwargs &kwargs) {
            Op_ptr op;
            std::vector<unsigned> args;
            if (arg0_in == arg_out) {
              op = OrWithOp();
              args = {arg1_in, arg_out};
            } else if (arg1_in == arg_out) {
              op = OrWithOp();
              args = {arg0_in, arg_out};
            } else {
              op = OrOp();
              args = {arg0_in, arg1_in, arg_out};
            }
            return add_gate_method<unsigned>(&circ, op, args, kwargs);
          },
          "Appends a binary OR operation to the end of the circuit."
          "\n\n:param arg0_in: first input bit"
          "\n:param arg1_in: second input bit"
          "\n:param arg_out: output bit"
          "\n:param kwargs: additional arguments passed to `add_gate_method` ."
          " Allowed parameters are `opgroup`,  `condition` , `condition_bits`,"
          " `condition_value`"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("arg0_in"), py::arg("arg1_in"), py::arg("arg_out"))
      .def(
          "add_c_or",
          [](Circuit &circ, Bit arg0_in, Bit arg1_in, Bit arg_out,
             const py::kwargs &kwargs) {
            Op_ptr op;
            std::vector<Bit> args;
            if (arg0_in == arg_out) {
              op = OrWithOp();
              args = {arg1_in, arg_out};
            } else if (arg1_in == arg_out) {
              op = OrWithOp();
              args = {arg0_in, arg_out};
            } else {
              op = OrOp();
              args = {arg0_in, arg1_in, arg_out};
            }
            return add_gate_method<Bit>(&circ, op, args, kwargs);
          },
          "See :py:meth:`add_c_or`.", py::arg("arg0_in"), py::arg("arg1_in"),
          py::arg("arg_out"))
      .def(
          "add_c_xor",
          [](Circuit &circ, unsigned arg0_in, unsigned arg1_in,
             unsigned arg_out, const py::kwargs &kwargs) {
            Op_ptr op;
            std::vector<unsigned> args;
            if (arg0_in == arg_out) {
              op = XorWithOp();
              args = {arg1_in, arg_out};
            } else if (arg1_in == arg_out) {
              op = XorWithOp();
              args = {arg0_in, arg_out};
            } else {
              op = XorOp();
              args = {arg0_in, arg1_in, arg_out};
            }
            return add_gate_method<unsigned>(&circ, op, args, kwargs);
          },
          "Appends a binary XOR operation to the end of the circuit."
          "\n\n:param arg0_in: first input bit"
          "\n:param arg1_in: second input bit"
          "\n:param arg_out: output bit"
          "\n:param kwargs: additional arguments passed to `add_gate_method` ."
          " Allowed parameters are `opgroup`,  `condition` , `condition_bits`,"
          " `condition_value`"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("arg0_in"), py::arg("arg1_in"), py::arg("arg_out"))
      .def(
          "add_c_xor",
          [](Circuit &circ, Bit arg0_in, Bit arg1_in, Bit arg_out,
             const py::kwargs &kwargs) {
            Op_ptr op;
            std::vector<Bit> args;
            if (arg0_in == arg_out) {
              op = XorWithOp();
              args = {arg1_in, arg_out};
            } else if (arg1_in == arg_out) {
              op = XorWithOp();
              args = {arg0_in, arg_out};
            } else {
              op = XorOp();
              args = {arg0_in, arg1_in, arg_out};
            }
            return add_gate_method<Bit>(&circ, op, args, kwargs);
          },
          "See :py:meth:`add_c_xor`.", py::arg("arg0_in"), py::arg("arg1_in"),
          py::arg("arg_out"))
      .def(
          "add_c_not",
          [](Circuit &circ, unsigned arg_in, unsigned arg_out,
             const py::kwargs &kwargs) {
            return add_gate_method<unsigned>(
                &circ, NotOp(), {arg_in, arg_out}, kwargs);
          },
          "Appends a NOT operation to the end of the circuit."
          "\n\n:param arg_in: input bit"
          "\n:param arg_out: output bit"
          "\n:param kwargs: additional arguments passed to `add_gate_method` ."
          " Allowed parameters are `opgroup`,  `condition` , `condition_bits`,"
          " `condition_value`"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("arg_in"), py::arg("arg_out"))
      .def(
          "add_c_not",
          [](Circuit &circ, Bit arg_in, Bit arg_out, const py::kwargs &kwargs) {
            return add_gate_method<Bit>(
                &circ, NotOp(), {arg_in, arg_out}, kwargs);
          },
          "See :py:meth:`add_c_not`.", py::arg("arg_in"), py::arg("arg_out"))
      .def(
          "add_c_range_predicate",
          [](Circuit &circ, _tket_uint_t a, _tket_uint_t b,
             const std::vector<unsigned> &args_in, unsigned arg_out,
             const py::kwargs &kwargs) {
            unsigned n_args_in = args_in.size();
            std::shared_ptr<RangePredicateOp> op =
                std::make_shared<RangePredicateOp>(n_args_in, a, b);
            std::vector<unsigned> args = args_in;
            args.push_back(arg_out);
            return add_gate_method<unsigned>(&circ, op, args, kwargs);
          },
          "Appends a range-predicate operation to the end of the circuit."
          "\n\n:param minval: lower bound of input in little-endian encoding"
          "\n:param maxval: upper bound of input in little-endian encoding"
          "\n:param args_in: input bits"
          "\n:param arg_out: output bit (distinct from input bits)"
          "\n:param kwargs: additional arguments passed to `add_gate_method` ."
          " Allowed parameters are `opgroup`,  `condition` , `condition_bits`,"
          " `condition_value`"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("minval"), py::arg("maxval"), py::arg("args_in"),
          py::arg("arg_out"))
      .def(
          "add_c_range_predicate",
          [](Circuit &circ, _tket_uint_t a, _tket_uint_t b,
             const std::vector<Bit> &args_in, Bit arg_out,
             const py::kwargs &kwargs) {
            unsigned n_args_in = args_in.size();
            std::shared_ptr<RangePredicateOp> op =
                std::make_shared<RangePredicateOp>(n_args_in, a, b);
            std::vector<Bit> args = args_in;
            args.push_back(arg_out);
            return add_gate_method<Bit>(&circ, op, args, kwargs);
          },
          "Appends a range-predicate operation to the end of the circuit."
          "\n\n:param minval: lower bound of input in little-endian encoding"
          "\n:param maxval: upper bound of input in little-endian encoding"
          "\n:param args_in: input bits"
          "\n:param arg_out: output bit (distinct from input bits)"
          "\n:param kwargs: additional arguments passed to `add_gate_method` ."
          " Allowed parameters are `opgroup`,  `condition` , `condition_bits`,"
          " `condition_value`"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("minval"), py::arg("maxval"), py::arg("args_in"),
          py::arg("arg_out"))
      .def(
          "add_c_and_to_registers",
          [](Circuit &circ, const BitRegister &reg0_in,
             const BitRegister &reg1_in, const BitRegister &reg_out,
             const py::kwargs &kwargs) {
            if (reg0_in == reg_out) {
              apply_classical_op_to_registers(
                  circ, AndWithOp(), {reg1_in, reg_out}, kwargs);
            } else if (reg1_in == reg_out) {
              apply_classical_op_to_registers(
                  circ, AndWithOp(), {reg0_in, reg_out}, kwargs);
            } else {
              apply_classical_op_to_registers(
                  circ, AndOp(), {reg0_in, reg1_in, reg_out}, kwargs);
            }
            return circ;
          },
          "Applies bitwise AND to linear registers."
          "\n\nThe operation is applied to the bits with indices 0, 1, 2, ... "
          "in "
          "each register, up to the size of the smallest register."
          "\n\n:param reg0_in: first input register"
          "\n:param reg1_in: second input register"
          "\n:param reg_out: output register"
          "\n:param kwargs: additional arguments passed to `add_gate_method` ."
          " Allowed parameters are `opgroup`,  `condition` , `condition_bits`,"
          " `condition_value`"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("reg0_in"), py::arg("reg1_in"), py::arg("reg_out"))
      .def(
          "add_c_or_to_registers",
          [](Circuit &circ, const BitRegister &reg0_in,
             const BitRegister &reg1_in, const BitRegister &reg_out,
             const py::kwargs &kwargs) {
            if (reg0_in == reg_out) {
              apply_classical_op_to_registers(
                  circ, OrWithOp(), {reg1_in, reg_out}, kwargs);
            } else if (reg1_in == reg_out) {
              apply_classical_op_to_registers(
                  circ, OrWithOp(), {reg0_in, reg_out}, kwargs);
            } else {
              apply_classical_op_to_registers(
                  circ, OrOp(), {reg0_in, reg1_in, reg_out}, kwargs);
            }
            return circ;
          },
          "Applies bitwise OR to linear registers."
          "\n\nThe operation is applied to the bits with indices 0, 1, 2, ... "
          "in "
          "each register, up to the size of the smallest register."
          "\n\n:param reg0_in: first input register"
          "\n:param reg1_in: second input register"
          "\n:param reg_out: output register"
          "\n:param kwargs: additional arguments passed to `add_gate_method` ."
          " Allowed parameters are `opgroup`,  `condition` , `condition_bits`,"
          " `condition_value`"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("reg0_in"), py::arg("reg1_in"), py::arg("reg_out"))
      .def(
          "add_c_xor_to_registers",
          [](Circuit &circ, const BitRegister &reg0_in,
             const BitRegister &reg1_in, const BitRegister &reg_out,
             const py::kwargs &kwargs) {
            if (reg0_in == reg_out) {
              apply_classical_op_to_registers(
                  circ, XorWithOp(), {reg1_in, reg_out}, kwargs);
            } else if (reg1_in == reg_out) {
              apply_classical_op_to_registers(
                  circ, XorWithOp(), {reg0_in, reg_out}, kwargs);
            } else {
              apply_classical_op_to_registers(
                  circ, XorOp(), {reg0_in, reg1_in, reg_out}, kwargs);
            }
            return circ;
          },
          "Applies bitwise XOR to linear registers."
          "\n\nThe operation is applied to the bits with indices 0, 1, 2, ... "
          "in "
          "each register, up to the size of the smallest register."
          "\n\n:param reg0_in: first input register"
          "\n:param reg1_in: second input register"
          "\n:param reg_out: output register"
          "\n:param kwargs: additional arguments passed to `add_gate_method` ."
          " Allowed parameters are `opgroup`,  `condition` , `condition_bits`,"
          " `condition_value`"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("reg0_in"), py::arg("reg1_in"), py::arg("reg_out"))
      .def(
          "add_c_not_to_registers",
          [](Circuit &circ, const BitRegister &reg_in,
             const BitRegister &reg_out, const py::kwargs &kwargs) {
            apply_classical_op_to_registers(
                circ, NotOp(), {reg_in, reg_out}, kwargs);
            return circ;
          },
          "Applies bitwise NOT to linear registers."
          "\n\nThe operation is applied to the bits with indices 0, 1, 2, ... "
          "in "
          "each register, up to the size of the smallest register."
          "\n\n:param reg_in: input register"
          "\n:param reg_out: name of output register"
          "\n:param kwargs: additional arguments passed to `add_gate_method` ."
          " Allowed parameters are `opgroup`,  `condition` , `condition_bits`,"
          " `condition_value`"
          "\n:return: the new :py:class:`Circuit`",
          py::arg("reg_in"), py::arg("reg_out"));
}

}  // namespace tket
