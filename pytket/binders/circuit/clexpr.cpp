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

#include "tket/Ops/ClExpr.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>
#include <sstream>
#include <stdexcept>
#include <variant>
#include <vector>

#include "UnitRegister.hpp"
#include "binder_json.hpp"
#include "deleted_hash.hpp"
#include "py_operators.hpp"

namespace py = pybind11;

namespace tket {

static std::string qasm_bit_repr(
    const ClExprTerm &term, const std::map<int, Bit> &input_bits) {
  if (const uint64_t *n = std::get_if<uint64_t>(&term)) {
    switch (*n) {
      case 0:
        return "0";
      case 1:
        return "1";
      default:
        throw std::logic_error("Invalid integer in bit operation");
    }
  } else {
    ClExprVar var = std::get<ClExprVar>(term);
    if (const ClBitVar *bvar = std::get_if<ClBitVar>(&var)) {
      const Bit b = input_bits.at(bvar->index);
      return b.repr();
    } else {
      throw std::logic_error("Expected bit variable, found register variable");
    }
  }
}

static std::string qasm_reg_repr(
    const ClExprTerm &term, const std::map<int, BitRegister> &input_regs) {
  if (const uint64_t *n = std::get_if<uint64_t>(&term)) {
    std::stringstream ss;
    ss << *n;
    return ss.str();
  } else {
    ClExprVar var = std::get<ClExprVar>(term);
    if (const ClRegVar *rvar = std::get_if<ClRegVar>(&var)) {
      const BitRegister r = input_regs.at(rvar->index);
      return r.name();
    } else {
      throw std::logic_error("Expected register variable, found bit variable");
    }
  }
}

enum class ArgValueType { Bit, Reg };

static std::string qasm_expr_repr(
    const ClExpr &expr, const std::map<int, Bit> &input_bits,
    const std::map<int, BitRegister> &input_regs);

static std::string qasm_arg_repr(
    const ClExprArg &arg, const std::map<int, Bit> &input_bits,
    const std::map<int, BitRegister> &input_regs, const ArgValueType typ) {
  if (const ClExpr *expr = std::get_if<ClExpr>(&arg)) {
    return qasm_expr_repr(*expr, input_bits, input_regs);
  } else {
    if (typ == ArgValueType::Bit) {
      return qasm_bit_repr(std::get<ClExprTerm>(arg), input_bits);
    } else {
      return qasm_reg_repr(std::get<ClExprTerm>(arg), input_regs);
    }
  }
}

static std::string qasm_expr_repr(
    const ClExpr &expr, const std::map<int, Bit> &input_bits,
    const std::map<int, BitRegister> &input_regs) {
  const ClOp op = expr.get_op();
  const std::vector<ClExprArg> args = expr.get_args();
  const unsigned n_args = args.size();
  std::stringstream ss;
  ss << "(";
  switch (op) {
    case ClOp::INVALID:
      throw std::logic_error("Invalid expression.");

    case ClOp::BitAnd:
      if (n_args == 0) {
        ss << "1";
      } else {
        for (unsigned i = 0; i < n_args; i++) {
          ss << qasm_arg_repr(
              args[i], input_bits, input_regs, ArgValueType::Bit);
          if (i + 1 < n_args) {
            ss << " & ";
          }
        }
      }
      break;

    case ClOp::BitOr:
      if (n_args == 0) {
        ss << "0";
      } else {
        for (unsigned i = 0; i < n_args; i++) {
          ss << qasm_arg_repr(
              args[i], input_bits, input_regs, ArgValueType::Bit);
          if (i + 1 < n_args) {
            ss << " | ";
          }
        }
      }
      break;

    case ClOp::BitXor:
      if (n_args == 0) {
        ss << "0";
      } else {
        for (unsigned i = 0; i < n_args; i++) {
          ss << qasm_arg_repr(
              args[i], input_bits, input_regs, ArgValueType::Bit);
          if (i + 1 < n_args) {
            ss << " ^ ";
          }
        }
      }
      break;

    case ClOp::BitEq:
      if (n_args != 2) {
        throw std::logic_error("BitEq with != 2 arguments");
      }
      ss << qasm_arg_repr(args[0], input_bits, input_regs, ArgValueType::Bit);
      ss << " == ";
      ss << qasm_arg_repr(args[1], input_bits, input_regs, ArgValueType::Bit);
      break;

    case ClOp::BitNeq:
      if (n_args != 2) {
        throw std::logic_error("BitNeq with != 2 arguments");
      }
      ss << qasm_arg_repr(args[0], input_bits, input_regs, ArgValueType::Bit)
         << " != "
         << qasm_arg_repr(args[1], input_bits, input_regs, ArgValueType::Bit);
      break;

    case ClOp::BitNot:
      if (n_args != 1) {
        throw std::logic_error("BitNot with != 1 argument");
      }
      ss << "~"
         << qasm_arg_repr(args[0], input_bits, input_regs, ArgValueType::Bit);
      break;

    case ClOp::BitZero:
      if (n_args != 0) {
        throw std::logic_error("BitZero with != 0 arguments");
      }
      ss << "0";
      break;

    case ClOp::BitOne:
      if (n_args != 0) {
        throw std::logic_error("BitOne with != 0 arguments");
      }
      ss << "1";
      break;

    case ClOp::RegAnd:
      if (n_args == 0) {
        ss << "-1";
      } else {
        for (unsigned i = 0; i < n_args; i++) {
          ss << qasm_arg_repr(
              args[i], input_bits, input_regs, ArgValueType::Reg);
          if (i + 1 < n_args) {
            ss << " & ";
          }
        }
      }
      break;

    case ClOp::RegOr:
      if (n_args == 0) {
        ss << "0";
      } else {
        for (unsigned i = 0; i < n_args; i++) {
          ss << qasm_arg_repr(
              args[i], input_bits, input_regs, ArgValueType::Reg);
          if (i + 1 < n_args) {
            ss << " | ";
          }
        }
      }
      break;

    case ClOp::RegXor:
      if (n_args == 0) {
        ss << "0";
      } else {
        for (unsigned i = 0; i < n_args; i++) {
          ss << qasm_arg_repr(
              args[i], input_bits, input_regs, ArgValueType::Reg);
          if (i + 1 < n_args) {
            ss << " ^ ";
          }
        }
      }
      break;

    case ClOp::RegEq:
      if (n_args != 2) {
        throw std::logic_error("RegEq with != 2 arguments");
      }
      ss << qasm_arg_repr(args[0], input_bits, input_regs, ArgValueType::Reg);
      ss << " == ";
      ss << qasm_arg_repr(args[1], input_bits, input_regs, ArgValueType::Reg);
      break;

    case ClOp::RegNeq:
      if (n_args != 2) {
        throw std::logic_error("RegNeq with != 2 arguments");
      }
      ss << qasm_arg_repr(args[0], input_bits, input_regs, ArgValueType::Reg)
         << " != "
         << qasm_arg_repr(args[1], input_bits, input_regs, ArgValueType::Reg);
      break;

    case ClOp::RegNot:
      if (n_args != 1) {
        throw std::logic_error("RegNot with != 1 argument");
      }
      ss << "~"
         << qasm_arg_repr(args[0], input_bits, input_regs, ArgValueType::Reg);
      break;

    case ClOp::RegZero:
      if (n_args != 0) {
        throw std::logic_error("RegZero with != 0 arguments");
      }
      ss << "0";
      break;

    case ClOp::RegOne:
      if (n_args != 0) {
        throw std::logic_error("RegOne with != 0 arguments");
      }
      ss << "-1";
      break;

    case ClOp::RegLt:
      if (n_args != 2) {
        throw std::logic_error("RegLt with != 2 arguments");
      }
      ss << qasm_arg_repr(args[0], input_bits, input_regs, ArgValueType::Reg);
      ss << " < ";
      ss << qasm_arg_repr(args[1], input_bits, input_regs, ArgValueType::Reg);
      break;

    case ClOp::RegGt:
      if (n_args != 2) {
        throw std::logic_error("RegGt with != 2 arguments");
      }
      ss << qasm_arg_repr(args[0], input_bits, input_regs, ArgValueType::Reg);
      ss << " > ";
      ss << qasm_arg_repr(args[1], input_bits, input_regs, ArgValueType::Reg);
      break;

    case ClOp::RegLeq:
      if (n_args != 2) {
        throw std::logic_error("RegLeq with != 2 arguments");
      }
      ss << qasm_arg_repr(args[0], input_bits, input_regs, ArgValueType::Reg);
      ss << " <= ";
      ss << qasm_arg_repr(args[1], input_bits, input_regs, ArgValueType::Reg);
      break;

    case ClOp::RegGeq:
      if (n_args != 2) {
        throw std::logic_error("RegGeq with != 2 arguments");
      }
      ss << qasm_arg_repr(args[0], input_bits, input_regs, ArgValueType::Reg);
      ss << " >= ";
      ss << qasm_arg_repr(args[1], input_bits, input_regs, ArgValueType::Reg);
      break;

    case ClOp::RegAdd:
      if (n_args == 0) {
        ss << "0";
      } else {
        for (unsigned i = 0; i < n_args; i++) {
          ss << qasm_arg_repr(
              args[i], input_bits, input_regs, ArgValueType::Reg);
          if (i + 1 < n_args) {
            ss << " + ";
          }
        }
      }
      break;

    case ClOp::RegSub:
      if (n_args != 2) {
        throw std::logic_error("RegSub with != 2 arguments");
      }
      ss << qasm_arg_repr(args[0], input_bits, input_regs, ArgValueType::Reg);
      ss << " - ";
      ss << qasm_arg_repr(args[1], input_bits, input_regs, ArgValueType::Reg);
      break;

    case ClOp::RegMul:
      if (n_args == 0) {
        ss << "1";
      } else {
        for (unsigned i = 0; i < n_args; i++) {
          ss << qasm_arg_repr(
              args[i], input_bits, input_regs, ArgValueType::Reg);
          if (i + 1 < n_args) {
            ss << " * ";
          }
        }
      }
      break;

    case ClOp::RegDiv:
      if (n_args != 2) {
        throw std::logic_error("RegDiv with != 2 arguments");
      }
      ss << qasm_arg_repr(args[0], input_bits, input_regs, ArgValueType::Reg);
      ss << " / ";
      ss << qasm_arg_repr(args[1], input_bits, input_regs, ArgValueType::Reg);
      break;

    case ClOp::RegPow:
      if (n_args != 2) {
        throw std::logic_error("RegPow with != 2 arguments");
      }
      ss << qasm_arg_repr(args[0], input_bits, input_regs, ArgValueType::Reg);
      ss << " ** ";
      ss << qasm_arg_repr(args[1], input_bits, input_regs, ArgValueType::Reg);
      break;

    case ClOp::RegLsh:
      if (n_args != 2) {
        throw std::logic_error("RegLsh with != 2 arguments");
      }
      ss << qasm_arg_repr(args[0], input_bits, input_regs, ArgValueType::Reg);
      ss << " << ";
      ss << qasm_arg_repr(args[1], input_bits, input_regs, ArgValueType::Reg);
      break;

    case ClOp::RegRsh:
      if (n_args != 2) {
        throw std::logic_error("RegRsh with != 2 arguments");
      }
      ss << qasm_arg_repr(args[0], input_bits, input_regs, ArgValueType::Reg);
      ss << " >> ";
      ss << qasm_arg_repr(args[1], input_bits, input_regs, ArgValueType::Reg);
      break;

    case ClOp::RegNeg:
      if (n_args != 1) {
        throw std::logic_error("RegNeg with != 1 argument");
      }
      ss << "-"
         << qasm_arg_repr(args[0], input_bits, input_regs, ArgValueType::Reg);
      break;
  }
  ss << ")";
  return ss.str();
}

void init_clexpr(py::module &m) {
  py::enum_<ClOp>(m, "ClOp", "A classical operation", py::arithmetic())
      .value("INVALID", ClOp::INVALID, "Invalid")
      .value("BitAnd", ClOp::BitAnd, "Bitwise AND")
      .value("BitOr", ClOp::BitOr, "Bitwise OR")
      .value("BitXor", ClOp::BitXor, "Bitwise XOR")
      .value("BitEq", ClOp::BitEq, "Bitwise equality")
      .value("BitNeq", ClOp::BitNeq, "Bitwise inequality")
      .value("BitNot", ClOp::BitNot, "Bitwise NOT")
      .value("BitZero", ClOp::BitZero, "Constant zero bit")
      .value("BitOne", ClOp::BitOne, "Constant one bit")
      .value("RegAnd", ClOp::RegAnd, "Registerwise AND")
      .value("RegOr", ClOp::RegOr, "Registerwise OR")
      .value("RegXor", ClOp::RegXor, "Registerwise XOR")
      .value("RegEq", ClOp::RegEq, "Registerwise equality")
      .value("RegNeq", ClOp::RegNeq, "Registerwise inequality")
      .value("RegNot", ClOp::RegNot, "Registerwise NOT")
      .value("RegZero", ClOp::RegZero, "Constant all-zeros register")
      .value("RegOne", ClOp::RegOne, "Constant all-ones register")
      .value("RegLt", ClOp::RegLt, "Integer less-than comparison")
      .value("RegGt", ClOp::RegGt, "Integer greater-than comparison")
      .value("RegLeq", ClOp::RegLeq, "Integer less-than-or-equal comparison")
      .value("RegGeq", ClOp::RegGeq, "Integer greater-than-or-equal comparison")
      .value("RegAdd", ClOp::RegAdd, "Integer addition")
      .value("RegSub", ClOp::RegSub, "Integer subtraction")
      .value("RegMul", ClOp::RegMul, "Integer multiplication")
      .value("RegDiv", ClOp::RegDiv, "Integer division")
      .value("RegPow", ClOp::RegPow, "Integer exponentiation")
      .value("RegLsh", ClOp::RegLsh, "Left shift")
      .value("RegRsh", ClOp::RegRsh, "Right shift")
      .value("RegNeg", ClOp::RegNeg, "Integer negation");

  py::class_<ClBitVar, std::shared_ptr<ClBitVar>>(
      m, "ClBitVar", "A bit variable within an expression")
      .def(
          py::init<unsigned>(),
          "Construct from an integer identifier.\n\n"
          ":param i: integer identifier for the variable",
          py::arg("i"))
      .def("__eq__", &py_equals<ClBitVar>)
      .def(
          "__str__",
          [](const ClBitVar &var) {
            std::stringstream ss;
            ss << var;
            return ss.str();
          })
      .def(
          "__repr__",
          [](const ClBitVar &var) {
            std::stringstream ss;
            ss << "ClBitVar(" << var.index << ")";
            return ss.str();
          })
      .def("__hash__", [](const ClBitVar &var) { return var.index; })
      .def_property_readonly(
          "index", [](const ClBitVar &var) { return var.index; },
          "integer identifier for the variable");

  py::class_<ClRegVar, std::shared_ptr<ClRegVar>>(
      m, "ClRegVar", "A register variable within an expression")
      .def(
          py::init<unsigned>(),
          "Construct from an integer identifier.\n\n"
          ":param i: integer identifier for the variable",
          py::arg("i"))
      .def("__eq__", &py_equals<ClRegVar>)
      .def(
          "__str__",
          [](const ClRegVar &var) {
            std::stringstream ss;
            ss << var;
            return ss.str();
          })
      .def(
          "__repr__",
          [](const ClRegVar &var) {
            std::stringstream ss;
            ss << "ClRegVar(" << var.index << ")";
            return ss.str();
          })
      .def("__hash__", [](const ClRegVar &var) { return var.index; })
      .def_property_readonly(
          "index", [](const ClRegVar &var) { return var.index; },
          "integer identifier for the variable");

  py::class_<ClExpr, std::shared_ptr<ClExpr>>(
      m, "ClExpr", "A classical expression")
      .def(
          py::init<ClOp, std::vector<ClExprArg>>(),
          "Construct from an operation type and a list of arguments.\n\n"
          ":param op: the operation type\n"
          ":param args: list of arguments to the expression (which may be "
          "integers, :py:class:`ClBitVar` variables, :py:class:`ClRegVar` "
          "variables, or other :py:class:`ClExpr`)",
          py::arg("op"), py::arg("args"))
      .def("__eq__", &py_equals<ClExpr>)
      .def(
          "__str__",
          [](const ClExpr &expr) {
            std::stringstream ss;
            ss << expr;
            return ss.str();
          })
      .def("__hash__", &deletedHash<ClExpr>, deletedHashDocstring)
      .def_property_readonly("op", &ClExpr::get_op, "main operation")
      .def_property_readonly("args", &ClExpr::get_args, "arguments")
      .def(
          "as_qasm",
          [](const ClExpr &expr, const std::map<int, Bit> input_bits,
             const std::map<int, BitRegister> input_regs) -> std::string {
            return qasm_expr_repr(expr, input_bits, input_regs);
          },
          "QASM-style string representation given corresponding bits and "
          "registers",
          py::arg("input_bits"), py::arg("input_regs"));

  py::class_<WiredClExpr, std::shared_ptr<WiredClExpr>>(
      m, "WiredClExpr",
      "An operation defined by a classical expression over a sequence of bits")
      .def(
          py::init<
              ClExpr, std::map<unsigned, unsigned>,
              std::map<unsigned, std::vector<unsigned>>,
              std::vector<unsigned>>(),
          "Construct from an expression with bit and register positions.\n\n"
          ":param expr: an abstract classical expression\n"
          ":param bit_posn: a map whose keys are the indices of the "
          ":py:class:`ClBitVar` occurring in the expression, and whose values "
          "are the positions of the corresponding bits in the arguments of the "
          "operation\n"
          ":param reg_posn: a map whose keys are the indices of the "
          ":py:class:`ClRegVar` occurring in the expression, and whose values "
          "are the sequences of positions of the corresponding bits in the "
          "arguments of the operation\n"
          ":param output_posn: a list giving the positions of the output bits "
          "in the arguments of the operation",
          py::arg("expr"), py::arg("bit_posn") = std::map<unsigned, unsigned>(),
          py::arg("reg_posn") = std::map<unsigned, std::vector<unsigned>>(),
          py::arg("output_posn") = std::vector<unsigned>())
      .def("__eq__", &py_equals<WiredClExpr>)
      .def(
          "__str__",
          [](const WiredClExpr &expr) {
            std::stringstream ss;
            ss << expr;
            return ss.str();
          })
      .def("__hash__", &deletedHash<WiredClExpr>, deletedHashDocstring)
      .def_property_readonly("expr", &WiredClExpr::get_expr, "expression")
      .def_property_readonly(
          "bit_posn", &WiredClExpr::get_bit_posn, "bit positions")
      .def_property_readonly(
          "reg_posn", &WiredClExpr::get_reg_posn, "register positions")
      .def_property_readonly(
          "output_posn", &WiredClExpr::get_output_posn, "output positions")
      .def(
          "to_dict",
          [](const WiredClExpr &wexpr) {
            return py::object(nlohmann::json(wexpr)).cast<py::dict>();
          },
          ":return: JSON-serializable dict representation")
      .def_static(
          "from_dict",
          [](const py::dict &wexpr_dict) {
            return nlohmann::json(wexpr_dict).get<WiredClExpr>();
          },
          "Construct from JSON-serializable dict representation");

  py::class_<ClExprOp, std::shared_ptr<ClExprOp>, Op>(
      m, "ClExprOp", "An operation defined by a classical expression")
      .def(
          py::init<WiredClExpr>(),
          "Construct from a wired classical expression")
      .def_property_readonly("type", &ClExprOp::get_type, "operation type")
      .def_property_readonly(
          "expr", &ClExprOp::get_wired_expr, "wired expression");
}

}  // namespace tket
