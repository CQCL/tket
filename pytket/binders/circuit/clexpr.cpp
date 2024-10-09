// Copyright 2019-2024 Cambridge Quantum Computing
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
#include <vector>

#include "py_operators.hpp"

namespace py = pybind11;

namespace tket {

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
          py::init<unsigned>(), "Construct from an integer identifier",
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
            ss << "ClBitVar(" << var.i << ")";
            return ss.str();
          })
      .def_property_readonly(
          "i", [](const ClBitVar &var) { return var.i; },
          ":return: integer identifier for the variable");

  py::class_<ClRegVar, std::shared_ptr<ClRegVar>>(
      m, "ClRegVar", "A register variable within an expression")
      .def(
          py::init<unsigned>(), "Construct from an integer identifier",
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
            ss << "ClRegVar(" << var.i << ")";
            return ss.str();
          })
      .def_property_readonly(
          "i", [](const ClRegVar &var) { return var.i; },
          ":return: integer identifier for the variable");

  py::class_<ClExpr, std::shared_ptr<ClExpr>>(
      m, "ClExpr", "A classical expression")
      .def(
          py::init<ClOp, std::vector<ClExprArg>>(),
          "Construct from an operation and a list of arguments", py::arg("op"),
          py::arg("args"))
      .def("__eq__", &py_equals<ClExpr>)
      .def(
          "__str__",
          [](const ClExpr &expr) {
            std::stringstream ss;
            ss << expr;
            return ss.str();
          })
      .def_property_readonly("op", &ClExpr::get_op, ":return: main operation")
      .def_property_readonly("args", &ClExpr::get_args, ":return: arguments");

  py::class_<WiredClExpr, std::shared_ptr<WiredClExpr>>(
      m, "WiredClExpr",
      "A classical expression defined over a sequence of bits")
      .def(
          py::init<
              ClExpr, std::map<unsigned, unsigned>,
              std::map<unsigned, std::vector<unsigned>>,
              std::vector<unsigned>>(),
          "Construct from an expression with bit and register positions",
          py::arg("expr"), py::arg("bit_posn") = std::map<unsigned, unsigned>(),
          py::arg("reg_posn") = std::map<unsigned, std::vector<unsigned>>(),
          py::arg("output_posn"))
      .def("__eq__", &py_equals<WiredClExpr>)
      .def(
          "__str__",
          [](const WiredClExpr &expr) {
            std::stringstream ss;
            ss << expr;
            return ss.str();
          })
      .def_property_readonly(
          "expr", &WiredClExpr::get_expr, ":return: expression")
      .def_property_readonly(
          "bit_posn", &WiredClExpr::get_bit_posn, ":return: bit positions")
      .def_property_readonly(
          "reg_posn", &WiredClExpr::get_reg_posn, ":return: register positions")
      .def_property_readonly(
          "output_posn", &WiredClExpr::get_output_posn,
          ":return: output positions");

  py::class_<ClExprOp, std::shared_ptr<ClExprOp>>(
      m, "ClExprOp", "An operation defined by a classical expression")
      .def(
          py::init<WiredClExpr>(),
          "Construct from a wired classical expression")
      .def_property_readonly(
          "expr", &ClExprOp::get_wired_expr, ":return: wired expression");
}

}  // namespace tket