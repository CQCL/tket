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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Circuit/ClassicalExpBox.hpp"
#include "Ops/ClassicalOps.hpp"
#include "Ops/Conditional.hpp"
#include "Ops/OpJsonFactory.hpp"
#include "Utils/Json.hpp"
#include "binder_json.hpp"

namespace py = pybind11;
using json = nlohmann::json;

namespace tket {

template <>
json ClassicalExpBox<py::object>::to_json(const Op_ptr& op) {
  const auto& box = static_cast<const ClassicalExpBox<py::object>&>(*op);
  json j = core_box_json(box);
  j["n_i"] = box.get_n_i();
  j["n_io"] = box.get_n_io();
  j["n_o"] = box.get_n_o();
  j["exp"] = box.get_exp().attr("to_dict")();
  return j;
}

template <>
Op_ptr ClassicalExpBox<py::object>::from_json(const json& j) {
  py::module logic_exp = py::module::import("pytket.circuit.logic_exp");
  ClassicalExpBox<py::object> box = ClassicalExpBox<py::object>(
      j.at("n_i").get<unsigned>(), j.at("n_io").get<unsigned>(),
      j.at("n_o").get<unsigned>(),
      logic_exp.attr("LogicExp")
          .attr("from_dict")(j.at("exp").get<py::dict>()));
  return set_box_id(
      box,
      boost::lexical_cast<boost::uuids::uuid>(j.at("id").get<std::string>()));
}

bool registered_cexpbox = OpJsonFactory::register_method(
    OpType::ClassicalExpBox, ClassicalExpBox<py::object>::from_json,
    ClassicalExpBox<py::object>::to_json);

void init_classical(py::module& m) {
  py::class_<Conditional, std::shared_ptr<Conditional>, Op>(
      m, "Conditional",
      "A wrapper for an operation to be applied conditionally on the "
      "value of some classical bits (following the nature of conditional "
      "operations in the OpenQASM specification).")
      .def_property_readonly(
          "op", &Conditional::get_op,
          "The operation to be applied conditionally")
      .def_property_readonly(
          "width", &Conditional::get_width,
          "The number of bits in the condition register")
      .def_property_readonly(
          "value", &Conditional::get_value,
          "The little-endian value the classical register must read "
          "in order to apply the operation (e.g. value 2 (10b) means "
          "bits[0] must be 0 and bits[1] must be 1)");
  py::class_<ClassicalOp, std::shared_ptr<ClassicalOp>, Op>(
      m, "ClassicalOp",
      "An operation to set the values of Bits to some constants.")
      .def_property_readonly(
          "n_inputs", &ClassicalOp::get_n_i, "Number of pure inputs.")
      .def_property_readonly(
          "n_input_outputs", &ClassicalOp::get_n_io,
          "Number of pure input/output arguments.")
      .def_property_readonly(
          "n_outputs", &ClassicalOp::get_n_o, "Number of pure outputs.");
  py::class_<SetBitsOp, std::shared_ptr<SetBitsOp>, ClassicalOp>(
      m, "SetBitsOp",
      "An operation to set the values of Bits to some constants.")
      .def_property_readonly(
          "values", &SetBitsOp::get_values, "The values to set bits to.");
  py::class_<MultiBitOp, std::shared_ptr<MultiBitOp>, ClassicalOp>(
      m, "MultiBitOp",
      "An operation to set the values of Bits to some constants.")
      .def_property_readonly(
          "basic_op", &MultiBitOp::get_op, "Underlying bitwise op.");
  py::class_<RangePredicateOp, std::shared_ptr<RangePredicateOp>, ClassicalOp>(
      m, "RangePredicateOp",
      "A predicate defined by a range of values in binary encoding.")
      .def_property_readonly(
          "upper", &RangePredicateOp::upper, "Inclusive upper bound.")
      .def_property_readonly(
          "lower", &RangePredicateOp::lower, "Inclusive lower bound.");
  py::class_<
      ClassicalExpBox<py::object>, std::shared_ptr<ClassicalExpBox<py::object>>,
      Op>(
      m, "ClassicalExpBox", "A box for holding classical expressions on Bits.")
      .def(
          "get_exp", &ClassicalExpBox<py::object>::get_exp,
          ":return: the classical expression")
      .def(
          "get_n_i", &ClassicalExpBox<py::object>::get_n_i,
          ":return: the number of pure inputs to the box.")
      .def(
          "get_n_io", &ClassicalExpBox<py::object>::get_n_io,
          ":return: the number of inputs/outputs to the box.")
      .def(
          "get_n_o", &ClassicalExpBox<py::object>::get_n_o,
          ":return: the number of pure outputs from the box.")
      .def(
          "content_equality", &ClassicalExpBox<py::object>::content_equality,
          "Check whether two ClassicalExpBox are equal in content");
}
}  // namespace tket
