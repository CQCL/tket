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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Circuit/ClassicalExpBox.hpp"
#include "Circuit/Conditional.hpp"
#include "Ops/ClassicalOps.hpp"
#include "Ops/OpJsonFactory.hpp"
#include "UnitRegister.hpp"
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
      .def(
          py::init<const Op_ptr&, unsigned, unsigned>(),
          "Construct from operation, bit width and (little-endian) value",
          py::arg("op"), py::arg("width"), py::arg("value"))
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
      m, "ClassicalOp", "Classical operation.")
      .def_property_readonly(
          "n_inputs", &ClassicalOp::get_n_i, "Number of pure inputs.")
      .def_property_readonly(
          "n_input_outputs", &ClassicalOp::get_n_io,
          "Number of pure input/output arguments.")
      .def_property_readonly(
          "n_outputs", &ClassicalOp::get_n_o, "Number of pure outputs.");
  py::class_<ClassicalEvalOp, std::shared_ptr<ClassicalEvalOp>, ClassicalOp>(
      m, "ClassicalEvalOp", "Evaluatable classical operation.");
  py::class_<SetBitsOp, std::shared_ptr<SetBitsOp>, ClassicalEvalOp>(
      m, "SetBitsOp",
      "An operation to set the values of Bits to some constants.")
      .def(
          py::init<const std::vector<bool>&>(),
          "Construct from a table of values.", py::arg("values"))
      .def_property_readonly(
          "values", &SetBitsOp::get_values, "The values to set bits to.");
  py::class_<CopyBitsOp, std::shared_ptr<CopyBitsOp>, ClassicalEvalOp>(
      m, "CopyBitsOp",
      "An operation to copy the values of Bits to other Bits.");
  py::class_<MultiBitOp, std::shared_ptr<MultiBitOp>, ClassicalEvalOp>(
      m, "MultiBitOp",
      "An operation to apply a classical op multiple times in parallel.")
      .def(
          py::init<std::shared_ptr<const ClassicalEvalOp>, unsigned>(),
          "Construct from a basic operation and a multiplier.", py::arg("op"),
          py::arg("multiplier"))
      .def_property_readonly(
          "basic_op", &MultiBitOp::get_op, "Underlying bitwise op.");
  py::class_<
      RangePredicateOp, std::shared_ptr<RangePredicateOp>, ClassicalEvalOp>(
      m, "RangePredicateOp",
      "A predicate defined by a range of values in binary encoding.")
      .def(
          py::init<unsigned, _tket_uint_t, _tket_uint_t>(),
          "Construct from a bit width, an upper bound and a lower bound.",
          py::arg("width"), py::arg("upper"), py::arg("lower"))
      .def_property_readonly(
          "upper", &RangePredicateOp::upper, "Inclusive upper bound.")
      .def_property_readonly(
          "lower", &RangePredicateOp::lower, "Inclusive lower bound.");
  py::class_<
      ClassicalExpBox<py::object>, std::shared_ptr<ClassicalExpBox<py::object>>,
      Op>(
      m, "ClassicalExpBox", "A box for holding classical expressions on Bits.")
      .def(
          py::init<unsigned, unsigned, unsigned, py::object>(),
          "Construct from signature (number of input, input/output, and output "
          "bits) and expression.",
          py::arg("n_i"), py::arg("n_io"), py::arg("n_o"), py::arg("exp"))
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
  py::class_<WASMOp, std::shared_ptr<WASMOp>, ClassicalOp>(
      m, "WASMOp",
      "An op holding an external classical call, defined by the external "
      "module id, the name of the function and the arguments. External calls "
      "can only act on entire registers (which will be interpreted as "
      "fixed-width integers).")
      .def(
          py::init<
              unsigned, std::vector<unsigned>, std::vector<unsigned>,
              const std::string&, const std::string&>(),
          "Construct from number of bits, bitwidths of inputs and outputs, "
          "function name and module id.",
          py::arg("num_bits"), py::arg("n_inputs"), py::arg("n_outputs"),
          py::arg("func_name"), py::arg("wasm_uid"))
      .def_property_readonly(
          "wasm_uid", &WASMOp::get_wasm_uid, "Wasm module id.")
      .def_property_readonly(
          "func_name", &WASMOp::get_func_name, "Name of function.")
      .def_property_readonly(
          "num_bits", &WASMOp::get_n, "Number of bits interacted with.")
      .def_property_readonly(
          "n_i32", &WASMOp::get_n_i32, "Number of integers acted on.")
      .def_property_readonly(
          "input_widths", &WASMOp::get_ni_vec, "Widths of input integers.")
      .def_property_readonly(
          "output_widths", &WASMOp::get_no_vec, "Widths of output integers.");
}
}  // namespace tket
