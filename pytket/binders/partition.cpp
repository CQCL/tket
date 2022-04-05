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

#include "MeasurementSetup/MeasurementReduction.hpp"
#include "binder_json.hpp"
#include "binder_utils.hpp"
#include "typecast.hpp"

namespace py = pybind11;
namespace tket {

PYBIND11_MODULE(partition, m) {
  py::enum_<PauliPartitionStrat>(
      m, "PauliPartitionStrat",
      "Enum for available strategies to partition Pauli tensors.")
      .value(
          "NonConflictingSets", PauliPartitionStrat::NonConflictingSets,
          "Build sets of Pauli tensors in which each qubit has the "
          "same Pauli or Pauli.I. Requires no additional CX gates.")
      .value(
          "CommutingSets", PauliPartitionStrat::CommutingSets,
          "Build sets of mutually commuting Pauli tensors. Requires "
          "O(n^2) CX gates to diagonalise.");

  py::enum_<GraphColourMethod>(
      m, "GraphColourMethod",
      "Enum for available methods to perform graph colouring.")
      .value(
          "Lazy", GraphColourMethod::Lazy,
          "Does not build the graph before performing the colouring; "
          "partitions while iterating through the Pauli tensors in "
          "the input order.")
      .value(
          "LargestFirst", GraphColourMethod::LargestFirst,
          "Builds the graph and then greedily colours by iterating "
          "through the vertices, with the highest degree first.")
      .value(
          "Exhaustive", GraphColourMethod::Exhaustive,
          "Builds the graph and then systematically checks all "
          "possibilities until it finds a colouring with the minimum "
          "possible number of colours. "
          "Such colourings need not be unique. "
          "Exponential time in the worst case, but often runs "
          "much faster.");

  py::class_<MeasurementSetup::MeasurementBitMap>(
      m, "MeasurementBitMap",
      "Maps Pauli tensors to Clifford circuit indices and bits required "
      "for measurement. A MeasurementBitMap belongs to a "
      "MeasurementSetup object, "
      "and dictates which bits are to be included in the measurement. As "
      "Clifford circuits may "
      "flip the parity of the corresponding Pauli tensor, the "
      "MeasurementBitMap optionally inverts "
      "the result.")
      .def(
          py::init<unsigned, std::vector<unsigned> &, bool>(),
          "Constructs a MeasurementBitMap for some Clifford circuit "
          "index and bits, with an option to invert the result."
          "\n\n:param circ_index: which measurement circuit the "
          "measurement map refers to"
          "\n:param bits: which bits are included in the measurement"
          "\n:param invert: whether to flip the parity of the result",
          py::arg("circ_index"), py::arg("bits"), py::arg("invert") = false)
      .def("__repr__", &MeasurementSetup::MeasurementBitMap::to_str)
      .def_property_readonly(
          "circ_index", &MeasurementSetup::MeasurementBitMap::get_circ_index,
          "Clifford circuit index")
      .def_property_readonly(
          "bits", &MeasurementSetup::MeasurementBitMap::get_bits,
          "Bits to measure")
      .def_property_readonly(
          "invert", &MeasurementSetup::MeasurementBitMap::get_invert,
          "Whether result is inverted or not")
      .def(
          "to_dict",
          [](const MeasurementSetup::MeasurementBitMap &map) {
            return nlohmann::json(map);
          },
          "JSON-serializable dict representation of the MeasurementBitMap."
          "\n\n:return: dict representation of the MeasurementBitMap")
      .def_static(
          "from_dict",
          [](const nlohmann::json &j) {
            return j.get<MeasurementSetup::MeasurementBitMap>();
          },
          "Construct MeasurementBitMap instance from dict representation.");

  py::class_<MeasurementSetup>(
      m, "MeasurementSetup",
      "Encapsulates an experiment in which the expectation value of an "
      "operator is to be measured via decomposition into "
      "QubitPauliStrings. Each tensor expectation value can be measured "
      "using shots. These values are then summed together with some "
      "weights to retrieve the desired operator expctation value.")
      .def(py::init<>(), "Constructs an empty MeasurementSetup object")
      .def("__repr__", &MeasurementSetup::to_str)
      .def_property_readonly(
          "measurement_circs", &MeasurementSetup::get_circs,
          "Clifford measurement circuits.")
      .def_property_readonly(
          "results", &MeasurementSetup::get_result_map,
          "Map from Pauli strings to MeasurementBitMaps")
      .def(
          "add_measurement_circuit", &MeasurementSetup::add_measurement_circuit,
          "Add a Clifford circuit that rotates into some Pauli basis",
          py::arg("circ"))
      .def(
          "add_result_for_term",
          (void(MeasurementSetup::*)(
              const QubitPauliString &,
              const MeasurementSetup::MeasurementBitMap &)) &
              MeasurementSetup::add_result_for_term,
          "Add a new Pauli string with a corresponding BitMap", py::arg("term"),
          py::arg("result"))
      .def(
          "verify", &MeasurementSetup::verify,
          "Checks that the strings to be measured correspond to the "
          "correct strings generated by the measurement circs. Checks "
          "for parity by comparing to the `invert` flag.\n\n"
          ":return: True or False")
      .def(
          "to_dict",
          [](const MeasurementSetup &setup) { return nlohmann::json(setup); },
          "JSON-serializable dict representation of the MeasurementSetup."
          "\n\n:return: dict representation of the MeasurementSetup")
      .def_static(
          "from_dict",
          [](const nlohmann::json &j) { return j.get<MeasurementSetup>(); },
          "Construct MeasurementSetup instance from dict representation.");

  m.def(
      "measurement_reduction", measurement_reduction,
      "Automatically performs graph colouring and diagonalisation to "
      "reduce measurements required for Pauli strings."
      "\n\n:param strings: A list of `QubitPauliString` objects to be "
      "partitioned."
      "\n:param strat: The `PauliPartitionStrat` to use."
      "\n:param method: The `GraphColourMethod` to use."
      "\n:param cx_config: Whenever diagonalisation is required, use "
      "this configuration of CX gates"
      "\n:return: a :py:class:`MeasurementSetup` object",
      py::arg("strings"), py::arg("strat"),
      py::arg("method") = GraphColourMethod::Lazy,
      py::arg("cx_config") = CXConfigType::Snake);

  m.def(
      "term_sequence", term_sequence,
      "Takes in a list of QubitPauliString objects and partitions them "
      "into mutually commuting sets according to some PauliPartitionStrat, "
      "then sequences in an arbitrary order."
      "\n\n:param tensors: A list of `QubitPauliString` objects to be "
      "sequenced. Assumes that each Pauli tensor is unique, and does not "
      "combine equivalent tensors."
      "\n:param strat: The `PauliPartitionStrat` to use. Defaults to "
      "`CommutingSets`."
      "\n:param method: The `GraphColourMethod` to use."
      "\n:return: a list of lists of " CLSOBJS(QubitPauliString),
      py::arg("strings"), py::arg("strat") = PauliPartitionStrat::CommutingSets,
      py::arg("method") = GraphColourMethod::Lazy);
}

}  // namespace tket
