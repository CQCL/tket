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

#include "tket/Architecture/Architecture.hpp"

#include <pybind11/pybind11.h>

#include "binder_json.hpp"
#include "deleted_hash.hpp"
#include "py_operators.hpp"
#include "typecast.hpp"

namespace py = pybind11;
using json = nlohmann::json;

namespace tket {

PYBIND11_MODULE(architecture, m) {
  py::module::import("pytket._tket.unit_id");
  py::class_<Architecture, std::shared_ptr<Architecture>>(
      m, "Architecture",
      "Class describing the connectivity of qubits on a general device.")
      .def(py::init<>(), "Produces an empty architecture")
      .def(
          py::init([](const py::tket_custom::SequenceVec<
                       std::pair<unsigned, unsigned>> &connections) {
            return Architecture(connections);
          }),
          "The constructor for an architecture with connectivity "
          "between qubits.\n\n:param connections: A list of pairs "
          "representing qubit indices that can perform two-qubit "
          "operations",
          py::arg("connections"))
      .def(
          py::init<
              const py::tket_custom::SequenceVec<std::pair<Node, Node>> &>(),
          "The constructor for an architecture with connectivity "
          "between qubits.\n\n:param connections: A list of pairs "
          "representing Nodes that can perform two-qubit operations",
          py::arg("connections"))
      .def(
          "__repr__",
          [](const Architecture &arc) {
            return "<tket::Architecture, nodes=" +
                   std::to_string(arc.n_nodes()) + ">";
          })
      .def(
          "get_distance", &Architecture::get_distance,
          "given two nodes in Architecture, "
          "returns distance between them",
          py::arg("node_0"), py::arg("node_1"))
      .def(
          "valid_operation",
          [](const Architecture &arch,
             const py::tket_custom::SequenceVec<Node> &ids) {
            return arch.valid_operation(ids);
          },
          "Returns true if the given operation acting on the given ",
          "nodes can be executed on the Architecture connectivity graph."
          "\n\n:param uids: list of UnitIDs validity is being checked for",
          py::arg("uids"))
      .def(
          "get_adjacent_nodes", &Architecture::get_neighbour_nodes,
          "given a node, returns adjacent nodes in Architecture.",
          py::arg("node"))
      .def_property_readonly(
          "nodes", &Architecture::get_all_nodes_vec,
          "Returns all nodes of architecture as Node objects.")
      .def_property_readonly(
          "coupling", &Architecture::get_all_edges_vec,
          "Returns the coupling map of the Architecture as "
          "UnitIDs. ")
      .def(
          "to_dict",
          [](const Architecture &arch) {
            return py::object(json(arch)).cast<py::dict>();
          },
          "Return a JSON serializable dict representation of "
          "the Architecture."
          "\n\n:return: dict containing nodes and links.")
      .def_static(
          "from_dict",
          [](const py::dict &architecture_dict) {
            return json(architecture_dict).get<Architecture>();
          },
          "Construct Architecture instance from JSON serializable "
          "dict representation of the Architecture.")
      // as far as Python is concerned, Architectures are immutable
      .def(
          "__deepcopy__", [](const Architecture &arc,
                             const py::dict & = py::dict()) { return arc; })
      .def("__eq__", &py_equals<Architecture>)
      .def("__hash__", &deletedHash<Architecture>, deletedHashDocstring);
  py::class_<SquareGrid, std::shared_ptr<SquareGrid>, Architecture>(
      m, "SquareGrid",
      "Inherited Architecture class for qubits arranged in a square lattice of "
      "given number of rows and columns. Qubits are arranged with qubits "
      "values increasing first along rows then along columns i.e. for a "
      "3 x 3 grid:\n\n 0 1 2\n\n 3 4 5\n\n 6 7 8")
      .def(
          py::init([](const unsigned &x, const unsigned &y,
                      const std::string &label) {
            return SquareGrid(x, y, 1, label);
          }),
          "The constructor for a Square Grid architecture with some "
          "undirected connectivity between qubits.\n\n:param n_rows: "
          "The number of rows in the grid\n:param n_columns: The number "
          "of columns in the grid\n:param label: Name for Node in SquareGrid "
          "Architecture",
          py::arg("n_rows"), py::arg("n_columns"),
          py::arg("label") = "gridNode")
      .def(
          py::init<
              const unsigned, const unsigned, const unsigned,
              const std::string>(),
          "The constructor for  a Square Grid architecture with some "
          "undirected connectivity between qubits.\n\n:param n_rows: "
          "The number of rows in the grid\n:param n_columns: The number "
          "of columns in the grid\n:param n_layers: The number of "
          "layers of grids\n:param label: Name for Node in SquareGrid "
          "Architecture",
          py::arg("n_rows"), py::arg("n_columns"), py::arg("n_layers") = 1,
          py::arg("label") = "gridNode")
      .def(
          "squind_to_qind",
          [](const SquareGrid &self, const unsigned row, const unsigned col) {
            return self.squind_to_qind(row, col);
          },
          "Converts a (row,column) index for a square grid to a "
          "single "
          "qubit index\n\n:param row: The given row index\n:param "
          "column: The given column index\n:return: the "
          "corresponding "
          "global qubit index",
          py::arg("row"), py::arg("column"))
      .def(
          "qind_to_squind", &SquareGrid::qind_to_squind,
          "Converts a single qubit index to a (row,column) index for a "
          "square grid.\n\n:param index: The global qubit "
          "index\n:return: the corresponding grid index as a pair "
          "(row,column)",
          py::arg("index"))
      // as far as Python is concerned, Architectures are immutable
      .def(
          "__deepcopy__",
          [](const SquareGrid &arc, py::dict = py::dict()) { return arc; })
      .def("__repr__", [](const SquareGrid &arc) {
        return "<tket::SquareGrid, rows=" + std::to_string(arc.get_rows()) +
               ", columns=" + std::to_string(arc.get_columns()) +
               ", layers=" + std::to_string(arc.get_layers()) + ">";
      });
  py::class_<RingArch, std::shared_ptr<RingArch>, Architecture>(
      m, "RingArch",
      "Inherited Architecture class for number of qubits arranged in a ring.")
      .def(
          py::init<const unsigned, const std::string>(),
          "The constructor for a RingArchitecture with some undirected "
          "connectivity between qubits.\n\n:param number of qubits:"
          "\n:param label: Name for Node in RingArch Architecture",
          py::arg("nodes"), py::arg("label") = "ringNode")
      .def(
          "__deepcopy__",
          [](const RingArch &arc, py::dict = py::dict()) { return arc; })
      .def("__repr__", [](const RingArch &arc) {
        return "<tket::RingArch, nodes=" + std::to_string(arc.n_nodes()) + ">";
      });
  py::class_<FullyConnected>(
      m, "FullyConnected",
      "A specialised non-Architecture object emulating an architecture with "
      "all qubits connected. "
      "Not compatible with Routing or Placement methods.")
      .def(
          py::init<unsigned, const std::string>(),
          "Construct a fully-connected architecture."
          "\n\n:param n: number of qubits"
          "\n:param label: Name for Node in "
          "FullyConnected Architecture",
          py::arg("n"), py::arg("label") = "fcNode")
      .def(
          "__deepcopy__", [](const FullyConnected &arc,
                             const py::dict & = py::dict()) { return arc; })
      .def(
          "__repr__",
          [](const FullyConnected &arc) {
            return "<tket::FullyConnected, nodes=" +
                   std::to_string(arc.n_nodes()) + ">";
          })
      .def("__eq__", &py_equals<FullyConnected>)
      .def("__hash__", &deletedHash<FullyConnected>, deletedHashDocstring)
      .def_property_readonly(
          "nodes", &FullyConnected::get_all_nodes_vec,
          "All nodes of the architecture as :py:class:`Node` objects.")
      .def(
          "to_dict",
          [](const FullyConnected &arch) {
            return py::object(json(arch)).cast<py::dict>();
          },
          "JSON-serializable dict representation of the architecture."
          "\n\n"
          ":return: dict containing nodes")
      .def_static(
          "from_dict",
          [](const py::dict &fully_connected_dict) {
            return json(fully_connected_dict).get<FullyConnected>();
          },
          "Construct FullyConnected instance from dict representation.");
}
}  // namespace tket
