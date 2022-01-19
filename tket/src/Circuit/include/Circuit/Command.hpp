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

#pragma once

#include <optional>
#include <sstream>
#include <string>

#include "DAGDefs.hpp"
#include "Ops/OpPtr.hpp"
#include "Utils/Json.hpp"
#include "Utils/TketLog.hpp"
#include "Utils/UnitID.hpp"

namespace tket {

class BadCommand : public std::logic_error {
 public:
  explicit BadCommand(const std::string &message) : std::logic_error(message) {}
};

class Command {
 public:
  Command() : op_ptr(nullptr) {}
  Command(
      const Op_ptr _gate, unit_vector_t _args,
      const std::optional<std::string> op_group = std::nullopt,
      const Vertex _vert = boost::graph_traits<DAG>::null_vertex())
      : op_ptr(_gate), args(_args), opgroup(op_group), vert(_vert) {}
  bool operator==(const Command other) const {
    return (
        *op_ptr == *other.op_ptr && args == other.args &&
        opgroup == other.opgroup);
  }

  Op_ptr get_op_ptr() const { return op_ptr; }
  std::optional<std::string> get_opgroup() const { return opgroup; }
  unit_vector_t get_args() const { return args; }
  qubit_vector_t get_qubits() const {
    qubit_vector_t qbs;
    op_signature_t sig = op_ptr->get_signature();
    for (unsigned i = 0; i < sig.size(); ++i) {
      if (sig[i] == EdgeType::Quantum) {
        qbs.push_back(Qubit(args[i]));
      }
    }
    return qbs;
  }
  bit_vector_t get_bits() const {
    bit_vector_t bs;
    op_signature_t sig = op_ptr->get_signature();
    for (unsigned i = 0; i < sig.size(); ++i) {
      if (sig[i] == EdgeType::Classical) {
        bs.push_back(Bit(args[i]));
      }
    }
    return bs;
  }
  Vertex get_vertex() const { return vert; }
  std::string to_str() const {
    std::stringstream cmd_str;
    if (opgroup) {
      cmd_str << "[" << opgroup.value() << "] ";
    }
    cmd_str << op_ptr->get_command_str(args);
    return cmd_str.str();
  }
  friend std::ostream &operator<<(std::ostream &out, const Command &c) {
    out << c.to_str();
    return out;
  };

 private:
  Op_ptr op_ptr;
  unit_vector_t args;  // indexed by port numbering
  std::optional<std::string> opgroup;
  Vertex vert;  // vertex in the DAG
};

JSON_DECL(Command)

}  // namespace tket
