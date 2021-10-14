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

#include "Program.hpp"

namespace tket {

bool Program::check_valid() const {
  bool valid = true;
  valid &= get_in_edges(entry_).empty();
  valid &= get_out_edges(entry_).size() == 1;
  valid &= get_out_edges(exit_).empty();
  valid &= entry_ != exit_;
  BGL_FORALL_VERTICES(block, flow_, FlowGraph) {
    FGEdgeVec outs = get_out_edges(block);
    if (block != exit_) {
      if (flow_[block].branch_condition) {
        valid &= outs.size() == 2;
        unsigned n_true = 0;
        unsigned n_false = 0;
        for (const FGEdge &e : outs) {
          if (get_branch(e))
            ++n_true;
          else
            ++n_false;
        }
        valid &= n_true == 1;
        valid &= n_false == 1;
      } else {
        valid &= outs.size() == 1;
        valid &= get_branch(outs.front()) == false;
      }
    }
    for (const Qubit &qb : flow_[block].circ.all_qubits()) {
      unit_lookup_t::iterator found = units_.get<TagID>().find(qb);
      valid &= found != units_.get<TagID>().end();
      valid &= found->type() == UnitType::Qubit;
    }
    for (const Bit &b : flow_[block].circ.all_bits()) {
      unit_lookup_t::iterator found = units_.get<TagID>().find(b);
      valid &= found != units_.get<TagID>().end();
      valid &= found->type() == UnitType::Bit;
    }
  }
  return valid;
}

void Program::to_graphviz_file(const std::string &filename) const {
  std::ofstream dot_file(filename);
  to_graphviz(dot_file);
}

void Program::to_graphviz(std::ostream &out) const {
  out << "digraph G {\n";

  std::map<FGVert, unsigned> index_map;
  unsigned i = 0;
  BGL_FORALL_VERTICES(v, flow_, FlowGraph) {
    index_map.insert({v, i});
    out << i << " [label = \"";
    if (flow_[v].label) {
      out << "LABEL " << *flow_[v].label << "\\n";
    }
    for (const Command c : flow_[v].circ) {
      out << c.to_str() << "\\n";
    }
    if (flow_[v].branch_condition) {
      out << "BRANCH " << flow_[v].branch_condition->repr();
    }
    out << "\"];\n";
    ++i;
  }
  BGL_FORALL_EDGES(e, flow_, FlowGraph) {
    FGVert source = get_source(e);
    FGVert target = get_target(e);
    out << index_map.at(source) << " -> " << index_map.at(target)
        << " [label = \"" << get_branch(e) << "\"];\n";
  }
  out << "}";
}

}  // namespace tket
