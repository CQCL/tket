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

#include "Program.hpp"

namespace tket {

Circuit &Program::get_circuit_ref(const FGVert &vert) {
  return flow_[vert].circ;
}

const Circuit &Program::get_circuit_ref(const FGVert &vert) const {
  return flow_[vert].circ;
}

std::optional<Bit> Program::get_condition(const FGVert &vert) const {
  return flow_[vert].branch_condition;
}

std::optional<std::string> Program::get_label(const FGVert &vert) const {
  return flow_[vert].label;
}

bool Program::get_branch(const FGEdge &edge) const {
  return flow_[edge].branch;
}

FGVertVec Program::get_successors(const FGVert &vert) const {
  FGEdgeVec outs = get_out_edges(vert);
  if (outs.size() == 1) {
    return {get_target(outs.front())};
  } else if (outs.size() == 2) {
    FGVertVec children(2);
    for (const Edge &e : outs) {
      if (get_branch(e))
        children[1] = get_target(e);
      else
        children[0] = get_target(e);
    }
    return children;
  } else
    throw ProgramError("Block does not have one or two successsors");
}

FGVertVec Program::get_predecessors(const FGVert &vert) const {
  FGEdgeVec ins = get_in_edges(vert);
  FGVertVec parents;
  std::unordered_set<FGVert> lookup;
  for (const Edge &e : ins) {
    FGVert pred = get_source(e);
    if (lookup.find(pred) == lookup.end()) {
      parents.push_back(pred);
      lookup.insert(pred);
    }
  }
  return parents;
}

FGVert Program::get_branch_successor(const FGVert &vert, bool branch) const {
  FGEdgeVec outs = get_out_edges(vert);
  for (const FGEdge &e : outs) {
    if (get_branch(e) == branch) {
      return get_target(e);
    }
  }
  throw ProgramError("Could not find successor on desired branch");
}

FGEdgeVec Program::get_in_edges(const FGVert &vert) const {
  FGEdgeVec ins;
  for (auto [it, end] = boost::in_edges(vert, flow_); it != end; ++it) {
    ins.push_back(*it);
  }
  return ins;
}

FGEdgeVec Program::get_out_edges(const FGVert &vert) const {
  FGEdgeVec outs;
  for (auto [it, end] = boost::out_edges(vert, flow_); it != end; ++it) {
    outs.push_back(*it);
  }
  return outs;
}

unsigned Program::n_in_edges(const FGVert &vert) const {
  return boost::in_degree(vert, flow_);
}

unsigned Program::n_out_edges(const FGVert &vert) const {
  return boost::out_degree(vert, flow_);
}

FGVert Program::get_source(const FGEdge &edge) const {
  return boost::source(edge, flow_);
}

FGVert Program::get_target(const FGEdge &edge) const {
  return boost::target(edge, flow_);
}

unsigned Program::get_n_vertices() const { return boost::num_vertices(flow_); }

}  // namespace tket
