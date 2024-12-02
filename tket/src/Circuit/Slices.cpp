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

#include "tket/Circuit/Slices.hpp"

#include "tket/Circuit/Circuit.hpp"

namespace tket {

/*SliceIterator related methods*/
SliceIterator::SliceIterator(const Circuit& circ) : cut_(), circ_(&circ) {
  cut_.init();
  // add qubits to u_frontier
  for (const Qubit& q : circ.all_qubits()) {
    Vertex in = circ.get_in(q);
    cut_.slice->push_back(in);
    cut_.u_frontier->insert({q, circ.get_nth_out_edge(in, 0)});
  }

  // add bits to u_frontier and b_frontier
  for (const Bit& b : circ.all_bits()) {
    Vertex in = circ.get_in(b);
    cut_.slice->push_back(in);
    cut_.b_frontier->insert({b, circ.get_nth_b_out_bundle(in, 0)});
    cut_.u_frontier->insert({b, circ.get_nth_out_edge(in, 0)});
  }

  // add WasmState to u_frontier
  for (unsigned i = 0; i < circ._number_of_wasm_wires; ++i) {
    Vertex in = circ.get_in(circ.wasmwire[i]);
    cut_.slice->push_back(in);
    cut_.u_frontier->insert({circ.wasmwire[i], circ.get_nth_out_edge(in, 0)});
  }

  prev_b_frontier_ = cut_.b_frontier;
  cut_ = circ.next_cut(cut_.u_frontier, cut_.b_frontier);

  // Add all vertices that have no Quantum or Classical edges (e.g. Phase) and
  // no Boolean inputs:
  BGL_FORALL_VERTICES(v, circ.dag, DAG) {
    if (circ.n_in_edges(v) == 0 &&
        circ.n_out_edges_of_type(v, EdgeType::Quantum) == 0 &&
        circ.n_out_edges_of_type(v, EdgeType::Classical) == 0 &&
        circ.n_out_edges_of_type(v, EdgeType::WASM) == 0) {
      cut_.slice->push_back(v);
    }
  }
}

SliceIterator::SliceIterator(
    const Circuit& circ, const std::function<bool(Op_ptr)>& skip_func)
    : cut_(), circ_(&circ) {
  cut_.init();
  // add qubits to u_frontier
  for (const Qubit& q : circ.all_qubits()) {
    Vertex in = circ.get_in(q);
    cut_.u_frontier->insert({q, circ.get_nth_out_edge(in, 0)});
  }

  // add bits to u_frontier and b_frontier
  for (const Bit& b : circ.all_bits()) {
    Vertex in = circ.get_in(b);
    cut_.b_frontier->insert({b, circ.get_nth_b_out_bundle(in, 0)});
    cut_.u_frontier->insert({b, circ.get_nth_out_edge(in, 0)});
  }

  // add WasmState to u_frontier
  for (unsigned i = 0; i < circ._number_of_wasm_wires; ++i) {
    Vertex in = circ.get_in(circ.wasmwire[i]);
    cut_.slice->push_back(in);
    cut_.u_frontier->insert({circ.wasmwire[i], circ.get_nth_out_edge(in, 0)});
  }

  prev_b_frontier_ = cut_.b_frontier;
  cut_ = circ.next_cut(cut_.u_frontier, cut_.b_frontier, skip_func);
}

SliceIterator Circuit::slice_begin() const { return SliceIterator(*this); }

SliceIterator Circuit::slice_end() { return nullsit; }
const SliceIterator Circuit::nullsit = SliceIterator();

SliceIterator::Sliceholder SliceIterator::operator++(int) {
  Sliceholder ret(*cut_.slice);
  ++*this;
  return ret;
}

SliceIterator& SliceIterator::operator++() {
  if (this->finished()) {
    *this = circ_->slice_end();
    return *this;
  }
  prev_b_frontier_ = cut_.b_frontier;
  cut_ = circ_->next_cut(cut_.u_frontier, cut_.b_frontier);
  return *this;
}

bool SliceIterator::finished() const {
  for (const std::pair<UnitID, Edge>& pair :
       this->cut_.u_frontier->get<TagKey>()) {
    if (!circ_->detect_final_Op(circ_->target(pair.second))) return false;
  }
  for (const std::pair<Bit, EdgeVec>& pair :
       this->cut_.b_frontier->get<TagKey>()) {
    if (!pair.second.empty()) return false;
  }
  return true;
}

}  // namespace tket
