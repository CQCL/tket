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

#include "Characterisation/FrameRandomisation.hpp"

#include <random>

#include "Ops/MetaOp.hpp"
#include "PauliGraph/ConjugatePauliFunctions.hpp"
#include "Utils/PauliStrings.hpp"

namespace tket {

class FrameRandomisationError : public std::logic_error {
 public:
  explicit FrameRandomisationError(const std::string& message)
      : std::logic_error(message) {}
};

std::string FrameRandomisation::to_string() const {
  std::string str = "<tket::FrameRandomisation, Cycle OpTypeSet: ";
  for (const OpType& ot : cycle_types_) {
    OpDesc od(ot);
    str.append(od.name() + " ");
  }
  str.append(", Frame OpTypeSet: ");
  for (const OpType& ot : frame_types_) {
    OpDesc od(ot);
    str.append(od.name() + " ");
  }
  str.append(">");
  return str;
}

// Wires Identity gates into each cycle edge. Identity gates then relabelled
// with Ops from OpTypeSet to create instances of Frame Randomisation
void add_noop_frames(std::vector<Cycle>& cycles, Circuit& circ) {
  std::map<Edge, Edge> replacement_rewiring_edges;
  for (Cycle& full_cycle : cycles) {
    std::vector<edge_pair_t> cycle = full_cycle.boundary_edges_;
    // wire "noop" vertices into each cycle edge
    // these are later relabelled with frame gates
    // also wire "barrier" vertices between "noop" vertices and cycle gates
    // i.e. something like this:
    // ---C---
    //    |
    // ---X---
    //    to
    // -noop--B--C--B--noop-
    //        A  |  A
    // -noop--R--X--R--noop-
    //        R     R
    //        I     I
    //        E     E
    //        R     R

    EdgeVec barrier_ins;
    EdgeVec barrier_outs;
    for (const std::pair<Edge, Edge>& boundary : cycle) {
      Vertex input_noop_vert = circ.add_vertex(OpType::noop);
      Vertex output_noop_vert = circ.add_vertex(OpType::noop);
      Edge in_edge = boundary.first;
      // boundary out edges can be equal to other boundary in edges
      // rewiring a vertex into these out edges replace the edge
      // first see if its been rewired and use new edge if necessary
      std::map<Edge, Edge>::iterator it =
          replacement_rewiring_edges.find(in_edge);
      if (it != replacement_rewiring_edges.end()) {
        in_edge = (*it).second;
      }

      circ.rewire(input_noop_vert, {in_edge}, {EdgeType::Quantum});
      circ.rewire(output_noop_vert, {boundary.second}, {EdgeType::Quantum});

      // Can guarantee both have one output edge only as just rewired
      Edge input_noop_out_edge =
          circ.get_out_edges_of_type(input_noop_vert, EdgeType::Quantum)[0];
      Edge output_noop_in_edge =
          circ.get_in_edges_of_type(output_noop_vert, EdgeType::Quantum)[0];

      barrier_ins.push_back(input_noop_out_edge);
      barrier_outs.push_back(output_noop_in_edge);

      replacement_rewiring_edges[boundary.second] =
          circ.get_out_edges_of_type(output_noop_vert, EdgeType::Quantum)[0];
      full_cycle.add_vertex_pair({input_noop_vert, output_noop_vert});
    }
    std::vector<EdgeType> sig(barrier_ins.size(), EdgeType::Quantum);
    Op_ptr o_ptr = std::make_shared<MetaOp>(OpType::Barrier, sig);
    Vertex input_barrier_vert = circ.add_vertex(o_ptr);
    Vertex output_barrier_vert = circ.add_vertex(o_ptr);
    circ.rewire(input_barrier_vert, barrier_ins, sig);
    circ.rewire(output_barrier_vert, barrier_outs, sig);
  }
}

std::vector<std::vector<OpTypeVector>> get_all_frame_permutations(
    const unsigned& max_frame_size, const OpTypeSet& frame_types) {
  OpTypeVector base_otv(frame_types.begin(), frame_types.end());
  std::sort(base_otv.begin(), base_otv.end());
  std::vector<OpTypeVector> base_frame;
  for (const OpType& ot : base_otv) {
    base_frame.push_back({ot});
  }
  std::sort(base_frame.begin(), base_frame.end());
  std::vector<std::vector<OpTypeVector>> out({base_frame});

  for (unsigned i = 1; i < max_frame_size; i++) {
    std::vector<OpTypeVector> new_ops;
    for (const OpTypeVector& base_ops : out[0]) {
      for (const OpTypeVector& merge_ops : out[i - 1]) {
        OpTypeVector new_base(base_ops);
        std::copy(
            merge_ops.begin(), merge_ops.end(), std::back_inserter(new_base));
        new_ops.push_back(new_base);
      }
    }
    out.push_back(new_ops);
  }
  return out;
}

std::vector<std::vector<OpTypeVector>> combine_vectors(
    const std::vector<std::vector<OpTypeVector>>& vec_0,
    const std::vector<OpTypeVector>& vec_1) {
  std::vector<std::vector<OpTypeVector>> new_vec;

  for (const std::vector<OpTypeVector>& base_frame : vec_0) {
    for (const OpTypeVector& new_frame : vec_1) {
      std::vector<OpTypeVector> copy = base_frame;
      copy.push_back(new_frame);
      new_vec.push_back(copy);
    }
  }
  return new_vec;
}

std::vector<std::vector<OpTypeVector>> get_all_permutation_combinations(
    const std::vector<unsigned>& frame_sizes,
    const std::vector<std::vector<OpTypeVector>>& frame_permutations) {
  std::vector<std::vector<OpTypeVector>> all_combinations;

  for (const OpTypeVector& frame : frame_permutations[frame_sizes[0] - 1]) {
    all_combinations.push_back({frame});
  }

  for (unsigned i = 1; i < frame_sizes.size(); i++) {
    all_combinations = combine_vectors(
        all_combinations, frame_permutations[frame_sizes[i] - 1]);
  }
  return all_combinations;
}

std::pair<OpTypeVector, std::vector<Vertex>> FrameRandomisation::get_out_frame(
    const OpTypeVector& in_frame, const Cycle& cycle) {
  OpTypeVector out_frame = in_frame;
  for (const CycleCom& cycle_op : cycle.coms_) {
    OpTypeVector current_frame;
    for (const unsigned& u : cycle_op.indices)
      current_frame.push_back(out_frame[u]);
    if (!is_initial_q_type(cycle_op.type)) {
      current_frame =
          (frame_cycle_conjugates_.at(cycle_op.type)).at(current_frame);
    }
    for (unsigned i = 0; i < current_frame.size(); i++) {
      out_frame[cycle_op.indices[i]] = current_frame[i];
    }
  }
  return {out_frame, {}};
}

std::pair<OpTypeVector, std::vector<Vertex>>
PauliFrameRandomisation::get_out_frame(
    const OpTypeVector& in_frame, const Cycle& cycle) {
  QubitPauliMap qpm;
  for (unsigned i = 0; i < in_frame.size(); i++) {
    switch (in_frame[i]) {
      case OpType::noop:
        qpm[Qubit("frame", i)] = Pauli::I;
        break;
      case OpType::X:
        qpm[Qubit("frame", i)] = Pauli::X;
        break;
      case OpType::Y:
        qpm[Qubit("frame", i)] = Pauli::Y;
        break;
      case OpType::Z:
        qpm[Qubit("frame", i)] = Pauli::Z;
        break;
      default: {
        throw FrameRandomisationError(std::string(
            "Frame OpType " + OpDesc(in_frame[i]).name() +
            " not a Pauli OpType."));
      }
    }
  }

  QubitPauliTensor qpt(qpm);

  for (const CycleCom& cycle_op : cycle.coms_) {
    switch (cycle_op.type) {
      case OpType::noop:
      case OpType::Input:
      case OpType::Create:
      case OpType::Output:
      case OpType::Discard:
        break;
      case OpType::H:
      case OpType::S:
      case OpType::X:
      case OpType::V:
      case OpType::Z:
      case OpType::Y:
        conjugate_PauliTensor(
            qpt, cycle_op.type, Qubit("frame", cycle_op.indices[0]), false);
        break;
      case OpType::Vdg:
      case OpType::Sdg:
        conjugate_PauliTensor(
            qpt, cycle_op.type, Qubit("frame", cycle_op.indices[0]), true);
        break;
      case OpType::CX:
        conjugate_PauliTensor(
            qpt, cycle_op.type, Qubit("frame", cycle_op.indices[0]),
            Qubit("frame", cycle_op.indices[1]));
        break;
      default: {
        throw FrameRandomisationError(std::string(
            "Cycle OpType " + OpDesc(cycle_op.type).name() +
            " not supported for PauliFrameRandomisation."));
      }
    }
  }

  OpTypeVector out_frame(in_frame.size());
  for (const auto& entry : qpt.string.map) {
    switch (entry.second) {
      case Pauli::I:
        out_frame[entry.first.index()[0]] = OpType::noop;
        break;
      case Pauli::X:
        out_frame[entry.first.index()[0]] = OpType::X;
        break;
      case Pauli::Y:
        out_frame[entry.first.index()[0]] = OpType::Y;
        break;
      case Pauli::Z:
        out_frame[entry.first.index()[0]] = OpType::Z;
        break;
    }
  }
  return {out_frame, {}};
}

std::pair<OpTypeVector, std::vector<Vertex>>
UniversalFrameRandomisation::get_out_frame(
    const OpTypeVector& in_frame, const Cycle& cycle) {
  QubitPauliMap qpm;
  std::vector<Vertex> to_dagger;
  for (unsigned i = 0; i < in_frame.size(); i++) {
    switch (in_frame[i]) {
      case OpType::noop:
        qpm[Qubit("frame", i)] = Pauli::I;
        break;
      case OpType::X:
        qpm[Qubit("frame", i)] = Pauli::X;
        break;
      case OpType::Y:
        qpm[Qubit("frame", i)] = Pauli::Y;
        break;
      case OpType::Z:
        qpm[Qubit("frame", i)] = Pauli::Z;
        break;
      default: {
        throw FrameRandomisationError(std::string(
            "Frame OpType " + OpDesc(in_frame[i]).name() +
            " not a Pauli OpType."));
      }
    }
  }

  QubitPauliTensor qpt(qpm);

  for (const CycleCom& cycle_op : cycle.coms_) {
    if (cycle_op.type == OpType::Rz) {
      Pauli frame_type = qpt.string.map[Qubit("frame", cycle_op.indices[0])];
      if (frame_type == Pauli::X || frame_type == Pauli::Y) {
        to_dagger.push_back(cycle_op.address);
      }
      // don't conjugate Pauli Tensor as no frame changes for Rz
    }
    if (cycle_op.type == OpType::H) {
      conjugate_PauliTensor(
          qpt, cycle_op.type, Qubit("frame", cycle_op.indices[0]), false);
    }
    if (cycle_op.type == OpType::CX) {
      conjugate_PauliTensor(
          qpt, cycle_op.type, Qubit("frame", cycle_op.indices[0]),
          Qubit("frame", cycle_op.indices[1]));
    }
  }

  OpTypeVector out_frame(in_frame.size());
  for (const auto& entry : qpt.string.map) {
    switch (entry.second) {
      case Pauli::I:
        out_frame[entry.first.index()[0]] = OpType::noop;
        break;
      case Pauli::X:
        out_frame[entry.first.index()[0]] = OpType::X;
        break;
      case Pauli::Y:
        out_frame[entry.first.index()[0]] = OpType::Y;
        break;
      case Pauli::Z:
        out_frame[entry.first.index()[0]] = OpType::Z;
        break;
    }
  }
  return {out_frame, to_dagger};
}

std::vector<Circuit> FrameRandomisation::label_frames(
    const std::vector<std::vector<OpTypeVector>>& all_frame_ops,
    const std::vector<Cycle>& cycles) {
  std::vector<Circuit> output_circuit_list;
  // make instances of all circuits using
  for (const std::vector<OpTypeVector>& cycle_frames : all_frame_ops) {
    std::vector<Vertex> dagger_vertices;
    if (cycle_frames.size() != cycles.size()) {
      throw FrameRandomisationError(
          std::string("Length of combination of Frame Permutations does not "
                      "equal number of Cycles."));
    }
    for (unsigned i = 0; i < cycle_frames.size(); i++) {
      if (cycle_frames[i].size() != cycles[i].size()) {
        throw FrameRandomisationError(
            std::string("Size of frame does not match the number "
                        "of qubits in Cycles."));
      }
      OpTypeVector in_frame = cycle_frames[i];
      std::pair<OpTypeVector, std::vector<Vertex>> frame_and_vertices =
          get_out_frame(in_frame, cycles[i]);
      dagger_vertices.insert(
          dagger_vertices.end(), frame_and_vertices.second.begin(),
          frame_and_vertices.second.end());
      assign_vertices(
          in_frame, frame_and_vertices.first, cycles[i].get_frame());
    }

    for (const Vertex& v : dagger_vertices) {
      circuit_.set_vertex_Op_ptr(
          v, circuit_.get_Op_ptr_from_Vertex(v)->dagger());
    }
    output_circuit_list.push_back(Circuit(circuit_));

    for (const Vertex& v : dagger_vertices) {
      circuit_.set_vertex_Op_ptr(
          v, circuit_.get_Op_ptr_from_Vertex(v)->dagger());
    }
  }
  return output_circuit_list;
}

std::pair<std::vector<unsigned>, unsigned> get_frame_sizes(
    const std::vector<Cycle>& cycles) {
  unsigned max_frame_size = 0;
  std::vector<unsigned> frame_sizes;
  for (const Cycle& cycle : cycles) {
    unsigned frame_size = cycle.size();
    frame_sizes.push_back(frame_size);
    if (frame_size > max_frame_size) {
      max_frame_size = frame_size;
    }
  }
  return {frame_sizes, max_frame_size};
}

std::vector<Circuit> FrameRandomisation::get_all_circuits(const Circuit& circ) {
  // Get all cycle boundaries
  // Get all sizes
  // Get all possible frames for all possible cycles
  // Apply every permutation of frames to cycles
  // Return as vector
  // get all cycles in circuit, wire noop vertices for relabelling into
  // boundaries
  this->circuit_ = circ;
  std::vector<Cycle> all_cycles = get_cycles(circuit_);
  if (all_cycles.size() == 0) {
    throw FrameRandomisationError(
        std::string("Circuit has no gates with OpType in Cycle OpTypes."));
  }
  add_noop_frames(all_cycles, circuit_);
  std::pair<std::vector<unsigned>, unsigned> frame_sizes =
      get_frame_sizes(all_cycles);

  // work out all possible permutations of given ops for all frame sizes
  std::vector<std::vector<OpTypeVector>> all_frame_perms =
      get_all_frame_permutations(frame_sizes.second, frame_types_);
  // combine all these permutations
  std::vector<std::vector<OpTypeVector>> all_permutation_combinations =
      get_all_permutation_combinations(frame_sizes.first, all_frame_perms);
  std::vector<Circuit> output_circuit_list =
      label_frames(all_permutation_combinations, all_cycles);
  return output_circuit_list;
}
OpTypeVector FrameRandomisation::sample_frame(const unsigned& size) const {
  OpTypeVector frame;
  for (unsigned i = 0; i < size; i++)
    std::sample(
        frame_types_.begin(), frame_types_.end(), std::back_inserter(frame), 1,
        std::mt19937{std::random_device{}()});
  return frame;
}

std::vector<std::vector<OpTypeVector>> FrameRandomisation::get_all_samples(
    const unsigned& samples, const std::vector<unsigned>& frame_sizes) const {
  std::vector<std::vector<OpTypeVector>> output_frames;
  for (unsigned i = 0; i < samples; i++) {
    std::vector<OpTypeVector> single_sample;
    for (const unsigned& size : frame_sizes)
      single_sample.push_back(sample_frame(size));
    output_frames.push_back(single_sample);
  }
  return output_frames;
}

void FrameRandomisation::assign_vertices(
    const OpTypeVector& in_frame, const OpTypeVector& out_frame,
    const std::vector<std::pair<Vertex, Vertex>>& frame_vertices) {
  if (in_frame.size() != out_frame.size() ||
      in_frame.size() != frame_vertices.size()) {
    throw FrameRandomisationError(
        std::string("Number of gates in sampled frame doesn't match "
                    "number of qubits in frame"));
  }
  for (unsigned i = 0; i < frame_vertices.size(); i++) {
    circuit_.set_vertex_Op_ptr(
        frame_vertices[i].first, get_op_ptr(in_frame[i]));
    circuit_.set_vertex_Op_ptr(
        frame_vertices[i].second, get_op_ptr(out_frame[i]));
  }
}

std::vector<Circuit> FrameRandomisation::sample_randomisation_circuits(
    const Circuit& circ, unsigned samples) {
  // get all cycles in circuit, wire noop vertices for relabelling into
  // boundaries
  this->circuit_ = circ;
  std::vector<Cycle> all_cycles = get_cycles(circuit_);
  if (all_cycles.size() == 0) {
    throw FrameRandomisationError(
        std::string("Circuit has no gates with OpType in Cycle OpTypes."));
  }
  add_noop_frames(all_cycles, circuit_);
  std::vector<unsigned> frame_sizes = get_frame_sizes(all_cycles).first;
  std::vector<std::vector<OpTypeVector>> all_samples =
      get_all_samples(samples, frame_sizes);
  std::vector<Circuit> output_circuit_list =
      label_frames(all_samples, all_cycles);
  return output_circuit_list;
}

std::vector<Cycle> FrameRandomisation::get_cycles(const Circuit& circ) const {
  CycleFinder cf(circ, this->cycle_types_);
  return cf.get_cycles();
}

std::vector<Circuit> PowerCycle::sample_cycles(
    const Circuit& circ, unsigned total_cycles, unsigned samples) {
  this->circuit_ = circ;
  std::vector<Circuit> out;
  std::vector<Cycle> all_cycles = get_cycles(circuit_);
  if (all_cycles.size() == 0) {
    throw FrameRandomisationError(
        std::string("Circuit has no gates with OpType in Clifford gates."));
  }
  if (all_cycles.size() > 1) {
    throw FrameRandomisationError(
        std::string("Circuit has non-Clifford gates."));
  }

  add_noop_frames(all_cycles, circuit_);
  std::vector<unsigned> frame_sizes = get_frame_sizes(all_cycles).first;
  std::vector<std::vector<OpTypeVector>> all_samples =
      get_all_samples(samples, frame_sizes);
  for (const std::vector<OpTypeVector>& sample : all_samples) {
    if (sample.size() > 1) {
      throw FrameRandomisationError(
          std::string("Frames have been sampled for more than one cycle."));
    }
    OpTypeVector in_frame = sample[0];
    OpTypeVector noop_frame;
    for (unsigned i = 0; i < in_frame.size(); i++) {
      noop_frame.push_back(OpType::noop);
    }
    std::pair<OpTypeVector, std::vector<Vertex>> frame_and_vertices =
        get_out_frame(in_frame, all_cycles[0]);
    assign_vertices(
        in_frame, frame_and_vertices.first, all_cycles[0].get_frame());
    Circuit full_circ = circuit_;
    for (unsigned i = 0; i < total_cycles - 1; i++) {
      frame_and_vertices =
          get_out_frame(frame_and_vertices.first, all_cycles[0]);
      assign_vertices(
          noop_frame, frame_and_vertices.first, all_cycles[0].get_frame());
      full_circ.append(circuit_);
    }
    out.push_back(full_circ);
  }
  return out;
}

}  // namespace tket
