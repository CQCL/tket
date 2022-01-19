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

#include "Utils/GraphHeaders.hpp"
#include "ZX/ZXDiagram.hpp"

namespace tket {

namespace zx {

WireProperties::WireProperties() {}

WireProperties::WireProperties(
    ZXWireType _type, QuantumType _qtype, std::optional<unsigned> _source_port,
    std::optional<unsigned> _target_port)
    : type(_type),
      qtype(_qtype),
      source_port(_source_port),
      target_port(_target_port) {}

ZXDiagram::ZXDiagram() : scalar(1.) { graph = std::make_unique<ZXGraph>(); }

ZXDiagram::ZXDiagram(
    unsigned in, unsigned out, unsigned classical_in, unsigned classical_out)
    : ZXDiagram() {
  for (unsigned i = 0; i < in; ++i) {
    ZXVert iv = add_vertex(ZXType::Input, QuantumType::Quantum);
    boundary.push_back(iv);
  }
  for (unsigned o = 0; o < out; ++o) {
    ZXVert ov = add_vertex(ZXType::Output, QuantumType::Quantum);
    boundary.push_back(ov);
  }
  for (unsigned i = 0; i < classical_in; ++i) {
    ZXVert iv = add_vertex(ZXType::Input, QuantumType::Classical);
    boundary.push_back(iv);
  }
  for (unsigned o = 0; o < classical_out; ++o) {
    ZXVert ov = add_vertex(ZXType::Output, QuantumType::Classical);
    boundary.push_back(ov);
  }
}

ZXDiagram::ZXDiagram(const ZXDiagram& other) : ZXDiagram() {
  this->copy_graph(other, true);
}

ZXDiagram::ZXDiagram(ZXDiagram&& other)
    : graph(std::move(other.graph)),
      boundary(std::move(other.boundary)),
      scalar(std::move(other.scalar)) {}

ZXDiagram& ZXDiagram::operator=(const ZXDiagram& other) {
  this->graph->clear();
  this->boundary.clear();
  this->scalar = 1.;

  this->copy_graph(other, true);

  return *this;
}

ZXDiagram& ZXDiagram::operator=(ZXDiagram&& other) {
  this->graph = std::move(other.graph);
  this->boundary = std::move(other.boundary);
  this->scalar = std::move(other.scalar);

  return *this;
}

std::pair<std::map<ZXVert, ZXVert>, std::map<Wire, Wire>> ZXDiagram::copy_graph(
    const ZXDiagram& other, bool merge_boundaries) {
  // The isomorphism from vertices of `other` to vertices of `this`
  std::map<ZXVert, ZXVert> iso;
  std::map<Wire, Wire> wiso;
  BGL_FORALL_VERTICES(v, *other.graph, ZXGraph) {
    ZXGen_ptr gen = other.get_vertex_ZXGen_ptr(v);
    ZXVert new_v = this->add_vertex(gen);
    iso.emplace(v, new_v);
  }

  BGL_FORALL_EDGES(w, *other.graph, ZXGraph) {
    WireProperties wp = other.get_wire_info(w);
    ZXVert source = other.source(w);
    ZXVert target = other.target(w);
    Wire new_w = this->add_wire(iso.at(source), iso.at(target), wp);
    wiso.emplace(w, new_w);
  }

  if (merge_boundaries) {
    for (const ZXVert& b : other.boundary) {
      this->boundary.push_back(iso.at(b));
    }
  }

  this->multiply_scalar(other.get_scalar());

  return {iso, wiso};
}

}  // namespace zx

}  // namespace tket
