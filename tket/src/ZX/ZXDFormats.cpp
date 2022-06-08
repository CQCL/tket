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

bool ZXDiagram::is_graphlike() const {
  BGL_FORALL_EDGES(w, *graph, ZXGraph) {
    if (is_boundary_type(get_zxtype(source(w))) ||
        is_boundary_type(get_zxtype(target(w)))) {
      if (get_wire_type(w) != ZXWireType::Basic) return false;
    } else {
      if (get_wire_type(w) != ZXWireType::H) return false;
    }
  }
  BGL_FORALL_VERTICES(v, *graph, ZXGraph) {
    ZXType type = get_zxtype(v);
    if (type != ZXType::ZSpider && !is_boundary_type(type)) return false;
  }
  return true;
}

bool ZXDiagram::is_MBQC() const {
  BGL_FORALL_EDGES(w, *graph, ZXGraph) {
    if (get_qtype(w) != QuantumType::Quantum) return false;
    if (is_boundary_type(get_zxtype(source(w))) ||
        is_boundary_type(get_zxtype(target(w)))) {
      if (get_wire_type(w) != ZXWireType::Basic) return false;
    } else {
      if (get_wire_type(w) != ZXWireType::H) return false;
    }
  }
  BGL_FORALL_VERTICES(v, *graph, ZXGraph) {
    ZXType type = get_zxtype(v);
    if (!is_MBQC_type(type) && type != ZXType::Input && type != ZXType::Output)
      return false;
    if (get_qtype(v) != QuantumType::Quantum) return false;
  }
  return true;
}

}  // namespace zx

}  // namespace tket
