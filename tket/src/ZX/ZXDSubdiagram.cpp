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

ZXDiagram::Subdiagram::Subdiagram() : boundary_(), verts_() {}

ZXDiagram::Subdiagram::Subdiagram(const std::vector<Wire, WireEnd>>& cut, const ZXVertSeqSet& verts) : boundary_(cut), verts_(verts) {}

void ZXDiagram::Subdiagram::check_validity(const ZXDiagram& diag) const {
  std::set<std::pair<Wire, WireEnd>> boundary_lookup;
  for (const std::pair<Wire, WireEnd>& w : boundary_) {
      auto inserted = boundary_lookup.insert(w);
      if (!inserted.second) throw ZXError("Malformed ZX Subdiagram: Wire appears multiple times in boundary");
      if (verts_.find(diag.vertex_at_end(w.first, w.second)) == verts_.end()) throw ZXError("Malformed ZX Subdiagram: Vertex adjacent to boundary is not in vertex set");
    }
  }
  for (const ZXVert& v : verts_) {
    if (is_boundary_type(diag.get_zxtype(v))) throw ZXError("Malformed ZX Subdiagram: Contains a boundary vertex");
    for (const Wire& w : diag.adj_wires(v)) {
      ZXVert n = diag.other_end(w, v);
      if (verts_.find(n) == verts_.end()) {
        if (boundary_lookup.find({w, diag.end_of(w, v)}) == boundary_lookup.end()) throw ZXError("Malformed ZX Subdiagram: subdiagram is not closed");
      }
      else {
        if (boundary_lookup.find({w, WireEnd::Source}) == boundary_lookup.end() ^ boundary_lookup.find({w, WireEnd::Target}) == boundary_lookup.end()) throw ZXError("Malformed ZX Subdiagram: wire between two interior vertices contains one boundary");
      }
    }
  }
}

ZXDiagram ZXDiagram::Subdiagram::to_diagram(const ZXDiagram& orig) {}

}  // namespace zx

}  // namespace tket
