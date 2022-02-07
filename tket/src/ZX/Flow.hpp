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

#include "Utils/BiMapHeaders.hpp"
#include "ZX/ZXDiagram.hpp"

namespace tket {

namespace zx {

/**
 * Data structure for flow in qubit MBQC.
 * Different classes of flow exist based on the types of measurements and
 * correction sets accepted, but the contents of the flow are the same. Causal <
 * XY gflow < 3Plane gflow < Pauli flow
 *
 * `c` defines the correction set for each measured vertex.
 * `d` approximates the partial order by giving the depth of the measurement
 * from the output, i.e. d(u) < d(v) => v is measured before u.
 */
class Flow {
 public:
  Flow(
      const std::map<ZXVert, ZXVertSeqSet>& c,
      const std::map<ZXVert, unsigned>& d);

  ZXVertSeqSet c(const ZXVert& v) const;
  ZXVertSeqSet odd(const ZXVert& v, const ZXDiagram& diag) const;
  unsigned d(const ZXVert& v) const;

  void verify(const ZXDiagram& diag) const;

  void focus(const ZXDiagram& diag);

  static Flow identify_causal_flow(const ZXDiagram& diag);
  static Flow identify_xy_gflow(const ZXDiagram& diag);
  static Flow identify_pauli_flow(const ZXDiagram& diag);

 private:
  std::map<ZXVert, ZXVertSeqSet> c_;
  std::map<ZXVert, unsigned> d_;

  static std::map<ZXVert, ZXVertSeqSet> gauss_solve_correctors(
      const ZXDiagram& diag, const boost::bimap<ZXVert, unsigned>& correctors,
      const boost::bimap<ZXVert, unsigned>& preserve, const ZXVertVec& to_solve,
      const boost::bimap<ZXVert, unsigned>& ys);
};

}  // namespace zx

}  // namespace tket
