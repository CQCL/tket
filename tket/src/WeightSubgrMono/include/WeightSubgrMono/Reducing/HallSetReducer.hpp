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
#include "HallSetDetector.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

class NodeWSM;

/** See the class "HallSetDetector", which tries to find sets of
 * pattern vertices {v(i)} such that
 *    |union Domain(v(i))| <= |{v(i)}|.
 * When this occurs, reduction is possible (but the detector class doesn't
 * know how to reduce).
 * This class actually handles the reductions, if any.
 */
class HallSetReducer {
 public:
  HallSetReducer();

  /** Call as soon as you reach a new node. */
  void clear();

  enum class Result { FINISHED, NOGOOD, NEW_ASSIGNMENTS };

  Result operator()(NodeWSM& node);

 private:
  bool m_awaiting_first_reduction;
  HallSetDetector m_detector;

  // Reuse, rather than reallocate.
  // The first will be the PV;
  // the second, the new Domain(PV) after reduction.
  std::vector<std::pair<VertexWSM, std::vector<VertexWSM>>> m_new_domains_data;

  struct FillResult {
    unsigned number_of_new_domains;
    bool new_assignments;
    bool nogood;
    bool changed;
  };

  FillResult fill_new_domains_data_from_hall_set(
      const PossibleAssignments& current_domains,
      const HallSetDetector::Result& detector_result);
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
