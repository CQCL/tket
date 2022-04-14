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
#include <limits>
#include <set>
#include <string>
#include <utility>

#include "../GraphTheoretic/NearNeighboursData.hpp"
#include "../GraphTheoretic/NeighboursData.hpp"
#include "../Reducing/DistancesReducer.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

/** Collect together all the data and algorithm components
 * that will be computed and needed BEFORE proper searching starts,
 * for initialisation,
 * which will also be reused during the search.
 * (Initialisation-only things will just be discarded after use).
 */
struct PreSearchComponents {
  const NeighboursData& pattern_ndata;
  const NeighboursData& target_ndata;
  NearNeighboursData pattern_near_ndata;
  NearNeighboursData target_near_ndata;

  // Contains a distance filter (checker) as well as a reducer.
  DistancesReducer distances_reducer;

  PreSearchComponents(
      const NeighboursData& pattern_ndata, const NeighboursData& target_ndata);
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
