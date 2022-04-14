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

#include "WeightSubgrMono/EndToEndWrappers/MainSolverParameters.hpp"

#include "WeightSubgrMono/Common/GeneralUtils.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

MainSolverParameters::MainSolverParameters(long long timeout_millisecs)
    : timeout_ms(timeout_millisecs),
      iterations_timeout(1000000),
      terminate_with_first_full_solution(false),
      for_multiple_full_solutions_the_max_number_to_obtain(0),
      max_distance_for_domain_initialisation_distance_filter(2),
      max_distance_for_distance_reduction_during_search(10) {}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
