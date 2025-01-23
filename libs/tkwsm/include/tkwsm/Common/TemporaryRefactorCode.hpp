// Copyright Quantinuum
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
#include <boost/dynamic_bitset.hpp>
#include <tkassert/Assert.hpp>

#include "../GraphTheoretic/GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

// As the name suggests, this whole file should be deleted when finished!
// Also, make an object of this next to any code needing to be refactored!!
struct TemporaryRefactorCode {
  static void set_bitset(
      const std::set<VertexWSM>& domain, boost::dynamic_bitset<>& bitset,
      unsigned number_of_tv) {
    bitset.resize(number_of_tv);
    bitset.reset();
    for (auto tv : domain) {
      TKET_ASSERT(tv < bitset.size());
      bitset.set(tv);
    }
  }

  static void set_domain_from_bitset(
      std::set<VertexWSM>& domain, const boost::dynamic_bitset<>& bitset) {
    domain.clear();
    for (auto tv = bitset.find_first(); tv < bitset.size();
         tv = bitset.find_next(tv)) {
      domain.insert(tv);
    }
  }
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
