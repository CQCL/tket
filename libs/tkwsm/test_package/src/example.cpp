// Copyright 2019-2023 Cambridge Quantum Computing
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

#include <iostream>

#include "tkwsm/Searching/NodesRawData.hpp"

using namespace tket;
using namespace WeightedSubgraphMonomorphism;

int main() {
  DomainInitialiser::InitialDomains initial_domains(4);
  for (auto& domain_bitset : initial_domains) {
    domain_bitset.resize(100);
  }
  initial_domains[0].set(0);
  initial_domains[0].set(1);
  initial_domains[3].set(2);
  initial_domains[1].set(17);
  initial_domains[2].set(77);
  initial_domains[2].set(88);

  NodesRawData nodes_raw_data(initial_domains, 100);
  std::cout << nodes_raw_data.nodes_data[0].str() << std::endl;

  return 0;
}
