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

#include <iostream>

#include "tkwsm/Searching/NodesRawData.hpp"

using namespace tket;
using namespace WeightedSubgraphMonomorphism;

int main() {
  DomainInitialiser::InitialDomains initial_domains(4);
  initial_domains[0] = {0, 1};
  initial_domains[3] = {2};
  initial_domains[1] = {17};
  initial_domains[2] = {77, 88};

  NodesRawData nodes_raw_data(initial_domains);
  std::cout << nodes_raw_data.nodes_data[0].str() << std::endl;

  return 0;
}
