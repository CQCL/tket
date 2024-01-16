// Copyright 2019-2024 Cambridge Quantum Computing
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

#include "TestSettings.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

// Uncomment to get a running printout for the tests.
/*
TestSettings::TestSettings()
    : os(OStreamWrapper::Target::CERR),
      print_solution_times(false),
      // print_solution_times(true),
      // print_verbose_solution_data(true) {}
      print_verbose_solution_data(false) {}
*/

// No printed output for the tests.
TestSettings::TestSettings()
    : os(OStreamWrapper::Target::NONE),
      print_solution_times(false),
      print_verbose_solution_data(false) {}

const TestSettings& TestSettings::get() {
  static const TestSettings settings;
  return settings;
}

OStreamWrapper::OStreamWrapper(Target target) : m_target(target) {}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
