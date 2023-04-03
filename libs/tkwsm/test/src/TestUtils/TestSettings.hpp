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

#pragma once
#include <iostream>

// To run longer or more verbose tests etc. manually,
// simply change the settings here.
namespace tket {
namespace WeightedSubgraphMonomorphism {

// An object behaving like a std::ostream,
// except that it can be switched off to produce no output.
class OStreamWrapper {
 public:
  enum class Target { COUT, CERR, NONE };

  OStreamWrapper(Target target = Target::NONE);

  // Accepts any object which std::cout or std::cerr
  // knows how to stream.
  template <class T>
  const OStreamWrapper& operator<<(const T& obj) const {
    switch (m_target) {
      case Target::COUT:
        std::cout << obj;
        break;
      case Target::CERR:
        std::cerr << obj;
        break;
      case Target::NONE:
        // fall through
      default:
        break;
    }
    return *this;
  }

 private:
  Target m_target;
};

struct TestSettings {
  OStreamWrapper os;
  OStreamWrapper os_null;

  // If false, we ONLY print number of iterations, which DOESN'T change
  // from run to run, unlike times.
  bool print_solution_times;

  bool print_verbose_solution_data;

  TestSettings();

  static const TestSettings& get();
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
