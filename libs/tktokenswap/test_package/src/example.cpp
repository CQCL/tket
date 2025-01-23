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

#include <iostream>

#include "tktokenswap/SwapFunctions.hpp"

using namespace tket;

int main() {
  SwapList swap_list;
  swap_list.push_front(get_swap(0, 1));
  std::cout << swap_list.to_vector().size() << std::endl;
  return 0;
}
