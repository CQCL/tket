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
#include <forward_list>

namespace tket {
namespace WeightedSubgraphMonomorphism {

template <class T>
class SimpleStorage {
 public:
  typedef typename std::forward_list<T>::iterator Iter;

  /** The purpose is, the iterator (and associated pointer) is easy to copy,
   * but also the T object it points to will remain valid
   * even as other elements are added.
   */
  Iter get_new_iter();

 private:
  std::forward_list<T> m_data;
};

// Implementations

template <class T>
typename std::forward_list<T>::iterator SimpleStorage<T>::get_new_iter() {
  m_data.emplace_front();
  return m_data.begin();
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
