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
#include <cstdint>
#include <vector>

namespace tket {
namespace WeightedSubgraphMonomorphism {

struct ReusableStorageId {
  std::size_t index;

  /** Necessary to allow storing IDs in a std::set. */
  bool operator<(const ReusableStorageId& other) const {
    return index < other.index;
  }
};

/** For storing objects of type T, which can be reused.
 * Objects are "released", i.e. marked as ready for reuse - like garbage
 * collection - rather than being erased, but the caller is responsible
 * for clearing such T objects when making use of them again.
 *
 * Access the objects by a std::size_t ID, in time O(1).
 * The IDs are allowed to be reused, and references are allowed to be
 * invalidated if other elements are added or "erased" (released).
 * An ID remains valid (unlike a reference) even as others are added
 * or released, until the object with that ID is released.
 *
 * So, this is useful if we frequently need objects like std::vector
 * which are much cheaper to clear and reuse than to construct afresh.
 */
template <class T>
class ReusableStorage {
 public:
  ReusableStorageId get_new_id();

  void release(ReusableStorageId id);

  const T& get_object(ReusableStorageId id) const;

  T& get_nonconst_object(ReusableStorageId id);

 private:
  std::vector<T> m_data;
  std::vector<std::size_t> m_released_indices;
};

// Implementations

template <class T>
ReusableStorageId ReusableStorage<T>::get_new_id() {
  ReusableStorageId id;
  if (m_released_indices.empty()) {
    id.index = m_data.size();
    m_data.resize(id.index + 1);
  } else {
    id.index = m_released_indices.back();
    m_released_indices.pop_back();
  }
  return id;
}

template <class T>
void ReusableStorage<T>::release(ReusableStorageId id) {
  m_released_indices.push_back(id.index);
}

template <class T>
const T& ReusableStorage<T>::get_object(ReusableStorageId id) const {
  return m_data.at(id.index);
}

template <class T>
T& ReusableStorage<T>::get_nonconst_object(ReusableStorageId id) {
  return m_data.at(id.index);
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
