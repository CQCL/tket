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
#include <vector>

namespace tket {
namespace WeightedSubgraphMonomorphism {

/** The interface is a mixture of std::stack and std::vector.
 * For convenience, still allows random access like vectors;
 * but no RESIZING operations other than push/pop (like a stack).
 * Doesn't actually resize the internal vector, to avoid memory reallocation.
 * So this is useful for T objects (e.g. vectors, ...) which are cheap
 * to clear and reuse, but expensive to construct.
 */
template <class T>
class LogicalStack {
 public:
  LogicalStack();

  /** The "logical" size, of course, not the actual size of the vector.
   * @return Would be the size of the stack, IF we had implemented it as an
   * actual stack, i.e. genuinely erasing elements whenever we pop.
   */
  std::size_t size() const;

  /** Provide random access to elements, like std::vector.
   * @param index The index of an element. No checks for validity.
   * @return x[index], as if x were a vector.
   */
  const T& operator[](unsigned index) const;

  /** Provide random access to elements, like std::vector.
   * @param index The index of an element. No checks for validity.
   * @return x[index], as if x were a vector.
   */
  T& operator[](unsigned index);

  /** Whether or not it's "logically" empty.
   * @return Equivalent to size() == 0.
   */
  bool empty() const;

  /** Like std::stack, assumes without checking that it's nonempty.
   * @return The top() element of the stack (i.e., back() element of the
   * vector).
   */
  T& top();

  /** Like std::stack, assumes without checking that it's nonempty.
   * @return The top() element of the stack (i.e., back() element of the
   * vector).
   */
  const T& top() const;

  /** Not in std::stack! Notice that it's an O(1) operation.
   * Only a logical clear, doesn't actually clear. After this, size() will be
   * zero.
   */
  void clear();

  /** Not in std::stack! When we push, we sometimes want
   * to access the old top element, to fill the new top.
   * Of course, we CANNOT just store the old top(),
   * as push() may invalidate the reference!
   * @return x[size()-2], treating it as a vector x. No checks on the size.
   */
  const T& one_below_top() const;

  /** Similar to std::stack (or emplace_back() for std::vector).
   * Increases size by 1.
   * However, very important to note that the new top will be "junk" data;
   * it's the caller's responsibility to fill it.
   */
  void push();

  /** Similar to std::stack::pop() (or std::vector::pop_back()).
   * Like std::stack, assumes without checking that it's nonempty.
   * However, note that no destructor is called; the element at top()
   * before calling is NOT erased; it just sits there as "junk" data.
   */
  void pop();

 private:
  std::vector<T> m_data;
  unsigned m_logical_size;
};

// Implementations below.

template <class T>
LogicalStack<T>::LogicalStack() : m_logical_size(0) {}

template <class T>
std::size_t LogicalStack<T>::size() const {
  return m_logical_size;
}

template <class T>
const T& LogicalStack<T>::operator[](unsigned index) const {
  return m_data[index];
}

template <class T>
T& LogicalStack<T>::operator[](unsigned index) {
  return m_data[index];
}

template <class T>
bool LogicalStack<T>::empty() const {
  return m_logical_size == 0;
}

template <class T>
T& LogicalStack<T>::top() {
  // Notice that m_data.back() would be wrong, because of course
  // the "logical" size is not the same as the actual size.
  return m_data[m_logical_size - 1];
}

template <class T>
const T& LogicalStack<T>::top() const {
  return m_data[m_logical_size - 1];
}

template <class T>
void LogicalStack<T>::clear() {
  m_logical_size = 0;
}

template <class T>
const T& LogicalStack<T>::one_below_top() const {
  return m_data[m_logical_size - 2];
}

template <class T>
void LogicalStack<T>::push() {
  ++m_logical_size;
  if (m_logical_size > m_data.size()) {
    m_data.resize(m_logical_size);
  }
}

template <class T>
void LogicalStack<T>::pop() {
  --m_logical_size;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
