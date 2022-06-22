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

#include <cstdint>
#include <string>
#include <vector>

namespace tket {
namespace tsa_internal {

/** This contains only support data and algorithms for VectorListHybrid,
 *  a data structure combining features of std::vector and linked lists.
 *  No checks for invalidated or never valid indices.
 *  This keeps track of indices for a std::vector of data,
 *  without actually holding any other data itself.
 *  Throughout, "after", "before", "next", "previous", "front", "back"
 *  refer to the logical ordering, AS IF elements were being inserted into
 *  and erased from a std::vector, but NOT the actual order in which elements
 *  are stored in an actual implementation (like VectorListHybrid).
 *  Erased elements are not actually erased, they are reused.
 */
class VectorListHybridSkeleton {
 public:
  /** Represents actual indices for a std::vector, which SHOULD store
   *  the objects we care about (templated on the object type; but this
   *  class stores no data except indexing information).
   */
  typedef size_t Index;

  VectorListHybridSkeleton();

  /** "Null" indices will always be represented by this value.
   *  @return An index value which is guaranteed NEVER to be valid.
   */
  static Index get_invalid_index();

  /** Indices will be valid until that element is erased,
   *  or clear() is called, regardless of other insertions/erasures.
   *  A "logical" clear; does not actually clear any data,
   *  "erased" elements will be reused.
   *  But, this is time O(n) because existing internal links will be
   *  reset to default values.
   */
  void clear();

  /** Time O(1), does not erase internal link indices. Identical erase/insert
   * calls after fast_clear() calls (i.e., respecting the ordering, but
   * ignoring the internal indices) will result in the same logical list,
   * BUT the returned Index values may be different.
   */
  void fast_clear();

  /** Reverses the logical order of the elements. Time O(n). */
  void reverse();

  /** The number of elements currently stored;
   *  NOT equal to the underlying vector size!
   *  @return The number of valid elements stored.
   */
  size_t size() const;

  /** The index of the front element (or the same index as returned by
   * get_invalid_index() if currently empty).
   *  @return The index of the front element.
   */
  Index front_index() const;

  /** The index of the back element (or the same index as returned by
   * get_invalid_index() if currently empty).
   *  @return The index of the back element.
   */
  Index back_index() const;

  // All input indices MUST be currently valid,
  // but this is not checked. (Checking would need O(log N) time,
  // since we'd have to use maps and sets).

  /** The index of the next element after the given one.
   *  @param index The index of a valid element (not checked).
   *  @return The index of the next element (or the same index as returned by
   * get_invalid_index() if no next element exists).
   */
  Index next(Index index) const;

  /** The index of the previous element before the given one.
   *  @param index The index of a valid element (not checked).
   *  @return The index of the previous element (or the same index as returned
   * by get_invalid_index() if no previous element exists).
   */
  Index previous(Index index) const;

  /** "Logical" erase of the element (the position is marked for reuse).
   *  @param index The index of a valid element (not checked).
   */
  void erase(Index index);

  /** Logical erase of an interval of linked elements (a, next(a),
   * next(next(a)), ...). Equivalent to looping with erase() and next(), but
   * more efficient. The list MUST contain enough elements to erase.
   * @param index The index of a valid element to start erasing at (not
   * checked).
   * @param number_of_elements Number of elements to erase; these MUST exist
   * (the list must be big enough).
   */
  void erase_interval(Index index, size_t number_of_elements);

  /** The list must currently be empty, but not checked. */
  void insert_for_empty_list();

  /** Insert a new element after the existing one.
   *  @param index The index of a valid element (not checked).
   */
  void insert_after(Index index);

  /** Insert a new element before the existing one.
   *  @param index The index of a valid element (not checked).
   */
  void insert_before(Index index);

  /** A platform-independent string which can be copied into tests.
   *  @return A string representing the current data, useful for testing.
   */
  std::string debug_str() const;

 private:
  struct Link {
    Index previous;
    Index next;
  };

  std::vector<Link> m_links;
  size_t m_size;
  Index m_front;
  Index m_back;

  // Deleted elements will form a second linked list for reuse
  // inside the data. TRICK: forward list only,
  // no need for doubly linked lists.
  Index m_deleted_front;

  /** Resizes m_links if necessary to ensure that the new index
   *  is valid (but will reuse erased elements if possible).
   *  However, DOESN'T set the "previous" and "next" data;
   *  the caller must do that (depending on what they're doing.
   *  Thus, it's initially an "orphan" link).
   *  @return A valid index for a new Link object (but with unset fields).
   */
  Index get_new_index();
};

}  // namespace tsa_internal
}  // namespace tket
