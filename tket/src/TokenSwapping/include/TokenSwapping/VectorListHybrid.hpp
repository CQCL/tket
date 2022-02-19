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

#include <optional>
#include <sstream>

#include "Utils/Assert.hpp"
#include "VectorListHybridSkeleton.hpp"

namespace tket {

struct OverwriteIntervalResult {
  size_t number_of_overwritten_elements;
  tsa_internal::VectorListHybridSkeleton::Index final_overwritten_element_id;
};

/** VectorListHybrid<T> combines some functionality of std::vector<T>
 *  and std::list<T>, with the following goals:
 *
 *  Objects are stored internally inside a std::vector.
 *
 *  UNLIKE STL linked lists: erasure/insertion does NOT cause dynamic
 *      memory allocation/deallocation (except when more space
 *       is needed, in which case a vector reallocation takes place).
 *
 *  All operations are O(1), except insertions which are amortised O(1)
 *     (because a vector reallocation may be needed for more storage space).
 *
 *  Objects are not actually destroyed, they are merely marked for later reuse.
 *  Thus this class is good when objects are expensive to construct,
 *  but cheap to reuse and clear, and will be reused many times.
 *  (E.g., imagine a std::vector<std::vector<A>> being repeatedly resized;
 *  all those inner std::vector<A> are repeatedly deallocated and reallocated).
 *
 *  Objects can be accessed at any position, via an ID (like a vector index).
 *
 *  Erasure/insertion does NOT invalidate other IDs, unless that element
 *  was erased (or the whole container cleared).
 *
 *  NOTE: "previous" and "next" directions, by analogy with std::vector,
 *          correspond to the logical order the elements are regarded to have,
 *          AS IF they sat in a vector which we iterated through in forwards
 * order (which, of course, is unrelated to where they are actually stored
 * internally). Thus, "next", "forward" moves go towards the BACK; "previous",
 * "backward" moves go towards the FRONT. This should not confuse if we remember
 * std::vector itself, with begin() and rbegin() iterators.
 *
 *  TODO: there are no O(log N) operations, and no checks for invalid indices.
 *      This could be achieved by wrapping this class and storing
 *      sets/maps of erased/inserted IDs, etc. etc. Then everything would become
 *      O(log N) or amortised O(log N) instead of O(1), but we'd also have
 * complete checks.
 *
 *  TODO: this class should have its own tests. Right now it is only used
 *      in other things (SwapListOptimiser) which do have end-to-end tests,
 *      so it's quite reliable but not as reliable as it could be.
 *
 *  TODO: Once this is well tested, move it to Utils for wider use.
 */
template <class T>
class VectorListHybrid {
 public:
  /** NOTE: the ID is NOT necessarily an actual vector index;
   *      that's an implementation detail.
   */
  typedef tsa_internal::VectorListHybridSkeleton::Index ID;

  VectorListHybrid();

  /** Returns an ID which is guaranteed NEVER to be valid.
   *  @return an ID value guaranteed NEVER to be valid.
   */
  static ID get_invalid_id();

  /** Logical clear: doesn't actually delete the elements,
   *  just relabels them for reuse. Time O(N).
   *  After this, all data - even IDs - will behave AS IF
   *  it were a new object.
   */
  void clear();

  /** Logical clear: doesn't actually delete the elements,
   *  just relabels them for reuse. Time O(1).
   *  After calling this function, IDs related to
   *  inserting/erasing elements may be different from
   *  those which would be obtained by the same sequence
   *  of operations on a new object.
   */
  void fast_clear();

  /** Like std::reverse, reverses the (logical) order of the elements. (Not the
   * physical order: the internal vector of T objects is unchanged, only the
   * links are changed). Existing ids may be invalidated. Time O(n).
   */
  void reverse();

  bool empty() const;

  /** The number of valid elements stored (not, of course, the actual
   *  internal number of elements, which is larger if some are waiting
   *  to be reused).
   *  @return The number of active elements stored.
   */
  size_t size() const;

  /** Exactly like std::vector push_back. Fine if T is lightweight.
   *  Otherwise, maybe better to reuse elements.
   *  @param elem The T object to be copied and stored.
   */
  void push_back(const T& elem);

  /** Like push_back, creates a new element after the current back,
   *  but returns the ID for the new element (which of course might not
   *  really be new; it is for reuse - it may be an old T object).
   *  Of course the returned ID is the same as would be obtained
   *  from back_id().
   *  @return The ID of the newly created (or reused) element.
   */
  ID emplace_back();

  /** Erase the element at the back, but no checks for validity. */
  void pop_back();

  /** Like push_back, but instead inserts the new element before
   *  the existing front element (so that it becomes the new front).
   *  @param elem The T object to be copied and stored.
   */
  void push_front(const T& elem);

  /** Like emplace_back(), but creates the new element at the front,
   *  like push_front. However, returns the ID of the new object
   *  at the front.
   *  @return The ID of the newly created (or reused) element, at the front.
   */
  ID emplace_front();

  /** Erase the element at the front, but no checks for validity. */
  void pop_front();

  /** Creates a new element after the existing one (not checked).
   *  @param id The ID of an existing element.
   *  @return The ID of the new element, inserted immediately after
   *    (i.e., "next"; towards the BACK) of the given element.
   */
  ID insert_after(ID id);

  /** Creates a new element before the existing one (not checked).
   *  @param id The ID of an existing element.
   *  @return The ID of the new element, inserted immediately before
   *    (i.e., "previous"; towards the FRONT) of the given element.
   */
  ID insert_before(ID id);

  /** Just like std::vector back().
   *  Retrieve the element for reuse; must exist!
   *  @return A reference to the existing element at the back.
   */
  T& back();

  /** Retrieve the element for reuse; must exist!
   *  @return A reference to the existing element at the front.
   */
  T& front();

  /** Retrieve the stored element at the existing ID (not checked!)
   *  @param id The ID of an existing element.
   *  @return A reference to the element.
   */
  T& at(ID id);

  /** Retrieve the stored element at the existing ID (not checked!)
   *  @param id The ID of an existing element.
   *  @return A reference to the element.
   */
  const T& at(ID) const;

  /** Get the element ID after the given one (which MUST be valid),
   *  or a null ID if we're already at the back.
   *  @param id The ID of an existing element.
   *  @return The ID of the element after it (towards the BACK),
   *    or null if it doesn't exist.
   */
  std::optional<ID> next(ID id) const;

  /** Get the ID of the element after the given one.
   *  @param id The ID of an existing element, OR null if none exists.
   *  @return The ID of the element after it (towards the BACK),
   *    OR null if it doesn't exist, or no ID was specified.
   */
  std::optional<ID> next(std::optional<ID> id) const;

  /** Like next. Get the element ID before the given one (which MUST be valid),
   *  or a null ID if we're already at the front.
   *  @param id The ID of an existing element.
   *  @return The ID of the element before it (towards the FRONT),
   *    or null if it doesn't exist.
   */
  std::optional<ID> previous(ID id) const;

  /** The ID of the back() element, if it exists.
   *  @return The ID of the element at back(), or null if there is none.
   */
  std::optional<ID> back_id() const;

  /** The ID of the front() element, if it exists.
   *  @return The ID of the element at front(), or null if there is none.
   */
  std::optional<ID> front_id() const;

  /** Erase the element with that ID, whilst updating other links
   *  (the ID must actually exist).
   *  @param id The ID of the existing element to erase.
   */
  void erase(ID id);

  /** Starting with the given ID, erase the given number of elements.
   * Equivalent to looping with erase() and next(), but more efficient.
   * The list MUST contain enough elements to erase.
   * @param id The ID of the initial existing element to erase. Must be valid.
   * @param number_of_elements The number of elements to erase. The list MUST
   * contain enough elements to be erased.
   */
  void erase_interval(ID id, size_t number_of_elements);

  /** Starting with the given ID, and given cbegin, cend iterators to a
   * container of T objects, overwrite whatever T objects are currently stored
   * in the list with the new T objects. The list MUST be big enough to allow
   * overwriting all of them. The container of T objects MUST be nonempty.
   * @param id The ID of the initial existing T element to overwrite. Must be
   * valid.
   * @param new_elements_cbegin Const iterator to the start of a sequence of new
   * T elements.
   * @param new_elements_cend Const iterator to the cend of a sequence of new T
   * elements.
   * @return The ID of the last T element that was overwritten; MUST be valid!
   */
  template <class CIter>
  OverwriteIntervalResult overwrite_interval(
      ID id, const CIter& new_elements_cbegin, const CIter& new_elements_cend);

  /** Returns an ordinary vector of the data (in the correct order,
   *  maybe not the same as the internal storage order of course).
   *  @return A copy of the valid T objects stored, in the correct LOGICAL
   *    order, AS IF they had been inserted into a vector object throughout.
   *    (Of course, probably not the same as the actual storage order).
   */
  std::vector<T> to_vector() const;

  /** Doesn't clear the vector, but copies all elements to the end of it.
   *  @param vect A vector, which will have all the valid elements in this
   *    object pushed back to it.
   */
  void append_to_vector(std::vector<T>& vect) const;

  /** Only for debugging purposes.
   *  @return A string giving further details of the internal data.
   */
  std::string debug_str() const;

 private:
  tsa_internal::VectorListHybridSkeleton m_links_data;

  /// The actual stored elements.
  std::vector<T> m_data;

  /** Returns the ID if valid, or null if not.
   *  @param id An ID, maybe invalid.
   *  @return The ID again, if valid, or null if not.
   */
  static std::optional<ID> optional_id(ID id);

  /** Checks if m_data is big enough for the ID (which is really an index,
   *  returned by m_links_data). If not, resizes m_data if necessary,
   *  and just returns the ID unchanged.
   *  @param id An ID, valid for m_links_data, but maybe not for m_data.
   *  @return The passed in ID, but now definitely valid.
   */
  ID get_checked_new_id(ID id);

  /** The list must currently be empty (but not checked). Creates a new
   *  element, resizes m_data if necessary, and returns the ID.
   *  @return The ID of the newly created (or reused) element.
   */
  ID insert_for_empty_list();
};

template <class T>
VectorListHybrid<T>::VectorListHybrid() {}

template <class T>
typename VectorListHybrid<T>::ID VectorListHybrid<T>::get_invalid_id() {
  return tsa_internal::VectorListHybridSkeleton::get_invalid_index();
}

template <class T>
std::optional<typename VectorListHybrid<T>::ID>
VectorListHybrid<T>::optional_id(ID id) {
  if (id == tsa_internal::VectorListHybridSkeleton::get_invalid_index()) {
    return {};
  }
  return id;
}

template <class T>
void VectorListHybrid<T>::clear() {
  m_links_data.clear();
}

template <class T>
void VectorListHybrid<T>::fast_clear() {
  m_links_data.fast_clear();
}

template <class T>
void VectorListHybrid<T>::reverse() {
  m_links_data.reverse();
}
template <class T>
bool VectorListHybrid<T>::empty() const {
  return m_links_data.size() == 0;
}

template <class T>
size_t VectorListHybrid<T>::size() const {
  return m_links_data.size();
}

template <class T>
void VectorListHybrid<T>::push_back(const T& elem) {
  emplace_back();
  back() = elem;
}

template <class T>
typename VectorListHybrid<T>::ID VectorListHybrid<T>::emplace_back() {
  if (empty()) {
    insert_for_empty_list();
  } else {
    insert_after(m_links_data.back_index());
  }
  return m_links_data.back_index();
}

template <class T>
void VectorListHybrid<T>::pop_back() {
  erase(m_links_data.back_index());
}

template <class T>
void VectorListHybrid<T>::push_front(const T& elem) {
  emplace_front();
  front() = elem;
}

template <class T>
typename VectorListHybrid<T>::ID VectorListHybrid<T>::emplace_front() {
  if (empty()) {
    insert_for_empty_list();
  } else {
    insert_before(m_links_data.front_index());
  }
  return m_links_data.front_index();
}

template <class T>
void VectorListHybrid<T>::pop_front() {
  erase(m_links_data.front_index());
}

template <class T>
typename VectorListHybrid<T>::ID VectorListHybrid<T>::insert_for_empty_list() {
  m_links_data.insert_for_empty_list();
  return get_checked_new_id(m_links_data.front_index());
}

template <class T>
typename VectorListHybrid<T>::ID VectorListHybrid<T>::insert_after(
    VectorListHybrid<T>::ID id) {
  m_links_data.insert_after(id);
  return get_checked_new_id(m_links_data.next(id));
}

template <class T>
typename VectorListHybrid<T>::ID VectorListHybrid<T>::insert_before(
    VectorListHybrid<T>::ID id) {
  m_links_data.insert_before(id);
  return get_checked_new_id(m_links_data.previous(id));
}

template <class T>
T& VectorListHybrid<T>::back() {
  return m_data[m_links_data.back_index()];
}

template <class T>
T& VectorListHybrid<T>::front() {
  return m_data[m_links_data.front_index()];
}

template <class T>
T& VectorListHybrid<T>::at(ID id) {
  return m_data[id];
}

template <class T>
const T& VectorListHybrid<T>::at(ID id) const {
  return m_data[id];
}

template <class T>
std::optional<typename VectorListHybrid<T>::ID> VectorListHybrid<T>::next(
    ID id) const {
  const ID index = m_links_data.next(id);
  return optional_id(index);
}

template <class T>
std::optional<typename VectorListHybrid<T>::ID> VectorListHybrid<T>::next(
    std::optional<ID> id) const {
  return next(id.value());
}

template <class T>
std::optional<typename VectorListHybrid<T>::ID> VectorListHybrid<T>::previous(
    ID id) const {
  return optional_id(m_links_data.previous(id));
}

template <class T>
std::optional<typename VectorListHybrid<T>::ID> VectorListHybrid<T>::back_id()
    const {
  return optional_id(m_links_data.back_index());
}

template <class T>
std::optional<typename VectorListHybrid<T>::ID> VectorListHybrid<T>::front_id()
    const {
  return optional_id(m_links_data.front_index());
}

template <class T>
void VectorListHybrid<T>::erase(ID id) {
  m_links_data.erase(id);
}

template <class T>
void VectorListHybrid<T>::erase_interval(
    typename VectorListHybrid<T>::ID id, size_t number_of_elements) {
  m_links_data.erase_interval(id, number_of_elements);
}

template <class T>
template <class CIter>
OverwriteIntervalResult VectorListHybrid<T>::overwrite_interval(
    typename VectorListHybrid<T>::ID id, const CIter& new_elements_cbegin,
    const CIter& new_elements_cend) {
  // The links are unchanged; only the elements need to be changed.
  OverwriteIntervalResult result;
  result.final_overwritten_element_id = id;
  CIter citer = new_elements_cbegin;
  TKET_ASSERT(citer != new_elements_cend);
  const auto max_number_of_elements = m_links_data.size();
  result.number_of_overwritten_elements = 0;
  for (;;) {
    m_data.at(result.final_overwritten_element_id) = *citer;
    ++result.number_of_overwritten_elements;
    // GCOVR_EXCL_START
    TKET_ASSERT(
        result.number_of_overwritten_elements <= max_number_of_elements);
    // GCOVR_EXCL_STOP
    ++citer;
    if (citer == new_elements_cend) {
      return result;
    }
    // There IS another element, where will it be overwritten?
    result.final_overwritten_element_id =
        m_links_data.next(result.final_overwritten_element_id);
  }
  // Should be impossible to reach here
  TKET_ASSERT(false);
}

template <class T>
void VectorListHybrid<T>::append_to_vector(std::vector<T>& vect) const {
  vect.reserve(vect.size() + size());
  for (ID current_index = m_links_data.front_index();
       current_index != m_links_data.get_invalid_index();
       current_index = m_links_data.next(current_index)) {
    vect.emplace_back(m_data[current_index]);
  }
}

template <class T>
std::vector<T> VectorListHybrid<T>::to_vector() const {
  std::vector<T> result;
  append_to_vector(result);
  return result;
}

template <class T>
typename VectorListHybrid<T>::ID VectorListHybrid<T>::get_checked_new_id(
    ID id) {
  if (m_data.size() <= id) {
    m_data.resize(id + 1);
  }
  return id;
}

template <class T>
std::string VectorListHybrid<T>::debug_str() const {
  std::stringstream ss;
  ss << "\nRaw stored elems:";
  for (size_t nn = 0; nn < m_data.size(); ++nn) {
    ss << "\nData[" << nn << "] = " << m_data[nn];
  }
  ss << "\n" << m_links_data.debug_str() << "\n";
  return ss.str();
}

}  // namespace tket
