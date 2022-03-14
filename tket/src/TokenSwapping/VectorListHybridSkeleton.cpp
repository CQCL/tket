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

#include "VectorListHybridSkeleton.hpp"

#include <limits>
#include <sstream>
#include <stdexcept>

#include "Utils/Assert.hpp"

namespace tket {
namespace tsa_internal {

using Index = VectorListHybridSkeleton::Index;

const Index INVALID_INDEX = std::numeric_limits<Index>::max();

Index VectorListHybridSkeleton::get_invalid_index() { return INVALID_INDEX; }

VectorListHybridSkeleton::VectorListHybridSkeleton()
    : m_size(0),
      m_front(INVALID_INDEX),
      m_back(INVALID_INDEX),
      m_deleted_front(INVALID_INDEX) {}

void VectorListHybridSkeleton::clear() {
  if (m_links.empty()) {
    TKET_ASSERT(m_size == 0);
    TKET_ASSERT(m_front == INVALID_INDEX);
    TKET_ASSERT(m_back == INVALID_INDEX);
    TKET_ASSERT(m_deleted_front == INVALID_INDEX);
    return;
  }
  m_size = 0;
  m_front = INVALID_INDEX;
  m_back = INVALID_INDEX;
  for (Index nn = 1; nn < m_links.size(); ++nn) {
    // Not strictly necessary, as deleted links are only a forward list;
    // but make absolutely sure no leakage of prior internal link data can
    // occur.
    m_links[nn].previous = nn - 1;
    m_links[nn - 1].next = nn;
  }
  m_links[0].previous = INVALID_INDEX;
  m_links.back().next = INVALID_INDEX;
  m_deleted_front = 0;
}

void VectorListHybridSkeleton::fast_clear() {
  if (m_back == INVALID_INDEX) {
    // No elements stored currently; nothing to do.
    TKET_ASSERT(m_size == 0);
    TKET_ASSERT(m_front == INVALID_INDEX);
    return;
  }
  TKET_ASSERT(m_size > 0);
  TKET_ASSERT(m_front != INVALID_INDEX);
  TKET_ASSERT(m_links[m_back].next == INVALID_INDEX);
  // There are some existing elements.
  // Recall that deleted elements are ONLY a forward list,
  // so we don't need to update "previous".
  // To combine existing active elements with
  // existing deleted elements,
  // the valid elements will be joined to
  // the start of the deleted list.
  if (m_deleted_front != INVALID_INDEX) {
    m_links[m_back].next = m_deleted_front;
  }
  // Convert "active" elements into deleted elements.
  m_deleted_front = m_front;
  m_front = INVALID_INDEX;
  m_back = INVALID_INDEX;
  m_size = 0;
}

void VectorListHybridSkeleton::reverse() {
  if (m_size <= 1) {
    // Nothing to do.
    return;
  }
  TKET_ASSERT(m_front != INVALID_INDEX);
  TKET_ASSERT(m_back != INVALID_INDEX);
  TKET_ASSERT(m_front != m_back);
  // The deleted element links don't need to change.
  {
    auto current_index = m_front;
    bool terminated_correctly = false;
    for (auto infinite_loop_guard = 1 + m_links.size(); infinite_loop_guard > 0;
         --infinite_loop_guard) {
      auto& link = m_links[current_index];
      const auto next_index = link.next;
      std::swap(link.next, link.previous);
      if (next_index >= m_links.size()) {
        TKET_ASSERT(next_index == INVALID_INDEX);
        terminated_correctly = true;
        break;
      }
      current_index = next_index;
    }
    TKET_ASSERT(terminated_correctly);
  }
  std::swap(m_front, m_back);
}

size_t VectorListHybridSkeleton::size() const { return m_size; }

Index VectorListHybridSkeleton::front_index() const { return m_front; }

Index VectorListHybridSkeleton::back_index() const { return m_back; }

Index VectorListHybridSkeleton::next(Index index) const {
  return m_links[index].next;
}

Index VectorListHybridSkeleton::previous(Index index) const {
  return m_links[index].previous;
}

void VectorListHybridSkeleton::erase(Index index) {
  --m_size;
  auto& current_link = m_links[index];
  if (current_link.previous == INVALID_INDEX) {
    // We're erasing the front.
    m_front = current_link.next;
  } else {
    m_links[current_link.previous].next = current_link.next;
  }
  if (current_link.next == INVALID_INDEX) {
    // We're erasing the back.
    m_back = current_link.previous;
  } else {
    m_links[current_link.next].previous = current_link.previous;
  }
  // Recall: deleted elements are a forward list ONLY.
  current_link.next = m_deleted_front;
  m_deleted_front = index;
}

void VectorListHybridSkeleton::erase_interval(
    Index index, size_t number_of_elements) {
  if (number_of_elements == 0) {
    return;
  }
  // First, find the index of the LAST element to be erased.
  // Notice that this is the only O(N) part; the rest are O(1).
  // We update only O(1) links in total, not O(N),
  // so slightly faster than a loop of next/erase calls.
  Index last_element_index = index;
  for (size_t nn = 1; nn < number_of_elements; ++nn) {
    last_element_index = m_links.at(last_element_index).next;

    // GCOVR_EXCL_START
    TKET_ASSERT(
        last_element_index < m_links.size() ||
        AssertMessage() << "erase_interval with start index " << index
                        << ", number_of_elements=" << number_of_elements
                        << ", size " << m_links.size()
                        << ", runs out of elements at N=" << nn
                        << " (got index " << last_element_index << ")");
    // GCOVR_EXCL_STOP
  }
  TKET_ASSERT(number_of_elements <= m_size);
  m_size -= number_of_elements;

  // Now, splice the soon-to-be-logically-erased interval into the deleted
  // elements. Start the new deleted list at the erased interval.
  const auto index_of_node_after_interval = m_links[last_element_index].next;

  // Correct whether or not m_deleted_front equals INVALID_INDEX.
  m_links[last_element_index].next = m_deleted_front;
  // No need to update previous, since the deleted nodes are only a forward
  // list.
  m_deleted_front = index;

  // Link the node BEFORE the interval to the new next node.
  const auto index_of_node_before_interval = m_links[index].previous;

  if (index_of_node_before_interval < m_links.size()) {
    // There IS a previous node to be dealt with.
    auto& next_node_index_ref = m_links[index_of_node_before_interval].next;
    TKET_ASSERT(next_node_index_ref == index);
    // This is correct even if index_of_node_after_interval is INVALID_INDEX.
    next_node_index_ref = index_of_node_after_interval;
    TKET_ASSERT(m_front != index);
  } else {
    // No previous node, we must have been at the start already.
    TKET_ASSERT(index_of_node_before_interval == INVALID_INDEX);
    TKET_ASSERT(m_front == index);
    m_front = index_of_node_after_interval;
  }
  // Link the node AFTER the interval to the new previous node.
  if (index_of_node_after_interval < m_links.size()) {
    // There are more unerased elements after the interval,
    // so the first one must be dealt with.
    auto& prev_node_index = m_links[index_of_node_after_interval].previous;
    TKET_ASSERT(prev_node_index == last_element_index);
    // Correct even if there IS no node before the interval.
    prev_node_index = index_of_node_before_interval;
    TKET_ASSERT(m_back != last_element_index);
  } else {
    // No node after, we have erased up to the back.
    TKET_ASSERT(index_of_node_after_interval == INVALID_INDEX);
    TKET_ASSERT(m_back == last_element_index);
    m_back = index_of_node_before_interval;
  }
  if (m_size == 0) {
    TKET_ASSERT(m_front == INVALID_INDEX);
    TKET_ASSERT(m_back == INVALID_INDEX);
  } else {
    TKET_ASSERT(m_front < m_links.size());
    TKET_ASSERT(m_back < m_links.size());
    if (m_size == 1) {
      TKET_ASSERT(m_front == m_back);
    }
  }
}

void VectorListHybridSkeleton::insert_for_empty_list() {
  const auto new_index = get_new_index();
  m_front = new_index;
  m_back = new_index;
  m_links[new_index].next = INVALID_INDEX;
  m_links[new_index].previous = INVALID_INDEX;
}

void VectorListHybridSkeleton::insert_after(Index index) {
  const auto new_index = get_new_index();
  const auto old_next = m_links[index].next;
  m_links[index].next = new_index;
  m_links[new_index].next = old_next;
  m_links[new_index].previous = index;
  if (old_next == INVALID_INDEX) {
    // The old element was already at the back.
    m_back = new_index;
  } else {
    m_links[old_next].previous = new_index;
  }
}

void VectorListHybridSkeleton::insert_before(Index index) {
  const auto new_index = get_new_index();
  const auto old_prev = m_links[index].previous;
  m_links[index].previous = new_index;
  m_links[new_index].next = index;
  m_links[new_index].previous = old_prev;
  if (old_prev == INVALID_INDEX) {
    // The old element was already at the front.
    m_front = new_index;
  } else {
    m_links[old_prev].next = new_index;
  }
}

Index VectorListHybridSkeleton::get_new_index() {
  ++m_size;
  if (m_deleted_front == INVALID_INDEX) {
    // We need to create a new element, it's full.
    m_links.emplace_back();
    return m_links.size() - 1;
  }
  // Reuse a deleted element.
  const auto old_deleted_front = m_deleted_front;
  m_deleted_front = m_links[old_deleted_front].next;
  return old_deleted_front;
}

std::string VectorListHybridSkeleton::debug_str() const {
  std::stringstream ss;
  const auto to_str = [](size_t ii) -> std::string {
    if (ii == INVALID_INDEX) {
      return "NULL";
    }
    return std::to_string(ii);
  };

  ss << "VLHS: size " << m_size << ", front " << to_str(m_front) << " back "
     << to_str(m_back) << ", del.front " << to_str(m_deleted_front);

  ss << "\nActive links: forward [";
  for (auto index = m_front; index != INVALID_INDEX;
       index = m_links[index].next) {
    ss << index << "->";
  }
  ss << "]\nBackward (";
  for (auto index = m_back; index != INVALID_INDEX;
       index = m_links[index].previous) {
    ss << index << "->";
  }
  ss << ")\nDel.links: {";
  for (auto index = m_deleted_front; index != INVALID_INDEX;
       index = m_links[index].next) {
    ss << index << "->";
  }
  ss << "}";
  return ss.str();
}

}  // namespace tsa_internal
}  // namespace tket
