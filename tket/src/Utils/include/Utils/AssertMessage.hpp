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

#include <sstream>
#include <stdexcept>

namespace tket {

/** This is for use with TKET_ASSERT, when we want to give a more detailed
 * error message than just the assertion code and location.
 */
class AssertMessage {
 public:
  /** Construct the object (the default non-verbose version) to begin writing to
   * the stream. */
  AssertMessage();

  /** Get a verbose object (not the default). */
  static AssertMessage verbose();

  /** Thrown when the message construction is finished, to store the necessary
   * data. */
  struct MessageData : public std::runtime_error {
    bool verbose;
    MessageData(const std::string& str, bool verbose);
  };

  /** Throws a MessageData object when called, with the message. */
  operator bool() const;

  /** Every streamable object can be written to the stream. */
  template <class T>
  AssertMessage& operator<<(const T& x) {
    m_ss << x;
    return *this;
  }

 private:
  bool m_verbose;
  std::stringstream m_ss;
};

}  // namespace tket
