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

namespace tket {

// GCOVR_EXCL_START
/** This is only for use with TKET_ASSERT.
 */
class AssertMessage {
 public:
  /** Construct the object, to begin writing to the stream. */
  AssertMessage();

  /** Always returns false, so that "|| AssertMessage() << a)" becomes
   * "|| false)".
   * Also, stores the error message for later use by TKET_ASSERT macros;
   * previously this information was passed on by exceptions, but that
   * generated lots of code coverage branching problems. */
  operator bool() const;

  /** Every streamable object x can be written to the stream.
   * @param x Any object which can be written to a stringstream.
   * @return This object, to allow chaining.
   */
  template <class T>
  AssertMessage& operator<<(const T& x) {
    m_ss << x;
    return *this;
  }

  /** Get the stored error message. Of course, if AssertMessage()
   * has not actually been called, just returns an empty string.
   * Also, clears the stored message, ready for the next time.
   */
  static std::string get_error_message();

 private:
  std::stringstream m_ss;

  static std::string& get_error_message_ref();
};
// GCOVR_EXCL_STOP

}  // namespace tket
