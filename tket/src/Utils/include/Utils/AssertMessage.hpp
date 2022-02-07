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
/** This is only for use with TKET_ASSERT, when we want to give a more detailed
 * error message than just the assertion code and location.
 * Also, some code might seem strange, but that's because exceptions
 * can generate many extra branches in test coverage, see
 *
 * https://stackoverflow.com/questions/42003783/
 * lcov-gcov-branch-coverage-with-c-producing-branches-all-over-the-place?rq=1
 *
 * Thus, we avoid exceptions.
 */
class AssertMessage {
 public:
  /** Construct the object, to begin writing to the stream. */
  AssertMessage();

  /** Always returns false, so that "... || AssertMessage() << a)"
   * becomes "... || false)".
   */
  operator bool() const;

  /** Every streamable object x can be written to the stream. */
  template <class T>
  const AssertMessage& operator<<(const T& x) const {
    get_error_stream() << x;
    return *this;
  }

  /** Get the stored error message. Of course, if AssertMessage()
   * has not actually been called, just returns an empty string.
   * Also, clears the stored message, ready for the next time.
   */
  static std::string get_error_message();

 private:
  /** Previously the error message for later use by TKET_ASSERT macros
   * was passed on by exceptions within operator bool(), but that
   * generated lots of code coverage branching problems.
   * So now we use a global variable. The AssertMessage object
   * will go out of scope, so there seems to be no other good way
   * to pass the information on.
   */
  static std::stringstream& get_error_stream();
};
// GCOVR_EXCL_STOP

}  // namespace tket
