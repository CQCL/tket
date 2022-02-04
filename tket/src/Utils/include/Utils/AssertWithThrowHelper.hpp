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

// GCOVR_EXCL_START
/** To be used only by the TKET_ASSERT_WITH_THROW macro.
 * Something like this is necessary to prevent exceptions
 * generating many extra branches in test coverage, see
 *
 * https://stackoverflow.com/questions/42003783/
 * lcov-gcov-branch-coverage-with-c-producing-branches-all-over-the-place?rq=1
 *
 * We want to hide the throws (or at least, have one single throw),
 * and also provide a stringstream to avoid having to construct one
 * every time the assert is checked (asserts must have
 * almost zero performance impact if they are not triggered).
 */
class AssertWithThrowHelper {
 public:
  /** Get a stored stream, to write errors to.
   * The caller should only call this if they are certain
   * that an error has occurred.
   * The caller can write to this multiple times
   * before calling throw_upon_error().
   */
  static std::stringstream& get_error_stream();

  /** If get_error_stream() was previously called,
   * throw an exception with the contents of the stream
   * (even if an empty string),
   * and clear the stream ready for the next use.
   * Otherwise does nothing.
   */
  static void throw_upon_error();

 private:
  bool m_get_error_stream_called;
  std::stringstream m_ss;

  AssertWithThrowHelper();

  static AssertWithThrowHelper& get();
};
// GCOVR_EXCL_STOP
