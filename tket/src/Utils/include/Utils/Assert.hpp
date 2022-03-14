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

#include <cstdlib>

#include "AssertMessage.hpp"
#include "TketLog.hpp"

/**
 * If `condition` is not satisfied, log a diagnostic message and abort.
 * You can abort with a fixed string:
 *
 *    TKET_ASSERT(!"Some error message...");
 *
 * For a simple statement like:
 *
 *    TKET_ASSERT(x<y);
 *
 * ...the final message will include the raw C++ code "x<y", but not
 * the actual values of x,y. To include more information in the message,
 * use the special "AssertMessage" class, like so:
 *
 *    TKET_ASSERT(x<y ||
 *        AssertMessage() << "The values are x=" << x << ", y=" << y);
 *
 * You can have more complicated conditions:
 *
 *    TKET_ASSERT(x<y || y<z || x<z ||
 *        AssertMessage() << "The values are " << x << "," << y << "," << z);
 *
 * ...but ensure that AssertMessage() occurs only at the end,
 * always OR-ed with the condition you're asserting.
 *
 * You can also abort without any condition, with a dynamically
 * constructed runtime message:
 *
 *    TKET_ASSERT(AssertMessage() << "Error: x=" << x
 *          << ", y=" << y << ", sum=" << x+y);
 *
 * The message construction with AssertMessage() will NOT begin
 * if `condition` is true (short circuit evaluation),
 * so there is no performance penalty.
 *
 * This also checks if evaluating `condition` itself throws an exception.
 *
 * Notes: (1) single line code is ignored by test code coverage.
 * However, multiline code is not. Currently we just manually surround the
 * worst multiline offenders with start/stop tags.
 *
 * (2) We tried putting the start/stop tags within the macro, to make test
 * coverage ignore the code; unfortunately this did NOT work.
 * We suspect that it's because comments are stripped from the macro before
 * test coverage sees the code, but we don't know.
 *
 * (3) It is known that exceptions cause problems by generating numerous
 * extra branches:
 * https://stackoverflow.com/questions/42003783/
 * lcov-gcov-branch-coverage-with-c-producing-branches-all-over-the-place?rq=1
 * Removing exceptions (or "hiding" them, by tricks)
 * did cut down the number of extra branches, but did not remove them
 * completely. Exceptions are not the sole cause of branching problems!
 */
#define TKET_ASSERT(condition)                                          \
  do {                                                                  \
    try {                                                               \
      if (!(condition)) {                                               \
        std::stringstream ss;                                           \
        ss << "Assertion '" << #condition << "' (" << __FILE__ << " : " \
           << __func__ << " : " << __LINE__ << ") failed. "             \
           << AssertMessage::get_error_message() << " Aborting.";       \
        tket::tket_log()->critical(ss.str());                           \
        std::abort();                                                   \
      }                                                                 \
    } catch (const std::exception& ex) {                                \
      std::stringstream ss;                                             \
      ss << "Evaluating assertion condition '" << #condition << "' ("   \
         << __FILE__ << " : " << __func__ << " : " << __LINE__          \
         << ") threw unexpected exception: '" << ex.what() << "'. "     \
         << AssertMessage::get_error_message() << " Aborting.";         \
      tket::tket_log()->critical(ss.str());                             \
      std::abort();                                                     \
    } catch (...) {                                                     \
      std::stringstream ss;                                             \
      ss << "Evaluating assertion condition '" << #condition << "' ("   \
         << __FILE__ << " : " << __func__ << " : " << __LINE__          \
         << ") Threw unknown exception. "                               \
         << AssertMessage::get_error_message() << " Aborting.";         \
      tket::tket_log()->critical(ss.str());                             \
      std::abort();                                                     \
    }                                                                   \
  } while (0)
