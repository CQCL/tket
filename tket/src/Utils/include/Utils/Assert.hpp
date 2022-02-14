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
#include <sstream>

#include "TketLog.hpp"

/**
 * If `condition` is not satisfied, log a diagnostic message and abort,
 * including the extra message "msg".
 * "msg" is passed directly to a stringstream, so you can write:
 *
 *  TKET_ASSERT_WITH_MESSAGE(x<y, "The values are x=" << x << ", y=" << y);
 *
 * Note that the macro would automatically generate the message "x<y",
 * i.e. the raw C++ code defining the condition.
 *
 * Note that the message construction (including streaming the second argument
 * to a stringstream) will NOT begin if `condition` is true,
 * so there is no performance penalty.
 *
 * This also checks if evaluating `condition` itself throws an exception.
 *
 * Notes: (1) single line code is ignored by test code coverage.
 * However, multiline code is not. Currently we just manually surround the
 * worst multiline offenders with start/stop tags.
 * (2) We tried putting the start/stop tags within the macro, to make test
 * coverage ignore the code; unfortunately this did NOT work.
 * We suspect that it's because comments are stripped from the macro before
 * test coverage sees the code, but we don't know. See TKET-1856.
 * (3) It is known that exceptions cause problems by generating numerous
 * extra branches:
 * https://stackoverflow.com/questions/42003783/
 * lcov-gcov-branch-coverage-with-c-producing-branches-all-over-the-place?rq=1
 * However, although removing exceptions (or "hiding" them, by tricks)
 * did cut down the number of extra branches, it did not remove them
 * completely, so exceptions are not the sole cause of branching problems.
 */
#define TKET_ASSERT_WITH_MESSAGE(condition, msg)                           \
  do {                                                                     \
    try {                                                                  \
      if (!(condition)) {                                                  \
        std::stringstream ss;                                              \
        ss << "Assertion '" << #condition << "' (" << __FILE__ << " : "    \
           << __func__ << " : " << __LINE__ << ") failed. " << msg         \
           << " Aborting.";                                                \
        tket::tket_log()->critical(ss.str());                              \
        std::abort();                                                      \
      }                                                                    \
    } catch (const std::exception& ex) {                                   \
      std::stringstream ss;                                                \
      ss << "Evaluating assertion condition '" << #condition << "' ("      \
         << __FILE__ << " : " << __func__ << " : " << __LINE__             \
         << ") threw unexpected exception: '" << ex.what() << "'. " << msg \
         << " Aborting.";                                                  \
      tket::tket_log()->critical(ss.str());                                \
      std::abort();                                                        \
    } catch (...) {                                                        \
      std::stringstream ss;                                                \
      ss << "Evaluating assertion condition '" << #condition << "' ("      \
         << __FILE__ << " : " << __func__ << " : " << __LINE__             \
         << ") Threw unknown exception. " << msg << " Aborting.";          \
      tket::tket_log()->critical(ss.str());                                \
      std::abort();                                                        \
    }                                                                      \
  } while (0)

#define TKET_ASSERT(condition)               \
  do {                                       \
    TKET_ASSERT_WITH_MESSAGE(condition, ""); \
  } while (0)
