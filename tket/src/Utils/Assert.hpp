// Copyright 2019-2021 Cambridge Quantum Computing
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

#include "AssertMessage.hpp"
#include "TketLog.hpp"

/**
 * If the condition `b` is not satisfied, log a diagnostic message and abort.
 * But note that the message includes only the raw C++ source code for b,
 * not the actual values of x,y in conditions like "x<y".
 * You can construct a dynamic error message
 * with extra information as follows:
 *
 * TKET_ASSERT(x<y ||
 *    AssertMessage() << "x=" << x << ", y=" << y << " and also z=" << z);
 *
 * A verbose version including the actual C++ source code
 * in the TKET_ASSERT statement with "AssertMessage", as well as
 * the constructed error message itself, can be achieved with
 * AssertMessage::verbose() instead of AssertMessage().
 * The message construction will NOT begin if `b` is true,
 * so there is no performance penalty.
 *
 * The code can be multiline, but should still be ignored by test code coverage.
 *
 * This should also work if evaluating the assertion condition
 * itself throws an exception.
 */
#define TKET_ASSERT(b)                                                         \
  /* GCOVR_EXCL_START */                                                       \
  do {                                                                         \
    try {                                                                      \
      if (!(b)) {                                                              \
        std::stringstream msg;                                                 \
        msg << "Assertion '" << #b << "' (" << __FILE__ << " : " << __func__   \
            << " : " << __LINE__ << ") failed: aborting.";                     \
        tket::tket_log()->critical(msg.str());                                 \
        std::abort();                                                          \
      }                                                                        \
    } catch (const AssertMessage::MessageData& e1) {                           \
      std::stringstream msg;                                                   \
      msg << "Assertion ";                                                     \
      if (e1.verbose) {                                                        \
        msg << "'" << #b << "' ";                                              \
      }                                                                        \
      msg << "(" << __FILE__ << " : " << __func__ << " : " << __LINE__         \
          << ") failed: '" << e1.what() << "': aborting.";                     \
      tket::tket_log()->critical(msg.str());                                   \
      std::abort();                                                            \
    } catch (const std::exception& e2) {                                       \
      std::stringstream msg;                                                   \
      msg << "Evaluating assertion condition '" << #b << "' (" << __FILE__     \
          << " : " << __func__ << " : " << __LINE__                            \
          << ") threw unexpected exception: '" << e2.what() << "': aborting."; \
      tket::tket_log()->critical(msg.str());                                   \
      std::abort();                                                            \
    } catch (...) {                                                            \
      std::stringstream msg;                                                   \
      msg << "Evaluating assertion condition '" << #b << "' (" << __FILE__     \
          << " : " << __func__ << " : " << __LINE__                            \
          << ") threw unknown exception. Aborting.";                           \
      tket::tket_log()->critical(msg.str());                                   \
      std::abort();                                                            \
    }                                                                          \
  } while (0) /* GCOVR_EXCL_STOP */
