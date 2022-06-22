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

#include "AssertMessage.hpp"

namespace tket {

// GCOVR_EXCL_START
AssertMessage::AssertMessage() {}

std::string& AssertMessage::get_error_message_ref() {
  static std::string error_string;
  return error_string;
}

std::string AssertMessage::get_error_message() {
  const std::string message = get_error_message_ref();
  // Asserts are SUPPOSED to lead to aborts, so clearing
  // shouldn't be necessary; but anyway, in case it's
  // called multiple times, clear ready for the next message.
  get_error_message_ref().clear();
  return message;
}

AssertMessage::operator bool() const {
  // Store the built up error message.
  get_error_message_ref() = m_ss.str();
  return false;
}
// GCOVR_EXCL_STOP

}  // namespace tket
