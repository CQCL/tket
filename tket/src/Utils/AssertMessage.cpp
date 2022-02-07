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

std::string AssertMessage::get_error_message() {
  const auto message = get_error_stream().str();

  // Clear the global stream, ready for the next message
  // (currently this isn't necessary, because tket assert
  // immediately aborts; but it may become necessary again in future,
  // if we have assert variants with throws and multiple try/catch).
  get_error_stream().str(std::string());
  return message;
}

AssertMessage::operator bool() const { return false; }

std::stringstream& AssertMessage::get_error_stream() {
  static std::stringstream ss;
  return ss;
}
// GCOVR_EXCL_STOP

}  // namespace tket
