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

#include "AssertWithThrowHelper.hpp"

// GCOVR_EXCL_START
AssertWithThrowHelper::AssertWithThrowHelper()
    : m_get_error_stream_called(false) {}

std::stringstream& AssertWithThrowHelper::get_error_stream() {
  auto& object = get();
  object.m_get_error_stream_called = true;
  return object.m_ss;
}

void AssertWithThrowHelper::throw_upon_error() {
  auto& object = get();
  if (!object.m_get_error_stream_called) {
    return;
  }
  const auto message = object.m_ss.str();
  // Clear the stream, ready for the next error
  // (since the caller might catch the expection and then throw others).
  object.m_get_error_stream_called = false;
  object.m_ss.str(std::string());
  throw std::runtime_error(message);
}

AssertWithThrowHelper& AssertWithThrowHelper::get() {
  static AssertWithThrowHelper object;
  return object;
}
// GCOVR_EXCL_STOP
