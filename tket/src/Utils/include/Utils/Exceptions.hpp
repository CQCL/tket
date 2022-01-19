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

#include <exception>
#include <stdexcept>
#include <string>

namespace tket {

/** Operation type not supported */
class NotImplemented : public std::logic_error {
 public:
  explicit NotImplemented(const std::string &message)
      : std::logic_error(message) {}
};

/** Operation invalid */
class NotValid : public std::logic_error {
 public:
  NotValid() : std::logic_error("Not a valid operation") {}
  explicit NotValid(const std::string &message) : std::logic_error(message) {}
};

/** Matrix not unitary */
class NotUnitary : public std::logic_error {
 public:
  NotUnitary() : std::logic_error("Not a unitary matrix") {}
  explicit NotUnitary(const std::string &message) : std::logic_error(message) {}
};

}  // namespace tket
