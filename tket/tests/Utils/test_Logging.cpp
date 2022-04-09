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

#include <Utils/TketLog.hpp>
#include <catch2/catch.hpp>

namespace tket {

SCENARIO("Logging") {
  tket_log()->set_level(LogLevel::Trace);

  tket_log()->trace("Test trace (output).");
  tket_log()->debug("Test debug (output).");
  tket_log()->info("Test info (output).");
  tket_log()->warn("Test warn (output).");
  tket_log()->error("Test error (output).");
  tket_log()->critical("Test critical (output).");

  tket_log()->set_level(LogLevel::Off);

  tket_log()->trace("Test trace (no output).");
  tket_log()->debug("Test debug (no output).");
  tket_log()->info("Test info (no output).");
  tket_log()->warn("Test warn (no output).");
  tket_log()->error("Test error (no output).");
  tket_log()->critical("Test critical (no output).");

  tket_log()->set_level(LogLevel::Err);
}

}  // namespace tket
