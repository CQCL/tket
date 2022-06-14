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
#include <catch2/catch_test_macros.hpp>
#include <sstream>
#include <string>

namespace tket {

SCENARIO("Logging") {
  std::stringstream ss;

  tket_log()->set_level(LogLevel::Trace);

  tket_log()->trace("Test trace (output).", ss);
  tket_log()->debug("Test debug (output).", ss);
  tket_log()->info("Test info (output).", ss);
  tket_log()->warn("Test warn (output).", ss);
  tket_log()->error("Test error (output).", ss);
  tket_log()->critical("Test critical (output).", ss);

  tket_log()->set_level(LogLevel::Off);

  tket_log()->trace("Test trace (no output).", ss);
  tket_log()->debug("Test debug (no output).", ss);
  tket_log()->info("Test info (no output).", ss);
  tket_log()->warn("Test warn (no output).", ss);
  tket_log()->error("Test error (no output).", ss);
  tket_log()->critical("Test critical (no output).", ss);

  tket_log()->set_level(LogLevel::Err);

  std::string s = ss.str();

  CHECK(s.find("Test trace (output).") != std::string::npos);
  CHECK(s.find("Test debug (output).") != std::string::npos);
  CHECK(s.find("Test info (output).") != std::string::npos);
  CHECK(s.find("Test warn (output).") != std::string::npos);
  CHECK(s.find("Test error (output).") != std::string::npos);
  CHECK(s.find("Test critical (output).") != std::string::npos);
  CHECK(s.find("(no output)") == std::string::npos);
}

}  // namespace tket
