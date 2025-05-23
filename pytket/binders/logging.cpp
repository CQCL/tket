// Copyright Quantinuum
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

#include <nanobind/nanobind.h>

#include <tklog/TketLog.hpp>

namespace nb = nanobind;

namespace tket {

NB_MODULE(logging, m) {
  nb::set_leak_warnings(false);
  nb::enum_<LogLevel>(m, "level")
      .value("trace", LogLevel::Trace, "all logs")
      .value("debug", LogLevel::Debug, "debug logs and above")
      .value("info", LogLevel::Info, "informational logs and above")
      .value("warn", LogLevel::Warn, "warnings and above")
      .value("err", LogLevel::Err, "error and critical logs only (default)")
      .value("critical", LogLevel::Critical, "critical logs only")
      .value("off", LogLevel::Off, "no logs");

  m.def(
      "set_level", [](LogLevel level) { tket_log()->set_level(level); },
      "Set the global logging level."
      "\n\n:param log_level: Desired logging level",
      nb::arg("log_level"));
}

}  // namespace tket
