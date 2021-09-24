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

#include <pybind11/pybind11.h>

#include "Utils/TketLog.hpp"

namespace py = pybind11;

namespace tket {

PYBIND11_MODULE(logging, m) {
  py::enum_<spdlog::level::level_enum>(m, "level")
      .value("trace", spdlog::level::trace, "all logs")
      .value("debug", spdlog::level::debug, "debug logs and above")
      .value("info", spdlog::level::info, "informational logs and above")
      .value("warn", spdlog::level::warn, "warnings and above")
      .value(
          "err", spdlog::level::err, "error and critical logs only (default)")
      .value("critical", spdlog::level::critical, "critical logs only")
      .value("off", spdlog::level::off, "no logs");

  m.def(
      "set_level",
      [](spdlog::level::level_enum level) { tket_log()->set_level(level); },
      "Set the global logging level."
      "\n\n:param level: Desired logging level");
}

}  // namespace tket
