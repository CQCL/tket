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

/**
 * @file
 * @brief Logging
 */

#include <iostream>
#include <memory>
#include <string>

namespace tket {

enum class LogLevel {
  Trace = 0,
  Debug = 1,
  Info = 2,
  Warn = 3,
  Err = 4,
  Critical = 5,
  Off = 6
};

class Logger {
 public:
  Logger(LogLevel level = LogLevel::Err);
  void trace(const std::string &s, std::ostream &os = std::cout);
  void debug(const std::string &s, std::ostream &os = std::cout);
  void info(const std::string &s, std::ostream &os = std::cout);
  void warn(const std::string &s, std::ostream &os = std::cout);
  void error(const std::string &s, std::ostream &os = std::cout);
  void critical(const std::string &s, std::ostream &os = std::cout);
  void set_level(LogLevel lev);

 private:
  LogLevel level;
  void log(const char *levstr, const std::string &s, std::ostream &os);
};

typedef std::shared_ptr<Logger> LogPtr_t;

/** Logger for tket messages */
LogPtr_t &tket_log();

}  // namespace tket
