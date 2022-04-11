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

#include "TketLog.hpp"

#include <ctime>
#include <iomanip>
#include <iostream>
#include <memory>

namespace tket {

Logger::Logger(LogLevel level) : level(level) {}
void Logger::trace(const std::string &s, std::ostream &os) {
  if (level <= LogLevel::Trace) {
    log("trace", s, os);
  }
}
void Logger::debug(const std::string &s, std::ostream &os) {
  if (level <= LogLevel::Debug) {
    log("debug", s, os);
  }
}
void Logger::info(const std::string &s, std::ostream &os) {
  if (level <= LogLevel::Info) {
    log("info", s, os);
  }
}
void Logger::warn(const std::string &s, std::ostream &os) {
  if (level <= LogLevel::Warn) {
    log("warn", s, os);
  }
}
void Logger::error(const std::string &s, std::ostream &os) {
  if (level <= LogLevel::Err) {
    log("error", s, os);
  }
}
void Logger::critical(const std::string &s, std::ostream &os) {
  if (level <= LogLevel::Critical) {
    log("critical", s, os);
  }
}

void Logger::log(const char *levstr, const std::string &s, std::ostream &os) {
  std::time_t t = std::time(nullptr);
  os << "[" << std::put_time(std::localtime(&t), "%Y-%m-%d %H:%M:%S")
     << "] [tket] [" << levstr << "] " << s << std::endl;
}

void Logger::set_level(LogLevel lev) { level = lev; }

LogPtr_t &tket_log() {
#ifdef ALL_LOGS
  static LogPtr_t logger = std::make_shared<Logger>(LogLevel::Trace);
#else
  static LogPtr_t logger = std::make_shared<Logger>(LogLevel::Err);
#endif
  return logger;
}
}  // namespace tket
