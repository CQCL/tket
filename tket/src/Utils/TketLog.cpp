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

#include "TketLog.hpp"

namespace tket {
LogPtr_t &tket_log() {
  static LogPtr_t logger = [](LogPtr_t l) -> LogPtr_t {
    l->set_pattern("%+");
#ifdef ALL_LOGS
    l->set_level(spdlog::level::trace);
#else
    l->set_level(spdlog::level::err);
#endif
    return l;
  }((LogPtr_t)spdlog::stdout_color_mt("tket"));
  return logger;
}
}  // namespace tket
