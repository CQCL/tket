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

#include "DAGProperties.hpp"

#include <algorithm>

#include "OpType/EdgeType.hpp"
#include "Utils/GraphHeaders.hpp"
#include "Utils/TketLog.hpp"

namespace tket {

#ifdef CHECK
#error "Macro already defined!"
#endif
#define CHECK(p)                                               \
  do {                                                         \
    if (!(p)) {                                                \
      tket_log()->warn("Invalid DAG: check (" #p ") failed."); \
      return false;                                            \
    }                                                          \
  } while (0)

enum class VertexType { Quantum, Classical, Measure };

bool is_valid(const DAG &G) {
  BGL_FORALL_VERTICES(v, G, DAG) {
    EdgeSet q_in, c_in, b_in;
    BGL_FORALL_INEDGES(v, e, G, DAG) {
      switch (G[e].type) {
        case EdgeType::Quantum:
          q_in.insert(e);
          break;
        case EdgeType::Classical:
          c_in.insert(e);
          break;
        case EdgeType::Boolean:
          b_in.insert(e);
          break;
        default:
          CHECK(!"unknown edge type");
      }
    }
    EdgeSet q_out, c_out, b_out;
    BGL_FORALL_OUTEDGES(v, e, G, DAG) {
      switch (G[e].type) {
        case EdgeType::Quantum:
          q_out.insert(e);
          break;
        case EdgeType::Classical:
          c_out.insert(e);
          break;
        case EdgeType::Boolean:
          b_out.insert(e);
          break;
        default:
          CHECK(!"unknown edge type");
      }
    }
    std::set<port_t> in_ports, q_in_ports, q_out_ports, c_in_ports, c_out_ports,
        b_in_ports;
    for (const auto &e : q_in) {
      port_t p = G[e].ports.second;
      in_ports.insert(p);
      q_in_ports.insert(p);
    }
    for (const auto &e : q_out) {
      q_out_ports.insert(G[e].ports.first);
    }
    for (const auto &e : c_in) {
      port_t p = G[e].ports.second;
      in_ports.insert(p);
      c_in_ports.insert(p);
    }
    for (const auto &e : c_out) {
      c_out_ports.insert(G[e].ports.first);
    }
    for (const auto &e : b_in) {
      port_t p = G[e].ports.second;
      in_ports.insert(p);
      b_in_ports.insert(p);
    }

    // Now check the required properties.

    CHECK(
        in_ports.size() ==
        q_in_ports.size() + c_in_ports.size() + b_in_ports.size());

    // Every Boolean out port matches a Classical out port.
    for (const Edge &e : b_out) {
      port_t p = G[e].ports.first;
      CHECK(std::any_of(c_out.cbegin(), c_out.cend(), [&](const Edge &f) {
        return G[f].ports.first == p;
      }));
    }

    if (c_in.empty() && c_out.empty()) {
      // It is a Quantum vertex.
      unsigned in_deg = q_in.size();
      unsigned out_deg = q_out.size();
      CHECK(q_in_ports.size() == in_deg);
      CHECK(q_out_ports.size() == out_deg);
      CHECK(
          (in_deg == 0 && out_deg == 1) || (in_deg == 1 && out_deg == 0) ||
          (q_in_ports == q_out_ports));  // bijection between in and out ports
      CHECK(b_out.empty());
    } else if (q_in.empty() && q_out.empty()) {
      // It is a Classical vertex.
      unsigned in_deg = c_in.size();
      unsigned out_deg = c_out.size();
      CHECK(c_in_ports.size() == in_deg);
      CHECK(c_out_ports.size() == out_deg);
      CHECK(
          (in_deg == 0 && out_deg == 1) || (in_deg == 1 && out_deg == 0) ||
          (c_in_ports == c_out_ports));  // bijection between in and out ports
    } else {
      // Check that it is a Measure vertex.
      CHECK(
          q_in.size() == 1 && q_out.size() == 1 && c_in.size() == 1 &&
          c_out.size() == 1);
      CHECK(q_in_ports == q_out_ports && c_in_ports == c_out_ports);
    }
  }
  return true;
}

#undef CHECK

}  // namespace tket
