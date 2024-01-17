// Copyright 2019-2024 Cambridge Quantum Computing
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
#include <tklog/TketLog.hpp>

#include "tket/Circuit/Circuit.hpp"
#include "tket/OpType/EdgeType.hpp"
#include "tket/Utils/GraphHeaders.hpp"

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
    EdgeSet q_in, c_in, b_in, w_in;
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
        case EdgeType::WASM:
          w_in.insert(e);
          break;
        default:
          CHECK(!"found unknown edge type in is_valid check");
      }
    }
    EdgeSet q_out, c_out, b_out, w_out;
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
        case EdgeType::WASM:
          w_out.insert(e);
          break;
        default:
          CHECK(!"found unknown edge type in is_valid check");
      }
    }
    std::set<port_t> in_ports, q_in_ports, q_out_ports, c_in_ports, c_out_ports,
        b_in_ports, w_in_ports, w_out_ports;
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
    for (const auto &e : w_in) {
      port_t p = G[e].ports.second;
      in_ports.insert(p);
      w_in_ports.insert(p);
    }
    for (const auto &e : w_out) {
      w_out_ports.insert(G[e].ports.first);
    }

    // Now check the required properties.

    CHECK(
        in_ports.size() == q_in_ports.size() + c_in_ports.size() +
                               b_in_ports.size() + w_in_ports.size());

    // Every Boolean out port matches a Classical out port.
    for (const Edge &e : b_out) {
      port_t p = G[e].ports.first;
      CHECK(std::any_of(c_out.cbegin(), c_out.cend(), [&](const Edge &f) {
        return G[f].ports.first == p;
      }));
    }

    // get size and empty only once

    unsigned q_in_deg = q_in.size();
    unsigned q_out_deg = q_out.size();

    unsigned c_in_deg = c_in.size();
    unsigned c_out_deg = c_out.size();

    unsigned b_in_deg = b_in.size();

    unsigned w_in_deg = w_in.size();
    unsigned w_out_deg = w_out.size();

    bool q_in_empty = q_in.empty();
    bool q_out_empty = q_out.empty();

    bool c_in_empty = c_in.empty();
    bool c_out_empty = c_out.empty();

    bool b_in_empty = b_in.empty();
    bool b_out_empty = b_out.empty();

    bool w_in_empty = w_in.empty();
    bool w_out_empty = w_out.empty();

    // check that ports maps to edges

    CHECK(q_in_ports.size() == q_in_deg);
    CHECK(q_out_ports.size() == q_out_deg);

    CHECK(w_in_ports.size() == w_in_deg);
    CHECK(w_out_ports.size() == w_out_deg);

    CHECK(c_in_ports.size() == c_in_deg);
    CHECK(c_out_ports.size() == c_out_deg);

    CHECK(b_in_ports.size() == b_in_deg);

    if (c_in_empty && c_out_empty && q_in_empty && !q_out_empty) {
      // quantum in vertex
      CHECK(b_in_empty);
      CHECK(b_out_empty);
      CHECK(w_in_empty);
      CHECK(w_out_empty);
      CHECK(q_in_deg == 0 && q_out_deg == 1);
    } else if (c_in_empty && c_out_empty && !q_in_empty && q_out_empty) {
      // quantum out vertex
      CHECK(b_in_empty);
      CHECK(b_out_empty);
      CHECK(w_in_empty);
      CHECK(w_out_empty);
      CHECK(q_in_deg == 1 && q_out_deg == 0);
    } else if (
        c_in_empty && c_out_empty && !q_in_empty && !q_out_empty &&
        !b_in_empty) {
      // conditional quantum vertex, can have bool edges
      CHECK(b_out_empty);
      CHECK(w_in_empty);
      CHECK(w_out_empty);
      CHECK(q_in_ports == q_out_ports);  // bijection between in and out ports
    } else if (
        c_in_empty && c_out_empty && !q_in_empty && !q_out_empty &&
        b_in_empty) {
      // quantum vertex
      CHECK(b_out_empty);
      CHECK(w_in_empty);
      CHECK(w_out_empty);
      CHECK(q_in_ports == q_out_ports);  // bijection between in and out ports
    } else if (w_in_empty && !w_out_empty) {
      // wasm out vertex
      CHECK(b_in_empty);
      CHECK(b_out_empty);
      CHECK(c_in_empty);
      CHECK(c_out_empty);
      CHECK(q_in_empty);
      CHECK(q_out_empty);
      CHECK(w_in_deg == 0 && w_out_deg == 1);
    } else if (!w_in_empty && w_out_empty) {
      // wasm in vertex
      CHECK(b_in_empty);
      CHECK(b_out_empty);
      CHECK(c_in_empty);
      CHECK(c_out_empty);
      CHECK(q_in_empty);
      CHECK(q_out_empty);
      CHECK(w_in_deg == 1 && w_out_deg == 0);
    } else if (!w_in_empty && !w_out_empty) {
      // wasm vertex
      CHECK(q_in_empty);
      CHECK(q_out_empty);
      CHECK(c_in_ports == c_out_ports);  // bijection between in and out ports
      CHECK(w_in_ports == w_out_ports);  // bijection between in and out ports
    } else if (c_in_empty && !c_out_empty) {
      // classical in vertex, can have bool out edges
      CHECK(b_in_empty);
      CHECK(w_in_empty);
      CHECK(w_out_empty);
      CHECK(q_in_empty);
      CHECK(q_out_empty);
      CHECK(c_in_deg == 0 && c_out_deg == 1);
    } else if (!c_in_empty && c_out_empty) {
      // classical out vertex, can't have bool out edges
      CHECK(b_out_empty);
      CHECK(w_in_empty);
      CHECK(w_out_empty);
      CHECK(q_in_empty);
      CHECK(q_out_empty);
      CHECK(c_in_deg == 1 && c_out_deg == 0);
    } else if (q_in_empty && q_out_empty && !c_in_empty && !c_out_empty) {
      // classical vertex, can have bool in and out edges
      CHECK(w_in_empty);
      CHECK(w_out_empty);
      CHECK(c_in_ports == c_out_ports);  // bijection between in and out ports
    } else if (!q_in_empty && !q_out_empty && !c_in_empty && !c_out_empty) {
      // mixed vertex, can have bool in and out edges
      CHECK(w_in_empty);
      CHECK(w_out_empty);
      CHECK(c_in_ports == c_out_ports);  // bijection between in and out ports
      CHECK(q_in_ports == q_out_ports);  // bijection between in and out ports
    } else if (q_in_empty && q_out_empty && c_in_empty && c_out_empty) {
      // unconected vertex, can't have bool in and out edges
      CHECK(w_in_empty);
      CHECK(w_out_empty);
      CHECK(b_in_empty);
      CHECK(b_out_empty);
    } else {
      TKET_ASSERT(!"Found invalid circuit");
    }
  }
  return true;
}

#undef CHECK

}  // namespace tket
