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

#include <catch2/catch_test_macros.hpp>
#include <sstream>
#include <tkrng/RNG.hpp>
#include <tkwsm/Searching/DomainsAccessor.hpp>
#include <tkwsm/Searching/NodeListTraversal.hpp>
#include <tkwsm/Searching/NodesRawData.hpp>

namespace tket {
namespace WeightedSubgraphMonomorphism {

// It's hard to test solving problems and debugging,
// because search trees can be so large.
// This provides another way to test directly.
SCENARIO("Test detailed complete search") {
  DomainInitialiser::InitialDomains initial_domains(3);
  const unsigned percentage_for_tv_in_domain = 30;
  RNG rng;

  for (unsigned pv = 0; pv < initial_domains.size(); ++pv) {
    initial_domains[pv].resize(4);
    for (unsigned tv = 0; tv < 3; ++tv) {
      initial_domains[pv].set(pv);
      initial_domains[pv].set(pv + 1);
      if (rng.check_percentage(percentage_for_tv_in_domain)) {
        initial_domains[pv].set(tv);
      }
    }
  }
  NodesRawDataWrapper wrapper(initial_domains, 4);
  DomainsAccessor accessor(wrapper);
  NodeListTraversal traversal(wrapper);

  // May change elsewhere, but the reference remains valid.
  const NodesRawData& raw_data = wrapper.get_raw_data_for_debug();

  std::stringstream ss;

  for (int count = 0; count < 100; ++count) {
    if (count != 0) {
      // Assign some PV; just choose the first.
      unsigned pv_to_assign = 0;
      bool move_up = true;
      for (; pv_to_assign < initial_domains.size(); ++pv_to_assign) {
        const auto size = accessor.get_domain_size(pv_to_assign);
        if (size == 0) {
          break;
        }
        if (size > 1) {
          move_up = false;
          break;
        }
      }
      ss << "\n" << raw_data.current_node_index() << ".";
      if (move_up) {
        ss << "u";
        if (!traversal.move_up()) {
          ss << "_FIN";
          break;
        }
        continue;
      }
      const auto tv = accessor.get_domain(pv_to_assign).find_first();

      ss << "d" << pv_to_assign << tv;
      accessor.clear_new_assignments();
      traversal.move_down(pv_to_assign, tv);
      if (!accessor.alldiff_reduce_current_node(0)) {
        ss << "_ng";
      }
    }
    ss << "_nn" << raw_data.nodes_data.size() << "_{";
    for (unsigned pv = 0; pv < raw_data.domains_data.size(); ++pv) {
      ss << pv << ":";
      const auto& domain = accessor.get_domain(pv);
      for (auto tv = domain.find_first(); tv < domain.size();
           tv = domain.find_next(tv)) {
        ss << tv;
      }
      const auto size = raw_data.domains_data[pv].entries.size();
      REQUIRE(size > 0);
      const auto& last_entry = raw_data.domains_data[pv].entries[size - 1];
      REQUIRE(domain == last_entry.domain);
      REQUIRE(last_entry.node_index <= raw_data.current_node_index());
      ss << "(s" << size << ",ni" << last_entry.node_index << ") ";
    }
    ss << "}";
  }
  // clang-format off
  CHECK(ss.str() == "_nn1_{0:01(s1,ni0) 1:12(s1,ni0) 2:0123(s1,ni0) }"
    "\n0.d00_nn2_{0:0(s2,ni1) 1:12(s1,ni0) 2:123(s2,ni1) }"
    "\n1.d11_nn3_{0:0(s2,ni1) 1:1(s3,ni2) 2:23(s3,ni2) }"
    "\n2.d22_nn4_{0:0(s2,ni1) 1:1(s3,ni2) 2:2(s4,ni3) }"
    "\n3.u"
    "\n2.u"
    "\n1.d21_nn3_{0:0(s2,ni1) 1:2(s2,ni1) 2:1(s3,ni2) }"
    "\n2.u"
    "\n1.d22_ng_nn3_{0:0(s2,ni1) 1:2(s2,ni1) 2:2(s3,ni2) }"
    "\n2.u"
    "\n1.u"
    "\n0.d11_ng_nn2_{0:1(s1,ni0) 1:1(s2,ni1) 2:0123(s1,ni0) }"
    "\n1.d20_nn3_{0:1(s1,ni0) 1:1(s2,ni1) 2:0(s3,ni2) }"
    "\n2.u"
    "\n1.d21_ng_nn3_{0:1(s1,ni0) 1:1(s2,ni1) 2:1(s3,ni2) }"
    "\n2.u"
    "\n1.d22_nn3_{0:1(s1,ni0) 1:1(s2,ni1) 2:2(s3,ni2) }"
    "\n2.u"
    "\n1.u"
    "\n0.d20_nn2_{0:1(s1,ni0) 1:2(s1,ni0) 2:0(s2,ni1) }"
    "\n1.u"
    "\n0.d21_ng_nn2_{0:1(s1,ni0) 1:2(s1,ni0) 2:1(s2,ni1) }"
    "\n1.u"
    "\n0.d22_ng_nn2_{0:1(s1,ni0) 1:2(s1,ni0) 2:2(s2,ni1) }"
    "\n1.u"
    "\n0.u_FIN");
  // clang-format on
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
