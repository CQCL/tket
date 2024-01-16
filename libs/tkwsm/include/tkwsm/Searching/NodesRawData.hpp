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

#pragma once
#include <boost/dynamic_bitset.hpp>
#include <string>

#include "../Common/LogicalStack.hpp"
#include "../GraphTheoretic/DomainInitialiser.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

// The "logical" data structure for a search tree would be a std::vector
// of "node" objects, with each node containing an object like
// std::map<VertexWSM, std::set<VertexWSM>> (with key being the PV,
// and value being Domain(pv)).
//
// Think of this as a 2D structure,
// with x-axis being the PV, and y-axis being the node index
// (i.e., the index in the std::vector of "node" objects).
//
// (Thus, logically, we have a single vertical strand - the nodes vector -
// with each node having a horizontal bar - the domains).
//
// However, for better performance we reduce copying
// by only storing Domain(pv) when it CHANGES.
//
// Thus, the actual data structure is quite different; each PV
// has its own std::vector of DOMAINS. Thus, we have a whole collection
// of "vertical strands" (columns) - one for each PV.
//
// We need two kinds of movement along the data: "horizontal"
// (i.e., for the current logical node, get information about the domains),
// provided by the class DomainsAccessor;
// and "vertical" (i.e., as we search, we move up and down the std::vector
// of node objects), provided by the class NodeListTraversal.
//
// Thus, this raw data is SHARED between these two classes, each of which
// understands how to do its own kind of motion.

/** This is the data structure containing all the data about nodes
 * and domains necessary for searching, but without many manipulation
 * functions; it's difficult to split the data up further, so instead
 * we divide up the data manipulation tasks between different classes.
 */
struct NodesRawData {
  // Think of DomainData as being "vertical"; each vertical strand
  // is information about one particular PV.
  // NodeData is "horizontal": each one gives all additional information
  // about one particular Node in a list,
  // which cuts roughly horizontally across the DomainData strands.

  /** DomainData is "vertical"; each particular object is a "vertical strand"
   * giving information about one particular PV (all along the search nodes).
   */
  struct DomainData {
    struct Entry {
      boost::dynamic_bitset<> domain;

      // The index in the std::vector of node objects
      // when this new domain was first set (i.e., changed from
      // the domain in the previous node).
      unsigned node_index;
    };

    // Not all domains are stored; only those which have just changed.
    // This is for performance reasons; otherwise we would have to copy
    // the domain for EVERY pv, at EVERY node, even if only a few domains
    // changed from one node to the next.
    LogicalStack<Entry> entries;

    std::string str() const;
  };

  /** "vertical strands": element[pv] gives information about all Domain(pv).
   */
  std::vector<DomainData> domains_data;

  /** This is the "horizontal" data; each one has information about the node
   * which does NOT depend on the PV.
   */
  struct NodeData {
    /** If true, this node has been invalidated (is impossible). */
    bool nogood;

    /** The total weight (scalar product) over all edges currently assigned. */
    WeightWSM scalar_product;

    /** The sum of the pattern edge weights which have currently been assigned.
     * (To decide between incomplete solutions, we may care about which one has
     * most assigned total weight).
     */
    WeightWSM total_p_edge_weights;

    /** Consists of PV->TV assignments which only occurred in this node;
     * some may be deleted again after processing.
     */
    std::vector<std::pair<VertexWSM, VertexWSM>> new_assignments;

    /** For performance: rather than searching through EVERY Dom(PV),
     * which will include many assigned PV (for which Dom(PV)={y}),
     * we have a list of all those PV known to have |Dom(PV)| > 1.
     * However, we don't need a strictly accurate list.
     * This is EITHER empty, OR guarantees that every PV with
     * size Dom(PV) > 1 is included in here.
     * However, some extra PV with  |Dom(PV)| = 1  may also be included.
     *
     * Note that this only works because, as we move down the search tree,
     * the domain sets reduce (so we can restrict to a previous list of PVs).
     */
    std::vector<VertexWSM> unassigned_vertices_superset;

    std::string str() const;
  };

  LogicalStack<NodeData> nodes_data;

  explicit NodesRawData(
      const DomainInitialiser::InitialDomains& initial_domains,
      std::size_t number_of_tv);

  const std::size_t number_of_tv;

  // The TV in "initial_domains" passed into the constructor
  // already have been relabelled to be {0,1,2,...,M} for some M.
  // But some might be unused. We COULD relabel AGAIN, but completely unused
  // TV are quite rare, so it's not worth the extra complication.
  // HOWEVER, the "number_of_tv" is set simply by looking at the max TV
  // occurring; so it's possible that some unused TV, e.g. M, from the end of
  // {0,1,2,...,M} might get removed and not included in this.
  // const std::size_t number_of_tv;

  const NodeData& get_current_node() const;
  NodeData& get_current_node_nonconst();

  unsigned current_node_index() const;

  DomainData& get_most_recent_domain_data_for_pv(VertexWSM);
};

/** We are careful to restrict access to the raw data. */
class NodesRawDataWrapper {
 public:
  explicit NodesRawDataWrapper(
      const DomainInitialiser::InitialDomains& initial_domains,
      std::size_t number_of_tv);

  /** For debugging and testing, it's helpful to access the raw data.
   * @return A const reference to the internal NodesRawData object.
   */
  const NodesRawData& get_raw_data_for_debug() const;

 private:
  NodesRawData m_raw_data;
  friend class DomainsAccessor;
  friend class NodeListTraversal;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
