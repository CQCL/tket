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
#include <string>

#include "../GraphTheoretic/GeneralStructs.hpp"

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
 * functions.
 */
struct NodesRawData {
  // Think of DomainData as being "vertical"; each vertical strand
  // is information about one particular PV.
  // NodeData is "horizontal": each one gives all additional information
  // about one particular Node in a list,
  // which cuts roughly horizontally across the DomainData strands.

  // A sorted list of all PV.
  const std::vector<VertexWSM> pattern_vertices;

  /** DomainData is "vertical"; each particular object is a "vertical strand"
   * giving information about one particular PV.
   */
  struct DomainData {
    struct Entry {
      std::set<VertexWSM> domain;

      // The index in the std::vector of node objects.
      unsigned node_level;
    };

    // Not all domains are stored; only those which have just changed.
    std::vector<Entry> entries;

    // The index of the "logical back()" of "entries". All data beyond
    // is junk to be reused (faster than reallocating).
    unsigned entries_back_index;

    std::string str() const;
  };

  /** Key: pv  Value: A "vertical strand": information about all Domain(pv).
   */
  std::map<VertexWSM, DomainData> domains_data;

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

    /** When searching, we first look at vertices adjacent
     * (in the pattern graph) to a vertex already assigned.
     * TODO: is this actually a good idea? Seems natural, but we should test!
     */
    std::set<VertexWSM> pvs_adjacent_to_newly_assigned_vertices;

    /** For performance: rather than searching through EVERY Dom(PV),
     * which will include many assigned PV (for which Dom(PV)={y}),
     * we have a set giving all those PV with |Dom(PV)| > 1.
     * However, it would be too slow to maintain an accurate list,
     * so instead this is EITHER an empty set,
     * OR every PV with size Dom(PV) > 1 is included in here.
     * However, some extra PV with  |Dom(PV)| = 1  may also be included.
     */
    std::set<VertexWSM> unassigned_vertices_superset;

    std::string str() const;
  };

  std::vector<NodeData> nodes_data;

  /** The index of the current search node in "nodes_data".
   * Always valid; but any NodeData objects beyond this
   * are junk data to be reused (much faster than resizing and
   * reallocating as we move vertically).
   */
  unsigned current_node_level;

  explicit NodesRawData(const PossibleAssignments& possible_assignments);

  const NodeData& get_current_node() const;
  NodeData& get_current_node_nonconst();

  DomainData& get_most_recent_domain_data_for_pv(VertexWSM);
};

/** We are careful to restrict access to the raw data. */
class NodesRawDataWrapper {
 public:
  explicit NodesRawDataWrapper(const PossibleAssignments& possible_assignments);

 private:
  NodesRawData m_raw_data;
  friend class DomainsAccessor;
  friend class NodeListTraversal;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
