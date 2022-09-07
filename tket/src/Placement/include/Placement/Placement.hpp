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

#include "Architecture/Architecture.hpp"

namespace tket {

class Placement {
 public: 

  typedef std:::shared_ptr<Placement> Ptr;

    explicit Placement(const Architecture& _arc) : arc_(_arc) {}
    
    /**
     * Reassigns some UnitID in circ_ as UnitID arc_
     * 
     * @param circ_ Circuit to be relabelled
     * 
     */
    bool place(Circuit& circ_) const;

    /**
     * 
     * For some Circuit, returns a map between Circuit UnitID and
     * Architecture UnitID that can be used for reassigning UnitID in
     * Circuit. Map is expected to give best performance for given method.
     * 
     * @param circ_ Circuit relabelling map is constructed from
     * 
     * @return Map between Circuit and Architecture UnitID
     */
    std::map<UnitID, UnitID> get_placement_map(const Circuit& circ_) const;


    /**
     * 
     * For some Circuit, returns maps between Circuit UnitID and
     * Architecture UnitID that can be used for reassigning UnitID in
     * Circuit. Maps expected to give similiar performance for given method.
     * 
     * @param circ_ Circuit relabelling map is constructed from
     * 
     * @return Map between Circuit and Architecture UnitID
     */
    virtual std::vector<std::map<UnitID, UnitID>> get_all_placement_maps(
      const Circuit& circ_) const;

 protected:
    Architecture arc_;

};

JSON_DECL(Placement::Ptr);



class GraphPlacement {
 public:
    /**
     * Holds information for constructing a weighted edge in a QubitGraph.
     * @param node0 UnitID for first node in edge
     * @param node1 UnitID for second node in edge
     * @param weight Unsigned giving a weight for implied edge
     */
    struct WeightedEdge {
        UnitID node0;
        UnitID node1;
        unsigned weight;
    };

    static const std::vector<WeightedEdge> default_weighting(const Circuit& circuit);

    explicit GraphPlacement(const Architecture& _arc, unsigned _maximum_matches, unsigned _timeout, const std::function<std::vector<WeightedEdge>(const Circuit&)> _weight_pattern_graph) : arc_(_arc), maximum_matches_(_maximum_matches), timeout_(_timeout), weight_pattern_graph_(_weight_pattern_graph) {};

    /**
     * 
     * For some Circuit, returns maps between Circuit UnitID and
     * Architecture UnitID that can be used for reassigning UnitID in
     * Circuit. Maps are constructed by running a Weighted Subgraph Monomorphism
     * for the given problem and returning maximum_matches_ number of
     * potential solutions, ranked.
     * 
     * @param circ_ Circuit relabelling map is constructed from
     * 
     * @return Map between Circuit and Architecture UnitID
     */
    std::vector<std::map<UnitID, UnitID>> get_all_placement_maps(
      const Circuit& circ_) const override;



 protected:

    const std::function<std::vector<WeightedEdge>(const Circuit&)> weight_pattern_graph_;
    unsigned maximum_matches_;
    unsigned timeout_;

    QubitGraph::UndirectedConnGraph construct_pattern_graph(const std::vector<WeightedEdge>& edges) const;

};


} // namespace tket