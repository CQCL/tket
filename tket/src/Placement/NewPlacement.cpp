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

#include "Placement.hpp"


namespace tket {

bool Placement::place(Circuit& circ_) const {
    std::map<UnitID, UnitID> map_ = this->get_placement_map(circ_);
    std::map<Qubit, Node> recast_map;
    for(const std::pair<UnitID, UnitID>& entry : map_){
        recast_map.insert({Qubit(entry.first), Node(entry.second)});
    }
    return circ_.rename_units(recast_map);
}


std::map<UnitID, UnitID> Placement::get_placement_map(const Circuit& circ_) const{
    return this->get_all_placement_maps(circ_)[0];
}   

std::vector<std::map<UnitID, UnitID> Placement::get_all_placement_maps(const Circuit& circ_) const {
    return {};
}

QubitGraph::UndirectedConnGraph construct_pattern_graph(const std::vector<WeightedEdge>& edges) const {
    QubitGraph q_graph();
    for(const WeightedEdge& weighted_edge : edges){
        Node node0 = Node(weighted_edge.node0);
        Node node1 = Node(weighted_edge.node1);
        bool e_01_exists = q_graph.edge_exists(node0, node1);
        bool e_10_exists = q_graph.edge_exists(node1, node0);
        if(e_01_exists || e_10_exists){
            throw std::invalid_argument("Graph can only have on edge between a pair of Node.");
        }
        if(!node_exists(node0)){
            q_graph.add_node(node0);
        }
        if(!node_exists(node1)){
            q_graph.add_node(node1);
        }
        q_graph.add_connection(node0, node1, weighted_edge.weight);
    }
    return q_graph;
}


std::vector<std::map<UnitID, UnitID>> GraphPlacement::get_all_placement_maps(const Circuit& circ_) const {
    std::vector<WeightedEdge> weighted_edges = this->weight_pattern_graph_(circ_);
    QubitGraph::UndirectedConnGraph pattern_graph = this->construct_pattern_graph(weighted_edges);
    Architecture::UndirectedConnGraph target_graph = this->arc_.get_undirected_connectivity();
    std::vector<qubit_bimap_t> all_maps = get_weighted_subgraph_monomorphism(pattern_graph, target_graph, this->maximum_matches_, this->timeout_);
    return all_maps;
}

    
}