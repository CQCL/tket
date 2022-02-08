#include "Mapping/LexiRoute.hpp"

#include "Mapping/MappingFrontier.hpp"
#include "Utils/Json.hpp"

namespace tket {

LexiRoute::LexiRoute(
    const ArchitecturePtr& _architecture,
    std::shared_ptr<MappingFrontier>& _mapping_frontier)
    : architecture_(_architecture), mapping_frontier_(_mapping_frontier) {
  this->set_interacting_uids();
  // set initial logical->physical labelling
  for (const Qubit& qb : this->mapping_frontier_->circuit_.all_qubits()) {
    this->labelling_.insert({qb, qb});
    Node n(qb);
    // store which Node have been asigned to Circuit already
    if (this->architecture_->node_exists(n)) {
      this->assigned_nodes_.insert(n);
    }
  }
}

bool LexiRoute::assign_at_distance(
    const UnitID& assignee, const Node& root, unsigned distances) {
  node_set_t valid_nodes;
  for (const Node& neighbour :
       this->architecture_->nodes_at_distance(root, distances)) {
    if (this->assigned_nodes_.find(neighbour) == this->assigned_nodes_.end() ||
        this->mapping_frontier_->ancilla_nodes_.find(neighbour) !=
            this->mapping_frontier_->ancilla_nodes_.end()) {
      valid_nodes.insert(neighbour);
    }
  }
  if (valid_nodes.size() == 1) {
    auto it = valid_nodes.begin();
    if (this->mapping_frontier_->ancilla_nodes_.find(*it) !=
        this->mapping_frontier_->ancilla_nodes_.end()) {
      // => node *it is already present in circuit, but as an ancilla
      this->mapping_frontier_->merge_ancilla(assignee, *it);
      this->mapping_frontier_->ancilla_nodes_.erase(*it);
      this->labelling_.erase(*it);
      this->labelling_[assignee] = *it;
    } else {
      this->labelling_[assignee] = *it;
      this->assigned_nodes_.insert(*it);
    }
    return true;
  }
  if (valid_nodes.size() > 1) {
    auto it = valid_nodes.begin();
    lexicographical_distances_t winning_distances =
        this->architecture_->get_distances(*it);
    Node preserved_node = *it;
    ++it;
    for (; it != valid_nodes.end(); ++it) {
      lexicographical_distances_t comparison_distances =
          this->architecture_->get_distances(*it);
      if (comparison_distances < winning_distances) {
        preserved_node = *it;
        winning_distances = comparison_distances;
      }
    }
    if (this->mapping_frontier_->ancilla_nodes_.find(preserved_node) !=
        this->mapping_frontier_->ancilla_nodes_.end()) {
      // => node *it is already present in circuit, but as an ancilla
      this->mapping_frontier_->merge_ancilla(assignee, preserved_node);
      this->mapping_frontier_->ancilla_nodes_.erase(preserved_node);
      this->labelling_.erase(preserved_node);
      this->labelling_[assignee] = preserved_node;
    } else {
      // add ancilla case
      this->labelling_[assignee] = preserved_node;
      this->assigned_nodes_.insert(preserved_node);
    }
    return true;
  }
  return false;
}

bool LexiRoute::update_labelling() {
  // iterate through interacting qubits, assigning them to an Architecture Node
  // if they aren't already
  bool relabelled = false;
  for (const auto& pair : this->interacting_uids_) {
    bool uid_0_exist =
        this->architecture_->node_exists(Node(this->labelling_[pair.first]));
    bool uid_1_exist =
        this->architecture_->node_exists(Node(this->labelling_[pair.second]));
    if (!uid_0_exist || !uid_1_exist) {
      relabelled = true;
    }
    if (!uid_0_exist && !uid_1_exist) {
      // Place one on free unassigned qubit
      // Then place second later
      // condition => No ancilla qubits assigned, so don't checl
      if (this->assigned_nodes_.size() == 0) {
        // find nodes with best averaged distance to other nodes
        // place it there...
        std::set<Node> max_degree_nodes =
            this->architecture_->max_degree_nodes();
        auto it = max_degree_nodes.begin();
        lexicographical_distances_t winning_distances =
            this->architecture_->get_distances(*it);
        Node preserved_node = Node(*it);
        ++it;
        for (; it != max_degree_nodes.end(); ++it) {
          lexicographical_distances_t comparison_distances =
              this->architecture_->get_distances(*it);
          if (comparison_distances < winning_distances) {
            preserved_node = Node(*it);
            winning_distances = comparison_distances;
          }
        }
        this->labelling_[pair.first] = preserved_node;
        this->assigned_nodes_.insert(preserved_node);
        uid_0_exist = true;
        // given best node, do something
      } else {
        auto root_it = this->assigned_nodes_.begin();
        while (!uid_0_exist && root_it != this->assigned_nodes_.end()) {
          Node root = *root_it;
          uid_0_exist = this->assign_at_distance(pair.first, root, 1);
          ++root_it;
        }
        if (!uid_0_exist) {
          throw LexiRouteError(
              "Unable to assign physical qubit - no free qubits remaining.");
        }
      }
    }
    if (!uid_0_exist && uid_1_exist) {
      Node root(this->labelling_[pair.second]);
      for (unsigned k = 1; k <= this->architecture_->get_diameter(); k++) {
        uid_0_exist = this->assign_at_distance(pair.first, root, k);
        if (uid_0_exist) {
          break;
        }
      }
      if (!uid_0_exist) {
        throw LexiRouteError(
            "Unable to assign physical qubit - no free qubits remaining.");
      }
    }
    if (uid_0_exist && !uid_1_exist) {
      Node root(this->labelling_[pair.first]);
      for (unsigned k = 1; k <= this->architecture_->get_diameter(); k++) {
        uid_1_exist = this->assign_at_distance(pair.second, root, k);
        if (uid_1_exist) {
          break;
        }
      }
      if (!uid_1_exist) {
        throw LexiRouteError(
            "Unable to assign physical qubit - no free qubits remaining.");
      }
    }
  }
  return relabelled;
}

/**
 * LexiRoute::set_interacting_uids
 * Updates this->interacting_uids_ with all "interacting" pairs
 * of UnitID in this->mapping_frontier_
 */
void LexiRoute::set_interacting_uids(bool assigned_only) {
  // return types
  this->interacting_uids_.clear();
  for (auto it =
           this->mapping_frontier_->quantum_boundary->get<TagKey>().begin();
       it != this->mapping_frontier_->quantum_boundary->get<TagKey>().end();
       ++it) {
    Edge e0 = this->mapping_frontier_->circuit_.get_nth_out_edge(
        it->second.first, it->second.second);
    Vertex v0 = this->mapping_frontier_->circuit_.target(e0);
    // should never be input vertex, so can always use in_edges
    if (this->mapping_frontier_->circuit_.get_OpType_from_Vertex(v0) !=
        OpType::Barrier) {
      int n_edges = this->mapping_frontier_->circuit_.n_in_edges_of_type(
          v0, EdgeType::Quantum);
      // make forwards = backwards
      if (n_edges == 2) {
        auto jt = it;
        ++jt;
        for (; jt !=
               this->mapping_frontier_->quantum_boundary->get<TagKey>().end();
             ++jt) {
          // i.e. if vertices match
          Edge e1 = this->mapping_frontier_->circuit_.get_nth_out_edge(
              jt->second.first, jt->second.second);
          Vertex v1 = this->mapping_frontier_->circuit_.target(e1);
          if (v0 == v1) {
            // we can assume a qubit will only be in one interaction
            // we can assume from how we iterate through pairs that each qubit
            // will only be found in one match
            if (!assigned_only ||
                (this->architecture_->node_exists(Node(it->first)) &&
                 this->architecture_->node_exists(Node(jt->first)))) {
              interacting_uids_.insert({it->first, jt->first});
              interacting_uids_.insert({jt->first, it->first});
            }
          }
        }
      } else if (
          n_edges > 2 &&
          this->mapping_frontier_->circuit_.get_OpType_from_Vertex(v0) !=
              OpType::Barrier) {
        throw LexiRouteError(
            "LexiRoute only supports non-Barrier vertices with 1 or 2 edges.");
      }
    }
  }
}

swap_set_t LexiRoute::get_candidate_swaps() {
  swap_set_t candidate_swaps;
  for (const auto& interaction : this->interacting_uids_) {
    Node assigned_first = Node(this->labelling_[interaction.first]);
    std::vector<Node> adjacent_uids_0 =
        this->architecture_->nodes_at_distance(assigned_first, 1);
    TKET_ASSERT(adjacent_uids_0.size() != 0);
    for (const Node& neighbour : adjacent_uids_0) {
      if (candidate_swaps.find({neighbour, assigned_first}) ==
          candidate_swaps.end()) {
        candidate_swaps.insert({assigned_first, neighbour});
      }
    }
    Node assigned_second = Node(this->labelling_[interaction.second]);
    std::vector<Node> adjacent_uids_1 =
        this->architecture_->nodes_at_distance(assigned_second, 1);
    TKET_ASSERT(adjacent_uids_1.size() != 0);
    for (const Node& neighbour : adjacent_uids_1) {
      if (candidate_swaps.find({neighbour, assigned_second}) ==
          candidate_swaps.end()) {
        candidate_swaps.insert({assigned_second, neighbour});
      }
    }
  }
  return candidate_swaps;
}

bool is_vertex_CX(const Circuit& circ_, const Vertex& v) {
  OpType ot = circ_.get_OpType_from_Vertex(v);
  if (ot != OpType::CX) {
    if (ot == OpType::Conditional) {
      const Conditional& b =
          static_cast<const Conditional&>(*circ_.get_Op_ptr_from_Vertex(v));
      if (b.get_op()->get_type() != OpType::CX) {
        return false;
      }
    } else {
      return false;
    }
  }
  return true;
}

std::pair<bool, bool> LexiRoute::check_bridge(
    const std::pair<Node, Node>& swap, unsigned lookahead) {
  std::pair<bool, bool> output = {false, false};
  // first confirm whether it even has an interaction
  auto it = this->interacting_uids_.find(swap.first);
  if (it != this->interacting_uids_.end()) {  // => in interaction
    if (this->architecture_->get_distance(swap.first, Node(it->second)) ==
        2) {  // => could be bridge
      // below should always return correct object given prior checks
      VertPort vp =
          (*this->mapping_frontier_->quantum_boundary->find(swap.first)).second;
      Edge out_edge = this->mapping_frontier_->circuit_.get_nth_out_edge(
          vp.first, vp.second);
      output.first = is_vertex_CX(
          this->mapping_frontier_->circuit_,
          this->mapping_frontier_->circuit_.target(out_edge));
    }
  }
  // repeat for second swap
  it = this->interacting_uids_.find(swap.second);
  if (it != this->interacting_uids_.end()) {
    if (this->architecture_->get_distance(swap.second, Node(it->second)) == 2) {
      VertPort vp =
          (*this->mapping_frontier_->quantum_boundary->find(swap.second))
              .second;
      Edge out_edge = this->mapping_frontier_->circuit_.get_nth_out_edge(
          vp.first, vp.second);
      output.second = is_vertex_CX(
          this->mapping_frontier_->circuit_,
          this->mapping_frontier_->circuit_.target(out_edge));
    }
  }
  if ((output.first && output.second) || (!output.first && !output.second)) {
    return {0, 0};
  }
  // implies conditions are set to at least check if BRIDGE is better
  swap_set_t candidate_swaps = {
      swap,
      {swap.first,
       swap.first}};  // second swap here will just compare the base case

  // as with best swap finder, we create a set of candidate swap gates and
  // then find best, except with only 2 swap (best swap and no swap)
  while (candidate_swaps.size() > 1 /*some lookahead parameter*/) {
    this->mapping_frontier_->advance_next_2qb_slice(lookahead);
    // true bool means it only sets interacting uids if both uids are in
    // architecture
    this->set_interacting_uids(true);
    // if 0, just take first swap rather than place
    if (this->interacting_uids_.size() == 0) {
      candidate_swaps = {*candidate_swaps.begin()};
    } else {
      interacting_nodes_t convert_uids;
      for (const auto& p : this->interacting_uids_) {
        convert_uids.insert(
            {Node(this->labelling_[p.first]),
             Node(this->labelling_[p.second])});
      }
      LexicographicalComparison lookahead_lc(this->architecture_, convert_uids);
      lookahead_lc.remove_swaps_lexicographical(candidate_swaps);
    }
  }
  // condition implies bridge is chosen
  // if both remained then lexicographically equivalent under given conditions
  // so either can be added with same consequences (for given hyper
  // parameters)
  if (*candidate_swaps.begin() == swap) {
    output = {0, 0};
  }
  return output;
}

// Returns the distance between n1 and p1 and the distance between n2 and p2,
// distance ordered (greatest first)
const std::pair<size_t, size_t> LexiRoute::pair_distances(
    const Node& p0_first, const Node& p0_second, const Node& p1_first,
    const Node& p1_second) const {
  {
    const bool valid = this->architecture_->node_exists(p0_first) &&
                       this->architecture_->node_exists(p0_second) &&
                       this->architecture_->node_exists(p1_first) &&
                       this->architecture_->node_exists(p1_second);
    TKET_ASSERT(valid);
  }
  size_t curr_dist1 = this->architecture_->get_distance(p0_first, p0_second);
  size_t curr_dist2 = this->architecture_->get_distance(p1_first, p1_second);
  return (curr_dist1 > curr_dist2) ? std::make_pair(curr_dist1, curr_dist2)
                                   : std::make_pair(curr_dist2, curr_dist1);
}

void LexiRoute::remove_swaps_decreasing(swap_set_t& swaps) {
  swap_set_t remaining_swaps;
  Node pair_first, pair_second;
  for (const auto& swap : swaps) {
    auto it = this->interacting_uids_.find(swap.first);
    // => swap.first is in interaction
    if (it != this->interacting_uids_.end()) {
      // find its pair
      pair_first = Node(it->second);
    } else {
      // => not interacting, assign pair to self (will give lexicographic
      // distance 0)
      pair_first = swap.first;
    }
    // => UnitID in SWAP are interacting
    if (pair_first == swap.second) {
      continue;
    }
    auto jt = this->interacting_uids_.find(swap.second);
    // => swap.second is in interaction
    if (jt != this->interacting_uids_.end()) {
      pair_second = Node(jt->second);
    } else {
      pair_second = swap.second;
    }
    // => UnitID in SWAP are interacting
    // Check should alrady be done with earlier continue
    TKET_ASSERT(pair_second != swap.first);

    const std::pair<size_t, size_t>& curr_dists =
        this->pair_distances(swap.first, pair_first, swap.second, pair_second);
    const std::pair<size_t, size_t>& news_dists =
        this->pair_distances(swap.second, pair_first, swap.first, pair_second);
    if (news_dists >= curr_dists) {
      continue;
    }
    remaining_swaps.insert(swap);
  }
}

void LexiRoute::solve_labelling() {
  this->update_labelling();
  this->mapping_frontier_->update_quantum_boundary_uids(this->labelling_);
  return;
}

void LexiRoute::solve(unsigned lookahead) {
  // store a copy of the original this->mapping_frontier_->quantum_boundray
  // this object will be updated and reset throughout the swap picking procedure
  // so need to return it to original setting at end
  unit_vertport_frontier_t copy;
  for (const std::pair<UnitID, VertPort>& pair :
       this->mapping_frontier_->quantum_boundary->get<TagKey>()) {
    copy.insert({pair.first, pair.second});
  }
  swap_set_t candidate_swaps = this->get_candidate_swaps();
  this->remove_swaps_decreasing(candidate_swaps);
  TKET_ASSERT(candidate_swaps.size() != 0);
  // Only want to substitute a single swap
  // check next layer of interacting qubits and remove swaps until only one
  // lexicographically superior swap is left
  unsigned counter = 0;
  while (candidate_swaps.size() > 1 && counter < lookahead) {
    // if 0, just take first swap rather than place
    if (this->interacting_uids_.size() == 0) {
      break;
    } else {
      interacting_nodes_t convert_uids;
      for (const auto& p : this->interacting_uids_) {
        convert_uids.insert(
            {Node(this->labelling_[p.first]),
             Node(this->labelling_[p.second])});
      }
      LexicographicalComparison lookahead_lc(this->architecture_, convert_uids);
      lookahead_lc.remove_swaps_lexicographical(candidate_swaps);
    }
    counter++;
    this->mapping_frontier_->advance_next_2qb_slice(lookahead);
    // true bool means it only sets interacting uids if both uids are in
    // architecture
    this->set_interacting_uids(true);
  }
  // find best swap
  auto it = candidate_swaps.end();
  --it;

  std::pair<Node, Node> chosen_swap = *it;
  this->mapping_frontier_->set_quantum_boundary(copy);

  this->set_interacting_uids();
  std::pair<bool, bool> check = this->check_bridge(chosen_swap, lookahead);
  // set for final time, to allow gates to be correctly inserted, but then leave
  // as is
  // insert gates
  this->mapping_frontier_->set_quantum_boundary(copy);
  if (!check.first && !check.second) {
    // update circuit with new swap
    // final_labelling is initial labelling permuted by single swap
    this->mapping_frontier_->add_swap(chosen_swap.first, chosen_swap.second);
  } else {
    // only need to reset in bridge case
    this->set_interacting_uids();

    auto add_ordered_bridge = [&](const Node& n) {
      auto it0 = this->mapping_frontier_->quantum_boundary->find(n);
      // this should implicitly be the case if this logic is reached
      TKET_ASSERT(it0 != this->mapping_frontier_->quantum_boundary->end());

      Node other_node = Node(this->interacting_uids_[n]);
      auto it1 = this->mapping_frontier_->quantum_boundary->find(other_node);
      // this should implicitly be the case if this logic is reached
      TKET_ASSERT(it1 != this->mapping_frontier_->quantum_boundary->end());

      auto path = this->architecture_->get_path(n, other_node);
      Node central = Node(path[1]);

      Edge n_edge = this->mapping_frontier_->circuit_.get_nth_out_edge(
          it0->second.first, it0->second.second);
      Edge other_edge = this->mapping_frontier_->circuit_.get_nth_out_edge(
          it1->second.first, it1->second.second);

      unsigned port0 =
          this->mapping_frontier_->circuit_.get_target_port(n_edge);
      unsigned port1 =
          this->mapping_frontier_->circuit_.get_target_port(other_edge);
      // compare port ordering to get control vs target
      TKET_ASSERT(port0 != port1);
      if (port0 < port1) {
        this->mapping_frontier_->add_bridge(n, central, other_node);
      } else {
        this->mapping_frontier_->add_bridge(other_node, central, n);
      }
    };

    if (check.first) {
      add_ordered_bridge(chosen_swap.first);
    }
    if (check.second) {
      add_ordered_bridge(chosen_swap.second);
    }
  }
  return;
}

LexiRouteRoutingMethod::LexiRouteRoutingMethod(unsigned _max_depth)
    : max_depth_(_max_depth){};

bool LexiRouteRoutingMethod::check_method(
    const std::shared_ptr<MappingFrontier>& mapping_frontier,
    const ArchitecturePtr& architecture) const {
  std::set<Vertex> unplaced;
  for (const std::pair<UnitID, VertPort>& pair :
       mapping_frontier->quantum_boundary->get<TagKey>()) {
    if ((mapping_frontier->circuit_.n_in_edges_of_type(
             pair.second.first, EdgeType::Quantum) > 2 &&
         mapping_frontier->circuit_.get_OpType_from_Vertex(pair.second.first) !=
             OpType::BRIDGE)) {
      return false;
    } else if (!architecture->node_exists(Node(pair.first))) {
      // if multi-qubit vertex doesn't have all edges in frontier then
      // won't be check in routing_method anyway
      if (unplaced.find(pair.second.first) == unplaced.end()) {
        unplaced.insert(pair.second.first);
      } else {
        return false;
      }
    }
  }
  return true;
}

unit_map_t LexiRouteRoutingMethod::routing_method(
    std::shared_ptr<MappingFrontier>& mapping_frontier,
    const ArchitecturePtr& architecture) const {
  LexiRoute lr(architecture, mapping_frontier);
  lr.solve(this->max_depth_);
  return {};
}

unsigned LexiRouteRoutingMethod::get_max_depth() const {
  return this->max_depth_;
}

nlohmann::json LexiRouteRoutingMethod::serialize() const {
  nlohmann::json j;
  j["depth"] = this->get_max_depth();
  j["name"] = "LexiRouteRoutingMethod";
  return j;
}

LexiRouteRoutingMethod LexiRouteRoutingMethod::deserialize(
    const nlohmann::json& j) {
  return LexiRouteRoutingMethod(j.at("depth").get<unsigned>());
}

}  // namespace tket
