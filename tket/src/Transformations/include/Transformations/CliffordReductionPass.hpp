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

#include "Transform.hpp"

namespace tket {

/**
 * Represents a position that one port of a 2qb Clifford gadget could be
 * commuted to (towards the end of the circuit).
 */
struct InteractionPoint {
  /**
   * Edge to which one point of the interaction could potentially be moved.
   * This is always in the causal future of the source vertex.
   */
  Edge e;

  /** Vertex representing the current interaction */
  Vertex source;

  /** What the Pauli term would be transformed to if moved to new location */
  Pauli p;

  /** Phase of source (true for -pi/2, false for +pi/2) */
  bool phase;  // For each source of an interaction, True for -pi/2, false for
               // +pi/2

  std::pair<Edge, Vertex> key() const { return {e, source}; }
};

/**
 * Represents a position that one port of a 2qb Clifford gadget could be
 * commuted to (towards the start of the circuit). We don't need to store the
 * source since we only consider commuting one gadget backwards at a time to
 * search for matches.
 */
struct RevInteractionPoint {
  /**
   * Edge to which one point of the interaction could potentially be moved.
   * This is always in the causal past of the source vertex.
   */
  Edge e;

  /** What the Pauli term would be transformed to if moved to new location */
  Pauli p;

  /** Phase of source (true for -pi/2, false for +pi/2) */
  bool phase;
};

/**
 * Represents a position to which two 2q Clifford gadgets could be commuted to
 * Assumes:
 * - point0.source == point1.source
 * - point0.e == rev0.e
 * - point1.e == rev1.e
 */
struct InteractionMatch {
  InteractionPoint point0;
  InteractionPoint point1;
  RevInteractionPoint rev0;
  RevInteractionPoint rev1;
};

struct TagEdge {};
struct TagSource {};

/**
 * Stores all current InteractionPoints encountered.
 * Each gadget can only commute to a given edge in a unique way.
 * Searching by edge is needed to search for possible matches.
 * Searching by vertex is needed to remove InteractionPoints when a pair merge.
 */
typedef boost::multi_index::multi_index_container<
    InteractionPoint,
    boost::multi_index::indexed_by<
        boost::multi_index::hashed_unique<
            boost::multi_index::tag<TagKey>,
            boost::multi_index::const_mem_fun<
                InteractionPoint, std::pair<Edge, Vertex>,
                &InteractionPoint::key>>,
        boost::multi_index::hashed_non_unique<
            boost::multi_index::tag<TagEdge>,
            boost::multi_index::member<
                InteractionPoint, Edge, &InteractionPoint::e>>,
        boost::multi_index::hashed_non_unique<
            boost::multi_index::tag<TagSource>,
            boost::multi_index::member<
                InteractionPoint, Vertex, &InteractionPoint::source>>>>
    interaction_table_t;

class CliffordReductionPass {
 public:
  /** Apply the transformation to a circuit */
  static bool reduce_circuit(Circuit &circ, bool allow_swaps = false);

 private:
  /** The circuit under transformation */
  Circuit &circ;

  /** Table of potential transports of 2qb interactions through the circuit */
  interaction_table_t itable;

  /** Map from vertex to depth in circuit */
  std::map<Vertex, unsigned> v_to_depth;

  /** Map from vertex to sets of units that it acts on */
  std::map<Vertex, unit_set_t> v_to_units;

  /** Map from edge to corresponding unit */
  std::map<Edge, UnitID> e_to_unit;

  /** Whether any changes have been made to the circuit */
  bool success;

  /** Current depth in forward pass through circuit. */
  unsigned current_depth;

  /** Whether to allow SWAP gates in the final circuit */
  bool allow_swaps;

  /**
   * Inserts a new interaction point into the table and tries to commute it as
   * far forwards as possible, taking into account commutations and
   * change-of-basis through 1q Cliffords. Does not change the circuit, but
   * extends the table of interaction points as it goes. Stops propagating
   * when it reaches an interaction point that is already in the table.
   */
  void insert_interaction_point(InteractionPoint ip);

  /**
   * Tries to commute an interaction backwards to search for a match,
   * consulting the current table of interaction points to find canditates.
   */
  std::optional<InteractionMatch> search_back_for_match(
      const RevInteractionPoint &rip0, const RevInteractionPoint &rip1) const;

  /**
   * Called when a new 2q Clifford gate is observed (either encountered by
   * advancing through the circuit or by substitution mid-circuit). Handles
   * scanning for pattern matches and substituting, and adding new
   * InteractionPoints to look for matches further into the circuit.
   */
  void process_new_interaction(const Vertex &v);

  /**
   * Inserts a circuit at the given point, updating itable, v_to_depth,
   * v_to_units, and e_to_unit appropriately.
   *
   * @pre Depth of `to_replace` is at most 1.
   */
  Subcircuit substitute(const Circuit &to_insert, const Subcircuit &to_replace);

  /**
   * Given a set of candidate edges, finds the earliest one that is in the
   * strict causal future of source. Returns nullopt if none are in the future
   * of source. Search threshold is set by only having v_to_depth defined for
   * the vertices to search.
   */
  std::optional<Edge> find_earliest_successor(
      const Edge &source, const EdgeSet &candidates) const;

  /**
   * Given two lists of edges, identifies a cut where a two-qubit operation
   * could feasibly be introduced without breaking the causal structure (i.e.
   * introducing loops in the DAG). Returns nullopt if there is no such cut,
   * otherwise will yield the future-most valid point with respect to each of
   * the sequences.
   */
  std::optional<std::pair<InteractionPoint, InteractionPoint>>
  valid_insertion_point(
      const std::list<InteractionPoint> &seq0,
      const std::list<InteractionPoint> &seq1) const;

  CliffordReductionPass(Circuit &c, bool swaps);

  friend class CliffordReductionPassTester;
};

/** a simple friend class of `CliffordReductionPass` used for testing
 * private member functions */
class CliffordReductionPassTester {
 public:
  explicit CliffordReductionPassTester(Circuit &circ);

  std::optional<std::pair<InteractionPoint, InteractionPoint>>
  valid_insertion_point(
      const std::list<InteractionPoint> &seq0,
      const std::list<InteractionPoint> &seq1) const;

 private:
  CliffordReductionPass context;
};

namespace Transforms {

Transform clifford_reduction(bool allow_swaps = false);

}

}  // namespace tket
