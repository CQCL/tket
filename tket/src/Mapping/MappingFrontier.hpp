#ifndef _TKET_MappingFrontier_H_
#define _TKET_MappingFrontier_H_

#include "Architecture/Architecture.hpp"
#include "Circuit/Circuit.hpp"
#include "Utils/BiMapHeaders.hpp"
#include "Utils/UnitID.hpp"

namespace tket {

typedef sequenced_map_t<UnitID, VertPort> unit_vertport_frontier_t;

// list of error types to throw out
class MappingFrontierError : public std::logic_error {
 public:
  explicit MappingFrontierError(const std::string& message)
      : std::logic_error(message) {}
};

struct MappingFrontier {
  /**
   * VertPort instead of Edge as Edge changes in substitution, but Vertex and
   * Port key information
   */
  std::shared_ptr<unit_vertport_frontier_t> quantum_boundary;

  std::shared_ptr<b_frontier_t> classical_boundary;

  /**
   * Circuit held by reference and directly modified with SWAP (or other
   * relevant) gates.
   */
  Circuit& circuit_;

  std::set<Node> ancilla_nodes_;

  MappingFrontier(Circuit& _circuit);

  /**
   * Given some Circuit Cut (or routed/unrouted boundary), advances the cut to
   * the next cut of just two-qubit vertices, not including the current
   * boundary.
   * @param max_advance maximum number of cuts checked before terminating
   */
  void advance_next_2qb_slice(unsigned max_advance);

  /**
   * mapping_frontier data members updated to reflect
   * the routed/non-routed boundary of mapping_frontier->circ
   * architecture.valid_gate confirms whether circuit vertices are physically
   * valid
   *
   * @param architecture Architecture governing physically allowed operations
   */
  void advance_frontier_boundary(const ArchitecturePtr& architecture);

  /**
   * Subcircuit produced from gates after held boundary.
   * @param _max_subcircuit_depth
   * @param _max_subcircuit_size
   *
   */
  Subcircuit get_frontier_subcircuit(
      unsigned _max_subcircuit_depth, unsigned _max_subcircuit_size) const;

  /**
   * update_quantum_boundary_uids
   * route_circuit has no constraint that passed circuits must have qubits
   * relabelled to architecture nodes route_subcircuit is allowed to either
   * permute labelled physical qubits, or label logical qubits if logical qubits
   * are labelled physical, update_quantum_boundary updates UnitID in
   * this->quantum_boundary to reflect this change Also updates this->circuit_
   * to reflect this relabelling
   *
   * @param relabel_map map between current UnitID's in quantum_boundary and new
   * UnitID's.
   */
  void update_quantum_boundary_uids(const unit_map_t& relabel_map);

  /**
   * permute_subcircuit_q_out_hole
   *
   * Given initial permutation of UnitIDs, finds final permutation via SWAPs in
   * circuit and updates mapping_frontier subcircuit q_out_hole to reflect this
   *
   * @param final_permutation map between initial and final physical qubits for
   * each logical qubit, used to permute subcircuit.q_out_hole
   * @param subcircuit Subcircuit for rearranging boundary
   */
  void permute_subcircuit_q_out_hole(
      const unit_map_t& final_permutation, Subcircuit& subcircuit);

  /**
   * get_default_to_quantum_boundary_unit_map
   * subcircuit circuits created with default q register
   * method returns map between default q register and physical qubit
   * permutation at frontier used for circuit.rename_units
   */
  unit_map_t get_default_to_quantum_boundary_unit_map() const;

  /**
   * add_swap
   * Inserts an OpType::SWAP gate into the uid_0 and uid_1 edges held in
   * quantum_boundary. This directly modifies circuit_.
   * Updates quantum_boundary to reflect new edges.
   *
   * @param uid_0 First Node in SWAP
   * @param uid_1 Second Node in SWAP
   */
  void add_swap(const UnitID& uid_0, const UnitID& uid_1);

  /**
   * add_bridge
   * Inserts an OpType::BRIDGE gate into edges relevant to passed UnitID.
   *
   * @param control First Node in BRIDGE
   * @param central Second Node in BRIDGE
   * @param target Third Node in BRIDGE
   */
  void add_bridge(
      const UnitID& control, const UnitID& central, const UnitID& target);

  /**
   * add_ancilla
   * Adds an Ancillary UnitID to Circuit and tracked information
   *
   * @param ancilla UnitID of added ancilla
   */
  void add_ancilla(const UnitID& ancilla);

  /**
   * merge_ancilla
   * Rewires this->circuit_.dag such that in wire to ancilla Output vertex
   * is now mapped to out wire of merge Input vertex.
   *
   * @param merge UnitID to which ancilla path is prepended
   * @param ancilla UnitID of ancilla opeartions
   */
  void merge_ancilla(const UnitID& merge, const UnitID& ancilla);

  /**
   * Assigns the quantum_boundary_ attribute to that passed to method.
   *
   * @param new_boundary Object to reassign with.
   */
  void set_quantum_boundary(const unit_vertport_frontier_t& new_boundary);
};

}  // namespace tket

#endif