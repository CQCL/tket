#ifndef _TKET_TokenSwapping_main_entry_functions_H_
#define _TKET_TokenSwapping_main_entry_functions_H_

#include <map>
#include <tuple>
#include <utility>
#include <vector>

#include "Architecture/Architectures.hpp"
#include "Circuit/Circuit.hpp"

namespace tket {

/** This specifies desired source->target vertex mappings.
 *  Any nodes not occurring as a key might be moved by the algorithm.
 */
typedef std::map<Node, Node> NodeMapping;

/** Version 1.1, not too bad.
 *  @param architecture The raw object containing the graph.
 *  @param node_mapping The desired source->target node mapping.
 *  @return The required list of node pairs to swap.
 */
std::vector<std::pair<Node, Node>> get_swaps(
    const Architecture& architecture, const NodeMapping& node_mapping);

/** An alternative interface, which just wraps the other "get_swaps" function.
 *  In the returned tuple, the Circuit implements using SWAP gates,
 *  and the unit_map_t objects are the initial and final mappings of
 *  logical qubits to architecture nodes.
 *  NOTE: the architecture may contain other nodes not mentioned in the
 *  input logical->physical maps, which may get moved.
 *  If you don't want this, you must include these nodes in the maps.
 *  @param architecture The architecture (containing nodes, and edges)
 *  @param initial_logical_to_physical_map The key is the initial logical qubit,
 *              the value is the existing physical node in the architecture
 *              which it currently maps to.
 *  @param desired_logical_to_physical_map The keys are the same logical qubits
 *              as in "initial_logical_to_physical_map", but the values are now
 *              the nodes where we want them to map AFTER the swaps.
 *  @return A circuit containing the swaps (SWAP gates only), plus the resultant
 *      logical to physical mappings before and after (necessarily the same as
 *      the input mappings, because the returned swaps should always result
 *      in the desired end-to-end mapping exactly).
 */
std::tuple<Circuit, unit_map_t, unit_map_t> get_swaps(
    const Architecture& architecture,
    const unit_map_t& initial_logical_to_physical_map,
    const unit_map_t& desired_logical_to_physical_map);

}  // namespace tket
#endif
