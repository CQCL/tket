#include <algorithm>
#include <array>
#include <iostream>
#include <iterator>
#include <memory>
#include <numeric>
#include <set>
#include <string>

#include "Circuit.hpp"
#include "Circuit/QInteraction.hpp"
#include "Circuit/CircUtils.hpp"
#include "Circuit/Circuit.hpp"
#include "Circuit/DAGDefs.hpp"
#include "Circuit/ThreeQubitConversion.hpp"
#include "OpType/EdgeType.hpp"
#include "OpType/OpType.hpp"
#include "Utils/Assert.hpp"
#include "Utils/GraphHeaders.hpp"

namespace tket {

std::vector<Subcircuit> Circuit::get_subcircuits(const unsigned& n_qubit, const unsigned& min_gate_count) {
	
	std::vector<Subcircuit>  subcircuits;

	// Step through the vertices in topological order.
	QISystem Is(*this, [](auto c) { return c;});  // set of "live" interactions
	for (const Vertex &v : vertices_in_order()) {
		const EdgeVec v_q_ins = get_in_edges_of_type(v, EdgeType::Quantum);
		const EdgeVec v_q_outs = get_out_edges_of_type(v, EdgeType::Quantum);
		unsigned n_q_ins = v_q_ins.size();
		unsigned n_q_outs = v_q_outs.size();

		// If v has no quantum wires, ignore it and move on.
		if (n_q_ins == 0 && n_q_outs == 0) continue;

		// If v is initial, create an interaction from its out-edge, and move on.
		if (n_q_ins == 0) {
			TKET_ASSERT(n_q_outs == 1);
			Is.create_new_interaction_from_edge(v_q_outs[0]);
			continue;
		}

		// If v is final, ignore it and move on.
		if (n_q_outs == 0) continue;

		// It's an internal operation with >0 quantum wires.
		TKET_ASSERT(n_q_ins == n_q_outs);

		Op_ptr op = get_Op_ptr_from_Vertex(v);
		OpType optype = op->get_type();

		// If there are any incoming classical wires, or if this is a Barrier or
		// Reset or Collapse operation, or if the operation contains symbols,
		// close all existing interactions meeting v and create new ones, then
		// move on.
		if (!get_in_edges_of_type(v, EdgeType::Classical).empty() ||
				!get_in_edges_of_type(v, EdgeType::Boolean).empty() ||
				optype == OpType::Barrier || optype == OpType::Reset ||
				optype == OpType::Collapse || !op->free_symbols().empty()) {
			std::vector<int> v_interactions = Is.interactions_feeding_vertex(v);
			std::map<int, iptr> * all_interactions = Is.get_interactions();
			for (const int& i : v_interactions) {
				if (all_interactions->at(i)->n_vertices() >= min_gate_count) {
					Subcircuit sub = all_interactions->at(i)->subcircuit();
					subcircuits.push_back(sub);
				}
			}
			Is.close_interactions_feeding_vertex(v, false);
			continue;
		}

		// Absorb v into existing interactions, closing or merging as necessary.
		bool done_with_v = false;
		while (!done_with_v) {
			std::vector<int> v_Is = Is.interactions_feeding_vertex(v);
			unsigned total_n_wires = Is.total_n_wires(v_Is);
			if (total_n_wires <= n_qubit) {
				Is.combine_and_append(v_Is, v);
				done_with_v = true;
			} else {
				// Close one of the interactions meeting v.
				int i = Is.largest_interaction(v_Is);
				// Append if it exceeds the size requirement
				std::map<int, iptr> * all_interactions = Is.get_interactions();
				if (all_interactions->at(i)->n_vertices() >= min_gate_count) {
					Subcircuit sub = all_interactions->at(i)->subcircuit();
					subcircuits.push_back(sub);
				}
				Is.close_interaction_and_spawn(i, false);
			}
		}
	}

	// Close all remaining interactions.
	std::map<int, iptr> * all_interactions = Is.get_interactions();
	for (auto it = all_interactions->begin(); it != all_interactions->end(); it++) {
		if (it->second->n_vertices() >= min_gate_count) {
			Subcircuit sub = it->second->subcircuit();
			subcircuits.push_back(sub);
		}
	}

	Is.close_all_interactions(false);

	// Delete removed vertices.
	Is.destroy_bin();

	return subcircuits;
}

}