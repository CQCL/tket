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

#include "ContextualReduction.hpp"

#include <algorithm>
#include <optional>
#include <sstream>

#include "Circuit/Circuit.hpp"
#include "Circuit/DAGDefs.hpp"
#include "Eigen/src/Core/Matrix.h"
#include "OpType/OpType.hpp"
#include "Ops/ClassicalOps.hpp"
#include "Ops/OpPtr.hpp"
#include "Transform.hpp"
#include "Utils/Assert.hpp"
#include "Utils/Exceptions.hpp"
#include "Utils/HelperFunctions.hpp"
#include "Utils/UnitID.hpp"

namespace tket {

namespace Transforms {

Transform remove_discarded_ops() {
  return Transform([](Circuit &circ) {
    // We want to keep all vertices that have an Output or ClOutput in their
    // causal future. Start by constructing this set, then remove the
    // remainder.
    VertexSet keep;
    for (auto v_end : circ.all_outputs()) {
      if (circ.get_OpType_from_Vertex(v_end) != OpType::Discard) {
        // Trace back from v, adding to keep-set. We can ignore vertices
        // that are already in the keep-set.
        VertexSet v_frontier;
        keep.insert(v_end);
        v_frontier.insert(v_end);
        while (!v_frontier.empty()) {
          VertexSet new_v_frontier;
          for (auto v : v_frontier) {
            for (const Vertex &v0 : circ.get_predecessors(v)) {
              if (keep.find(v0) == keep.end()) {
                keep.insert(v0);
                new_v_frontier.insert(v0);
              }
            }
          }
          v_frontier = std::move(new_v_frontier);
        }
      }
    }
    // Now remove all operations not in the keep-set.
    VertexList to_remove;
    BGL_FORALL_VERTICES(v, circ.dag, DAG) {
      if (keep.find(v) == keep.end()) {
        OpType optype = circ.get_OpType_from_Vertex(v);
        if (is_gate_type(optype) || is_box_type(optype)) {
          to_remove.push_back(v);
        }
      }
    }
    circ.remove_vertices(
        to_remove, Circuit::GraphRewiring::Yes, Circuit::VertexDeletion::Yes);
    return !to_remove.empty();
  });
}

static std::optional<Eigen::MatrixXcd> op_unitary(Op_ptr op) {
  try {
    return op->get_unitary();
  } catch (NotValid &) {
    return std::nullopt;
  } catch (NotImplemented &) {
    return std::nullopt;
  }
  // These cover the "expected" exceptions; anything else is unexpected.
}

/**
 * Return i s.t. |U[i,j]| = 1, if it exists. U is unitary, j a column index.
 */
static std::optional<unsigned> unique_unit_row(
    const Eigen::MatrixXcd U, unsigned j) {
  unsigned n = U.rows();
  for (unsigned i = 0; i < n; i++) {
    double a = std::abs(U(i, j));
    if (std::abs(a - 1) < EPS) {
      return i;
    } else if (a >= EPS) {
      return std::nullopt;
    }
  }
  // If we get here there is a logical error or numerical instability.
  std::stringstream ss;
  ss << U;
  throw NotUnitary(ss.str());
}

/**
 * Return the set of v in F all of whose quantum in-edges are in qvals.
 */
static VertexSet known_inputs_only(
    const Circuit &circ, const VertexSet &F,
    const std::map<Edge, bool> &qvals) {
  VertexSet v_frontier;
  for (const Vertex &v : F) {
    EdgeVec v_q_inedges = circ.get_in_edges_of_type(v, EdgeType::Quantum);
    if (std::all_of(
            v_q_inedges.begin(), v_q_inedges.end(),
            [&qvals](const Edge &e) { return qvals.find(e) != qvals.end(); })) {
      v_frontier.insert(v);
    }
  }
  return v_frontier;
}

Transform simplify_initial(
    AllowClassical allow_classical, CreateAllQubits create_all_qubits,
    std::shared_ptr<const Circuit> xcirc) {
  return Transform([allow_classical, create_all_qubits, xcirc](Circuit &circ) {
    if (create_all_qubits == CreateAllQubits::Yes) {
      circ.qubit_create_all();
    }

    // Find all Create and Reset vertices.
    VertexSet zeroing_vertices;
    BGL_FORALL_VERTICES(v, circ.dag, DAG) {
      OpType optype = circ.get_OpType_from_Vertex(v);
      if (optype == OpType::Create || optype == OpType::Reset) {
        zeroing_vertices.insert(v);
      }
    }

    // Partial map from quantum edges to values.
    std::map<Edge, bool> qvals;

    // Assign 0 values to all edges coming out of zeroing vertices, and
    // construct the set of their target vertices.
    VertexSet F;
    for (auto z : zeroing_vertices) {
      EdgeVec z_outedges = circ.get_all_out_edges(z);
      TKET_ASSERT(z_outedges.size() == 1);
      Edge z_out = z_outedges[0];
      TKET_ASSERT(circ.get_edgetype(z_out) == EdgeType::Quantum);
      qvals[z_out] = false;
      F.insert(circ.target(z_out));
    }

    // Construct frontier of vertices with all-known input values.
    VertexSet v_frontier = known_inputs_only(circ, F, qvals);

    // Partial map from vertices to sequences of X gates to replace them.
    std::map<Vertex, std::vector<bool>> reductions;

    // Partial map from Measure vertices to bits to set after measure.
    std::map<Vertex, bool> measurebits;

    while (!v_frontier.empty()) {
      // Simplify vertices in the frontier that we can, and move it on.
      F.clear();
      for (const Vertex &v : v_frontier) {
        // If there are any Boolean inputs to v, skip it.
        if (circ.n_in_edges_of_type(v, EdgeType::Boolean) != 0) {
          continue;
        }

        Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
        std::optional<Eigen::MatrixXcd> U = op_unitary(op);
        bool is_measure = (op->get_type() == OpType::Measure);

        if (!U && (!is_measure || allow_classical == AllowClassical::No)) {
          continue;
        }

        // Compute input state.
        EdgeVec v_q_inedges = circ.get_in_edges_of_type(v, EdgeType::Quantum);
        unsigned n_q = v_q_inedges.size();
        std::vector<bool> v_invals(n_q);
        for (unsigned i = 0; i < n_q; i++) {
          v_invals[i] = qvals[v_q_inedges[i]];
        }

        if (U) {
          // Compute relevant column J of U.
          unsigned J = 0, pow2 = 1u << n_q;
          for (unsigned i = 0; i < n_q; i++) {
            pow2 >>= 1;
            if (v_invals[i]) {
              J |= pow2;
            }
          }

          // Check if there is a unique I s.t. |U(I,J)| == 1.
          std::optional<unsigned> I = unique_unit_row(*U, J);
          if (!I) continue;

          // Convert I to a vector of bool
          std::vector<bool> v_outvals(n_q);
          for (unsigned i = 0; i < n_q; i++) {
            v_outvals[i] = (*I >> (n_q - 1 - i)) & 1;
          }

          // Label out-edges of v; construct equivalent X-gate rep.
          EdgeVec v_outedges = circ.get_all_out_edges(v);
          TKET_ASSERT(v_outedges.size() == n_q);
          std::vector<bool> x_gates(n_q);
          for (unsigned i = 0; i < n_q; i++) {
            bool outval = v_outvals[i];
            qvals[v_outedges[i]] = outval;
            x_gates[i] = v_invals[i] ^ outval;
          }

          // Record vertex for later replacement with X-gates.
          reductions[v] = x_gates;

          // Add all successors of v to potential next frontier F.
          for (const Edge &e : v_outedges) {
            F.insert(circ.target(e));
          }
        } else {
          TKET_ASSERT(allow_classical == AllowClassical::Yes);
          TKET_ASSERT(n_q == 1);
          measurebits[v] = v_invals[0];
        }
      }

      // Replace v_frontier with vertices from F having all-known inputs.
      v_frontier = known_inputs_only(circ, F, qvals);
    }

    // Perform substitutions.
    VertexList bin;
    for (const auto &[v, x_gates] : reductions) {
      unsigned n_q = x_gates.size();
      Circuit xs_circ(n_q);
      for (unsigned i = 0; i < n_q; i++) {
        if (x_gates[i]) {
          if (xcirc) {
            unit_map_t qm;
            qm.insert({Qubit(i), Qubit(0)});
            xs_circ.append_with_map(*xcirc, qm);
          } else {
            xs_circ.add_op<unsigned>(OpType::X, {i});
          }
        }
      }
      EdgeVec v_in_edges = circ.get_in_edges(v);        // all Quantum
      EdgeVec v_out_edges = circ.get_all_out_edges(v);  // all Quantum
      TKET_ASSERT(v_in_edges.size() == n_q);
      TKET_ASSERT(v_out_edges.size() == n_q);
      Subcircuit subc{v_in_edges, v_out_edges, {v}};
      circ.substitute(xs_circ, subc, Circuit::VertexDeletion::No);
      bin.push_back(v);
    }
    for (const auto &[v, bitval] : measurebits) {
      Circuit newc(1, 1);
      // Replace the measure with a set-bit.
      std::vector<bool> values = {bitval};
      Op_ptr setbitop = std::make_shared<SetBitsOp>(values);
      newc.add_op<Bit>(setbitop, {Bit(0)});
      EdgeVec q_in_edges = circ.get_in_edges_of_type(v, EdgeType::Quantum);
      TKET_ASSERT(q_in_edges.size() == 1);
      EdgeVec q_out_edges = circ.get_out_edges_of_type(v, EdgeType::Quantum);
      TKET_ASSERT(q_out_edges.size() == 1);
      EdgeVec c_in_edges = circ.get_in_edges_of_type(v, EdgeType::Classical);
      TKET_ASSERT(c_in_edges.size() == 1);
      EdgeVec c_out_edges = circ.get_out_edges_of_type(v, EdgeType::Classical);
      TKET_ASSERT(c_out_edges.size() == 1);
      EdgeVec b_out_edges = circ.get_out_edges_of_type(v, EdgeType::Boolean);
      Subcircuit subc{q_in_edges,  q_out_edges, c_in_edges,
                      c_out_edges, b_out_edges, {v}};
      circ.substitute(newc, subc, Circuit::VertexDeletion::No);
      bin.push_back(v);
    }
    bool changed = !bin.empty();
    circ.remove_vertices(
        bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);

    return changed;
  });
}

std::optional<std::shared_ptr<const ClassicalTransformOp>> classical_transform(
    Op_ptr op) {
  std::optional<Eigen::MatrixXcd> U = op_unitary(op);
  if (!U) return std::nullopt;
  unsigned n = op->get_desc().n_qubits().value();
  unsigned pow2n = 1u << n;
  TKET_ASSERT(U->cols() == pow2n);
  std::vector<uint32_t> values(pow2n);
  for (uint32_t J = 0; J < pow2n; J++) {
    // Look at the column U[*,J]. Is there a unique nonzero element?
    std::optional<unsigned> I = unique_unit_row(*U, J);
    if (!I) return std::nullopt;
    values[reverse_bits(J, n)] = reverse_bits(*I, n);
  }
  return std::make_shared<ClassicalTransformOp>(n, values);
}

Transform simplify_measured() {
  return Transform([](Circuit &circ) {
    // First construct the set of all Measure vertices that are followed by
    // Discard vertices (and no Boolean out-edges).
    VertexSet M;
    BGL_FORALL_VERTICES(v, circ.dag, DAG) {
      if (circ.get_OpType_from_Vertex(v) == OpType::Measure) {
        if (circ.n_out_edges_of_type(v, EdgeType::Boolean) == 0) {
          EdgeVec m_q_outs = circ.get_out_edges_of_type(v, EdgeType::Quantum);
          TKET_ASSERT(m_q_outs.size() == 1);
          Edge m_q_out = m_q_outs[0];
          Vertex v1 = circ.target(m_q_out);
          if (circ.get_OpType_from_Vertex(v1) == OpType::Discard) {
            M.insert(v);
          }
        }
      }
    }
    bool changed = false;
    bool carry_on;
    do {
      carry_on = false;
      VertexList bin;
      // Find all classical maps all of whose successors are in M
      for (const Vertex &v : M) {
        VertexVec preds = circ.get_predecessors(v);
        for (const Vertex &v0 : preds) {
          // Any Boolean inputs?
          if (circ.n_in_edges_of_type(v0, EdgeType::Boolean) != 0) {
            continue;
          }
          // No. Are all successors in M?
          VertexVec succs = circ.get_successors(v0);
          if (std::any_of(succs.begin(), succs.end(), [&](const Vertex &u) {
                return M.find(u) == M.end();
              })) {
            continue;
          }
          // Yes. Is it a classical map?
          Op_ptr op = circ.get_Op_ptr_from_Vertex(v0);
          std::optional<std::shared_ptr<const ClassicalTransformOp>> cm =
              classical_transform(op);
          if (!cm) continue;
          // Yes. Remove v0.
          unsigned n_qb = succs.size();
          circ.remove_vertex(
              v0, Circuit::GraphRewiring::Yes, Circuit::VertexDeletion::No);
          bin.push_back(v0);
          // Insert *cm on the classical target wires.
          EdgeVec cl_edges(n_qb);
          for (unsigned i = 0; i < n_qb; i++) {
            EdgeVec m_c_outs =
                circ.get_out_edges_of_type(succs[i], EdgeType::Classical);
            TKET_ASSERT(m_c_outs.size() == 1);
            cl_edges[i] = m_c_outs[0];
          }
          Subcircuit cl_subc{{}, {}, cl_edges, cl_edges, {}, {}};
          Circuit cl_circ(0, n_qb);
          std::vector<unsigned> args(n_qb);
          std::iota(args.begin(), args.end(), 0);
          cl_circ.add_op<unsigned>(*cm, args);
          circ.substitute(cl_circ, cl_subc, Circuit::VertexDeletion::No);
          changed = true;
          carry_on = true;
        }
      }
      circ.remove_vertices(
          bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
    } while (carry_on);
    return changed;
  });
}

std::pair<Circuit, Circuit> separate_classical(const Circuit &circ) {
  // Initialize the two circuits to return.
  qubit_vector_t qubits = circ.all_qubits();
  bit_vector_t bits = circ.all_bits();
  Circuit c0, c1;
  for (const Qubit &qb : qubits) {
    c0.add_qubit(qb);
  }
  for (const Bit &b : bits) {
    c0.add_bit(b);
    c1.add_bit(b);
  }

  // Get the command list.
  std::vector<Command> cmds = circ.get_commands();

  // Find the final Measure (or ClInput if no Measure) on each bit.
  VertexVec c_in = circ.c_inputs();
  unsigned n_bits = bits.size();
  TKET_ASSERT(n_bits == c_in.size());
  std::map<Bit, Vertex> final_vert_on_bit;
  for (unsigned i = 0; i < n_bits; i++) {
    final_vert_on_bit[bits[i]] = c_in[i];
  }
  for (const Command &cmd : cmds) {
    OpType optype = cmd.get_op_ptr()->get_type();
    if (optype == OpType::Measure) {
      bit_vector_t cmd_bits = cmd.get_bits();
      TKET_ASSERT(cmd_bits.size() == 1);
      final_vert_on_bit[cmd_bits[0]] = cmd.get_vertex();
    }
  }

  // Construct the set of final vertices.
  VertexSet final_verts;
  for (const Bit &b : bits) {
    final_verts.insert(final_vert_on_bit[b]);
  }

  // Construct the set of vertices to go in c1.
  // Initially include final_verts; later exclude them.
  VertexSet c1_verts(final_verts);
  for (const Command &cmd : cmds) {
    Vertex v = cmd.get_vertex();
    VertexVec preds = circ.get_predecessors(v);
    if (std::all_of(preds.begin(), preds.end(), [&c1_verts](const Vertex &v) {
          return c1_verts.find(v) != c1_verts.end();
        })) {
      c1_verts.insert(v);
    }
  }
  for (const Vertex &v : final_verts) {
    c1_verts.erase(v);
  }

  // Step through the circuit, filling in c0 and c1.
  for (const Command &cmd : cmds) {
    Op_ptr op = cmd.get_op_ptr();
    unit_vector_t args = cmd.get_args();
    if (c1_verts.find(cmd.get_vertex()) == c1_verts.end()) {
      c0.add_op(op, args);
    } else {
      c1.add_op(op, args);
    }
  }

  return {c0, c1};
}

}  // namespace Transforms

}  // namespace tket
