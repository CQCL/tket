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

#include <cmath>

#include "Architecture/Architecture.hpp"
#include "Circuit/Circuit.hpp"
#include "ComparisonFunctions.hpp"
#include "Predicates/PassGenerators.hpp"
#include "Predicates/PassLibrary.hpp"
#include "Predicates/Predicates.hpp"
#include "Simulation/CircuitSimulator.hpp"
#include "Transformations/ContextualReduction.hpp"
#include "Transformations/Transform.hpp"
#include "Utils/Exceptions.hpp"
#include "rapidcheck.h"

using namespace tket;

#define ALL_PASSES(DO)                    \
  DO(SynthesiseTK)                        \
  DO(SynthesiseTket)                      \
  DO(SynthesiseHQS)                       \
  DO(SynthesiseUMD)                       \
  DO(PeepholeOptimise2Q)                  \
  DO(FullPeepholeOptimise)                \
  DO(RemoveRedundancies)                  \
  DO(CommuteThroughMultis)                \
  DO(DecomposeArbitrarilyControlledGates) \
  DO(KAKDecomposition)                    \
  DO(ThreeQubitSquash)                    \
  DO(DecomposeMultiQubitsCX)              \
  DO(DecomposeSingleQubitsTK1)            \
  DO(DecomposeBoxes)                      \
  DO(DecomposeTK2)                        \
  DO(ComposePhasePolyBoxes)               \
  DO(SquashTK1)                           \
  DO(RebaseTket)                          \
  DO(DecomposeBridges)                    \
  DO(FlattenRegisters)                    \
  DO(RemoveBarriers)                      \
  DO(DelayMeasures)                       \
  DO(GlobalisePhasedX)

// Map from PassPtr to readable name
static const std::map<PassPtr, std::string> passes = {
#define NAMEPASS(p) {p(), #p},
    ALL_PASSES(NAMEPASS)
#undef NAMEPASS
};

static std::vector<unsigned> random_subset(
    const std::vector<unsigned> &v, unsigned k) {
  if (k > v.size()) throw std::logic_error("invalid subset size");
  std::vector<unsigned> rvec(k);
  std::set<unsigned> rset;
  unsigned i = 0;
  while (i < k) {
    unsigned x = *rc::gen::elementOf(v);
    bool added = rset.insert(x).second;
    if (added) {
      rvec[i] = x;
      i++;
    }
  }
  return rvec;
}

static std::vector<Expr> random_params(unsigned k) {
  std::vector<Expr> rvec(k);
  for (unsigned i = 0; i < k; i++) {
    double x = *rc::gen::arbitrary<double>();
    // Constrain to [0,2} to avoid rounding errors arising from enormous
    // values.
    x -= floor(x);
    rvec[i] = 2 * x;
  }
  return rvec;
}

// Generate a random circuit with no classical wires.
static Circuit random_circuit() {
  int n_qb = *rc::gen::inRange(1, 5);
  int n_g = *rc::gen::inRange(0, 16);
  Circuit c(n_qb);
  std::vector<unsigned> qbs(n_qb);
  std::iota(qbs.begin(), qbs.end(), 0);
  int i = 0;
  while (i < n_g) {
    OpType g = *rc::gen::elementOf(all_gate_types());
    if (g == OpType::Measure)  // invalid without classical output
    {
      continue;
    }
    const OpTypeInfo &opinfo = optypeinfo().at(g);
    const std::optional<op_signature_t> sig = opinfo.signature;
    int g_nq;
    if (sig) {
      g_nq =
          std::count(sig.value().begin(), sig.value().end(), EdgeType::Quantum);
    } else {
      g_nq = *rc::gen::inRange(1, n_qb + 1);
    }
    if (g_nq > n_qb) continue;
    unsigned g_np = opinfo.n_params();
    std::vector<unsigned> qb = random_subset(qbs, g_nq);
    std::vector<Expr> params = random_params(g_np);
    // For OpType::TK2, angles must currently be normalised.
    // TODO: Remove this when NormalisedTK2Predicate is implemented.
    if (g == OpType::TK2) {
      TKET_ASSERT(params.size() == 3);
      double p0 = *eval_expr(params[0]);
      double p1 = *eval_expr(params[1]);
      double p2 = *eval_expr(params[2]);
      p0 = std::fmod(p0, 0.5);
      p1 = p0 < EPS ? 0. : std::fmod(p1, p0);
      p2 = p1 < EPS ? 0. : std::fmod(p2, p1);
      p2 *= (*rc::gen::arbitrary<bool>()) ? 1 : -1;
      params[0] = p0;
      params[1] = p1;
      params[2] = p2;
    }
    c.add_op<unsigned>(g, params, qb);
    i++;
  }
  return c;
}

static bool sanity_check(const Circuit &c) {
  // Check that n_gates() gives the right number according to the command list.
  unsigned n_gates = c.n_gates();
  unsigned i = 0;
  for (Circuit::CommandIterator it = c.begin(); it != c.end(); ++it) i++;
  return i == n_gates;
}

static bool verify_n_qubits_for_ops(const Circuit &circ) {
  // Check that n_qubits() gives the right answer for all operations.
  for (const Command &com : circ) {
    Op_ptr op = com.get_op_ptr();
    if (op->n_qubits() != com.get_args().size()) {
      RC_LOG() << "Failure at command " << com << std::endl;
      RC_LOG() << "Op::n_qubits() = " << op->n_qubits() << std::endl;
      return false;
    }
  }
  return true;
}

static std::set<Node> random_node_set() {
  int n = *rc::gen::inRange(1, 10);
  std::set<Node> nodes;
  for (int i = 0; i < n; i++) {
    int idx = *rc::gen::inRange(0, 20);
    nodes.insert(Node("x", idx));
  }
  return nodes;
}

static std::set<std::pair<Node, Node>> random_connected_graph(
    std::set<Node> nodes) {
  int n_nodes = nodes.size();
  std::set<std::pair<Node, Node>> links;
  // First connect all the nodes, one by one. Connect each new nodes to a
  // random already-connected node. This creates a spanning tree.
  std::set<Node> connected;
  for (const auto &node : nodes) {
    if (connected.empty()) {
      connected.insert(node);
    } else {
      Node node0 = *rc::gen::elementOf(connected);
      links.insert({node, node0});
      connected.insert(node);
    }
  }
  // Now add a random selection of new links.
  int n_new = *rc::gen::inRange(0, 2 * n_nodes);
  for (int i = 0; i < n_new; i++) {
    Node node0 = *rc::gen::elementOf(nodes);
    Node node1 = *rc::gen::elementOf(nodes);
    if (node0 != node1) {
      links.insert({node0, node1});
    }
  }
  return links;
}

static Architecture random_architecture() {
  std::set<Node> nodes = random_node_set();
  std::set<std::pair<Node, Node>> links = random_connected_graph(nodes);
  Architecture arc(
      std::vector<std::pair<Node, Node>>(links.begin(), links.end()));
  return arc;
}

/**
 * Check correctness of a completed compilation pass.
 *
 * @param[in] c0 original circuit
 * @param[in] cu compliation pass having been applied to \p c0
 */
static void check_correctness(const Circuit &c0, const CompilationUnit &cu) {
  RC_LOG() << "In Check Correctness" << std::endl;

  const Circuit &c1 = cu.get_circ_ref();
  const unit_bimap_t &initial_map = cu.get_initial_map_ref();
  const unit_bimap_t &final_map = cu.get_final_map_ref();
  // Account for reordering in the initial and final maps.

  std::map<unsigned, unsigned> ini, inv_fin;
  std::map<UnitID, unsigned> c0_idx, c1_idx;
  unsigned i;

  i = 0;
  for (const auto &id : c0.all_units()) {
    c0_idx.insert({id, i});
    i++;
  }
  i = 0;
  for (const auto &id : c1.all_units()) {
    c1_idx.insert({id, i});
    i++;
  }
  TKET_ASSERT(c0_idx.size() <= c1_idx.size());
  Circuit c0_copy(c0);
  for (const auto &pair : initial_map.left) {
    auto it = c0_idx.find(pair.first);
    // qubit not in original circuit => ancilla added
    if (it == c0_idx.end()) {
      TKET_ASSERT(c1_idx.find(pair.first) != c1_idx.end());
      c0_idx.insert({pair.first, c0_idx.size()});
      c0_copy.add_qubit(Qubit(pair.first));
    }
    TKET_ASSERT(c1_idx.find(pair.second) != c1_idx.end());
    ini.insert({c0_idx[pair.first], c1_idx[pair.second]});
  }
  // all qubits should have been tracked from initial map
  for (const auto &pair : final_map.left) {
    inv_fin.insert({c1_idx[pair.second], c0_idx[pair.first]});
  }

  Eigen::PermutationMatrix<Eigen::Dynamic> m_ini = lift_perm(ini),
                                           m_inv_fin = lift_perm(inv_fin);

  try {
    const auto u1 = tket_sim::get_unitary(c1);
    const auto u0 = tket_sim::get_unitary(c0_copy);
    RC_ASSERT(tket_sim::compare_statevectors_or_unitaries(
        u0, m_inv_fin * u1 * m_ini));
  } catch (const Unsupported &) {
  } catch (const NotImplemented &) {
  }
  // If tket-sim doesn't recognise a gate, just ignore it.
  // If a gate is unknown, which exception SHOULD it be:
  // "Unsupported" or "NotImplemented" ?
  RC_ASSERT(sanity_check(c1));
}

namespace rc {

template <>
struct Arbitrary<Circuit> {
  static Gen<Circuit> arbitrary() {
    return gen::exec([] { return random_circuit(); });
  }
};

template <>
struct Arbitrary<PassPtr> {
  static Gen<PassPtr> arbitrary() {
    return gen::exec([] {
      auto p = *rc::gen::elementOf(passes);
      return p.first;
    });
  }
};

template <>
struct Arbitrary<Architecture> {
  static Gen<Architecture> arbitrary() {
    return gen::exec([] { return random_architecture(); });
  }
};

}  // namespace rc

bool check_n_qubits() {
  return rc::check("n_qubits is correct", [](const unsigned &n) {
    unsigned m = n % 20;
    Circuit c(m);
    RC_ASSERT(c.n_qubits() == m);
  });
}

bool check_passes() {
  return rc::check(
      "preconditions and postconditions of passes are correct",
      // Also perform some sanity checks on the circuits before and after
      // the transforms.
      [](const Circuit &c, const PassPtr &p) {
        verify_n_qubits_for_ops(c);
        PassConditions pcons = p->get_conditions();
        PredicatePtrMap precons = pcons.first;
        PredicatePtrMap postcons = pcons.second.specific_postcons_;
        if (std::all_of(precons.begin(), precons.end(), [&c](auto precon) {
              return precon.second->verify(c);
            })) {
          RC_LOG() << "\nCircuit (" << c.n_qubits() << " qubits, "
                   << c.n_gates() << " gates): " << c << std::endl;
          RC_LOG() << "Pass: " << passes.at(p) << std::endl;
          CompilationUnit cu(c);
          bool applied = (p->apply(cu));
          const Circuit &c1 = cu.get_circ_ref();
          verify_n_qubits_for_ops(c1);
          RC_LOG() << "\nNew circuit(" << c1.n_qubits() << " qubits, "
                   << c1.n_gates() << " gates): " << c1 << std::endl;
          if (applied) {
            for (auto postcon : postcons) {
              RC_ASSERT(postcon.second->verify(c1));
            }
            check_correctness(c, cu);
          } else {
            RC_ASSERT(c == c1);
          }
        }
      });
}

bool check_mapping() {
  return rc::check(
      "routing to different architectures",
      [](const Circuit &c, const Architecture &a) {
        // Exclude circuits with classical controls. TKET-235
        PredicatePtr pp1 = std::make_shared<NoClassicalControlPredicate>();
        if (!pp1->verify(c)) return;
        // Architecture must be big enough.
        PredicatePtr pp2 = std::make_shared<MaxNQubitsPredicate>(a.n_nodes());
        if (!pp2->verify(c)) return;
        // All gates must act on 1 or 2 qubits.
        PredicatePtr pp3 = std::make_shared<MaxTwoQubitGatesPredicate>();
        if (!pp3->verify(c)) return;
        PassPtr pass = gen_default_mapping_pass(a, true);
        CompilationUnit cu(c);
        bool applied = pass->apply(cu);
        const Circuit &c1 = cu.get_circ_ref();
        RC_LOG() << "Circuit (" << c.n_qubits() << " qubits, " << c.n_gates()
                 << " gates): " << c;
        RC_LOG() << "Architecture (" << a.n_nodes() << " nodes): ";
        const node_vector_t nodes = a.get_all_nodes_vec();
        for (Node node0 : nodes) {
          for (Node node1 : nodes) {
            if (a.edge_exists(node0, node1)) {
              RC_LOG() << node0.repr() << "-->" << node1.repr() << "; ";
            }
          }
        }
        RC_LOG() << std::endl;
        RC_LOG() << "Circuit (" << c1.n_qubits() << " qubits, " << c1.n_gates()
                 << " gates): " << c1;

        const unit_bimap_t &initial_map = cu.get_initial_map_ref();
        const unit_bimap_t &final_map = cu.get_final_map_ref();
        RC_LOG() << "Initial Map:" << std::endl;
        for (const auto &x : initial_map.left) {
          RC_LOG() << x.first.repr() << " " << x.second.repr() << std::endl;
        }
        RC_LOG() << "Final Map:" << std::endl;
        for (const auto &x : final_map.left) {
          RC_LOG() << x.first.repr() << " " << x.second.repr() << std::endl;
        }
        if (applied) {
          check_correctness(c, cu);
        } else {
          RC_ASSERT(c == c1);
        }
      });
}

bool check_initial_simplification() {
  return rc::check(
      "initial simplification produces equivalent final state",
      [](const Circuit &c) {
        Circuit c1 = c;
        Transforms::simplify_initial(
            Transforms::AllowClassical::No, Transforms::CreateAllQubits::Yes)
            .apply(c1);
        try {
          const auto s = tket_sim::get_statevector(c);
          const auto s1 = tket_sim::get_statevector(c1);
          RC_ASSERT(tket_sim::compare_statevectors_or_unitaries(
              s, s1, tket_sim::MatrixEquivalence::EQUAL_UP_TO_GLOBAL_PHASE));
        } catch (const Unsupported &) {
        } catch (const NotImplemented &) {
        }
        // If tket-sim doesn't recognise a gate, just ignore it.
        // But if a gate is unknown, which exception SHOULD it be:
        // "Unsupported" or "NotImplemented" ?
      });
}

int main() {
  bool ok = true;
  ok = ok && check_n_qubits();
  ok = ok && check_passes();
  ok = ok && check_mapping();
  ok = ok && check_initial_simplification();
  return ok ? 0 : 1;
}
