// Copyright 2019-2024 Cambridge Quantum Computing
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

#include <math.h>

#include <numeric>
#include <optional>

#include "tket/Circuit/CircPool.hpp"
#include "tket/Circuit/CircUtils.hpp"
#include "tket/Circuit/DAGDefs.hpp"
#include "tket/Gate/GatePtr.hpp"
#include "tket/Gate/Rotation.hpp"
#include "tket/OpType/OpType.hpp"
#include "tket/Utils/EigenConfig.hpp"
#include "tket/Utils/HelperFunctions.hpp"

namespace tket {

namespace CircPool {

/* all of these methods are from https://arxiv.org/pdf/quant-ph/9503016.pdf
or
https://algassert.com/circuits/2015/06/05/Constructing-Large-Controlled-Nots.html
*/

typedef std::vector<std::pair<Edge, Vertex>>
    candidate_t;  // each CnX candidate to decompose needs a spare wire to put
                  // some extra controls on

static Circuit lemma72(unsigned control_m);  // rule lemma 7.2
static void lemma73(
    Circuit& circ, const std::pair<Edge, Vertex>& pairy);  // rule lemma 7.3

// `n` = size of incrementer, circuit returned is of size `n+1`
/* this is slightly less efficient than perhaps it could be -- asymptotically it
   is still good. In an ideal world, this would decompose the incrementers
   smarter for the "even" case */
Circuit incrementer_borrow_1_qubit(unsigned n) {
  bool is_odd = n % 2;
  Circuit circ(n + 1);
  if (n < 6) {
    if (n > 4) circ.append_qubits(C4X_normal_decomp(), {0, 1, 2, 3, 4});
    if (n > 3) circ.append_qubits(C3X_normal_decomp(), {0, 1, 2, 3});
    if (n > 2) circ.add_op<unsigned>(OpType::CCX, {0, 1, 2});
    if (n > 1) circ.add_op<unsigned>(OpType::CX, {0, 1});
    if (n > 0) circ.add_op<unsigned>(OpType::X, {0});
    return circ;
  }
  unsigned j;
  unsigned k;
  // j is bottom qubits, k is top qubits
  // k + j = n + 1 (total no. of qbs)
  if (is_odd) {
    /* if the number of bits we are incrementing is odd, we can just split the
    incrementer in 2 and use the `incrementer_borrow_n_qubits` method twice */
    j = (n + 1) / 2, k = (n + 1) / 2;
  } else {
    /* otherwise, we will also have to pull out a cnx */
    j = n / 2 + 1, k = n / 2;
  }

  Circuit top_incrementer = incrementer_borrow_n_qubits(k);
  std::vector<unsigned> top_qbs(2 * k);
  for (unsigned i = 0; i != k; ++i) {
    top_qbs[2 * i] = i + k;  // borrowed qubits
    top_qbs[2 * i + 1] = i;  // qbs we are trying to increment
  }

  Circuit cnx_top;
  std::vector<unsigned> cnx1_qbs;
  if (k == 3) {  // code is unreachable if k<3
    cnx_top = C3X_normal_decomp();
    cnx1_qbs = {0, 1, 2, n};
  } else if (k == 4) {
    cnx_top = C4X_normal_decomp();
    cnx1_qbs = {0, 1, 2, 3, n};
  } else {
    cnx_top = lemma72(k);  // k controls on cnx
    cnx1_qbs.resize(
        2 * k - 2);  // size of replacement using borrowed qbs = 2*k-1
    std::iota(cnx1_qbs.begin(), cnx1_qbs.end(), 0);
    cnx1_qbs.push_back(n);  // target is last qubit
  }

  Circuit bottom_incrementer;
  std::vector<unsigned> bot_qbs;
  if (is_odd) {
    bottom_incrementer = incrementer_borrow_n_qubits(j);
    bot_qbs.resize(2 * j);
    for (unsigned i = 0; i != j; ++i) {
      bot_qbs[2 * i] = i;  // 0,2,4...n-1 //borrowed qubits
      if (i != 0)
        bot_qbs[2 * i + 1] =
            i + j -
            1;  // 3,5...n //other qbs we are actually trying to increment
    }
    bot_qbs[1] = n;  // incremented qubit 0 in incrementer is bottom one
  } else {
    if (j == 4) {  // code is unreachable if j<4
      bottom_incrementer.add_blank_wires(4);
      bottom_incrementer.append_qubits(C3X_normal_decomp(), {0, 1, 2, 3});
      bottom_incrementer.add_op<unsigned>(OpType::CCX, {0, 1, 2});
      bottom_incrementer.add_op<unsigned>(OpType::CX, {0, 1});
      bottom_incrementer.add_op<unsigned>(OpType::X, {0});
      bot_qbs = {n, n - 3, n - 2, n - 1};
    } else if (j == 5) {
      bottom_incrementer.add_blank_wires(5);
      bottom_incrementer.append_qubits(C4X_normal_decomp(), {0, 1, 2, 3, 4});
      bottom_incrementer.append_qubits(C3X_normal_decomp(), {0, 1, 2, 3});
      bottom_incrementer.add_op<unsigned>(OpType::CCX, {0, 1, 2});
      bottom_incrementer.add_op<unsigned>(OpType::CX, {0, 1});
      bottom_incrementer.add_op<unsigned>(OpType::X, {0});
      bot_qbs = {n, n - 4, n - 3, n - 2, n - 1};
    } else {
      // insert peeled-out cnx
      Circuit cnx_bot = lemma72(j - 1);
      std::vector<unsigned> cnx2_qbs(
          2 * j - 3);  // lemma 7.2 uses 2j-3 qubits for a (j-1)-controlled X
      for (unsigned i = 0; i < j - 2; ++i) {
        cnx2_qbs[i] = k + i;
      }
      cnx2_qbs[j - 2] = n;  // replace the last control
      for (unsigned i = 0; i < j - 3; ++i) {
        cnx2_qbs[j + i - 1] = i;
      }
      cnx2_qbs[2 * j - 4] = n - 1;  // the target of the peeled out cnx

      circ.append_qubits(cnx_bot, cnx2_qbs);
      bottom_incrementer = incrementer_borrow_n_qubits(
          j - 1);  // insert incrementer over remaining qubits
      bot_qbs.resize(2 * j - 2);
      for (unsigned i = 0; i != j - 1; ++i) {
        bot_qbs[2 * i] = i;  // 0,2,4...n-1 //borrowed qubits
        if (i != 0)
          bot_qbs[2 * i + 1] =
              i + k -
              1;  // 3,5...n //other qbs we are actually trying to increment
      }
      bot_qbs[1] = n;  // incremented qubit 0 in incrementer is bottom one
    }
  }

  circ.append_qubits(bottom_incrementer, bot_qbs);
  circ.add_op<unsigned>(
      OpType::X,
      {n});  // to convert controlled-incrementer to larger incrementer
  for (unsigned i = k; i != n; ++i) circ.add_op<unsigned>(OpType::CX, {n, i});
  circ.append_qubits(cnx_top, cnx1_qbs);
  if (!is_odd && j > 5) {
    // insert peeled-out cnx
    Circuit cnx_bot = lemma72(j - 1);
    std::vector<unsigned> cnx2_qbs(2 * j - 3);
    for (unsigned i = 0; i < j - 1; ++i) {
      cnx2_qbs[i] = k + i;
    }
    cnx2_qbs[j - 2] = n;  // the last control
    for (unsigned i = 0; i < j - 3; ++i) {
      cnx2_qbs[j + i - 1] = i;
    }
    cnx2_qbs[2 * j - 4] = n - 1;  // the target of the peeled out cnx
    circ.append_qubits(cnx_bot, cnx2_qbs);
  }
  circ.append_qubits(bottom_incrementer, bot_qbs);
  circ.add_op<unsigned>(OpType::X, {n});
  circ.append_qubits(cnx_top, cnx1_qbs);
  for (unsigned i = k; i != n; ++i) circ.add_op<unsigned>(OpType::CX, {n, i});
  circ.append_qubits(top_incrementer, top_qbs);
  return circ;
}

// an optimised version of
// https://algassert.com/circuits/2015/06/12/Constructing-Large-Increment-Gates.html
/* every second qubit (0,2,4...) is a borrowed qubit */
Circuit incrementer_borrow_n_qubits(unsigned n) {
  const unsigned N = 2 * n;
  Circuit circ(N);
  /* deal with small cases where borrowing qubits is unnecessary */
  if (n < 6) {
    if (n > 4) circ.append_qubits(C4X_normal_decomp(), {1, 3, 5, 7, 9});
    if (n > 3) circ.append_qubits(C3X_normal_decomp(), {1, 3, 5, 7});
    if (n > 2) circ.add_op<unsigned>(OpType::CCX, {1, 3, 5});
    if (n > 1) circ.add_op<unsigned>(OpType::CX, {1, 3});
    if (n > 0) circ.add_op<unsigned>(OpType::X, {1});
    return circ;
  }

  for (unsigned i = 1; i < N; ++i) {
    if (i % 2)
      circ.add_op<unsigned>(OpType::CX, {0, i});
    else
      circ.add_op<unsigned>(OpType::X, {i});
  }

  circ.add_op<unsigned>(OpType::X, {N - 1});

  for (unsigned i = 2; i < N; ++(++i)) {
    std::vector<unsigned> ladder_down_qbs = {i - 2, i - 1, i};
    circ.append_qubits(ladder_down(), ladder_down_qbs);
  }
  circ.add_op<unsigned>(OpType::CX, {N - 2, N - 1});
  for (unsigned i = N - 2; i > 1; --(--i)) {
    std::vector<unsigned> tof_up_qbs = {i - 2, i - 1, i};
    circ.add_op<unsigned>(OpType::CCX, tof_up_qbs);
  }

  for (unsigned i = 2; i < N; ++(++i)) {
    std::vector<unsigned> ladder_down_2_qbs = {i - 2, i - 1, i};
    circ.append_qubits(ladder_down_2(), ladder_down_2_qbs);
  }
  circ.add_op<unsigned>(OpType::CX, {N - 2, N - 1});
  for (unsigned i = N - 2; i > 1; --(--i)) {
    std::vector<unsigned> ladder_up_qbs = {i - 2, i - 1, i};
    circ.append_qubits(ladder_up(), ladder_up_qbs);
  }
  for (unsigned i = 1; i < N; ++(++i))
    circ.add_op<unsigned>(OpType::CX, {0, i});
  return circ;
}

// decompose CnX gate using
// https://algassert.com/circuits/2015/06/22/Using-Quantum-Gates-instead-of-Ancilla-Bits.html
// `n` = no. of controls
Circuit CnX_normal_decomp(unsigned n) {
  /* handle low qb edge cases */
  bool insert_c4xs;  // dictate whether to add C4Xs or n>4 controlled Xs
                     // when bootstrapping
  switch (n) {
    case 0: {
      return X();
    }
    case 1: {
      return CX();
    }
    case 2: {
      return CCX_normal_decomp();
    }
    case 3: {
      return C3X_normal_decomp();
    }
    case 4: {
      return C4X_normal_decomp();
    }
    case 5: {
      insert_c4xs = true;
      break;
    }
    default: {
      insert_c4xs = false;
      break;
    }
  }

  Circuit circ(n + 1);
  std::vector<unsigned> cnx_qbs(n - 1);
  std::iota(cnx_qbs.begin(), cnx_qbs.end(), 0);
  cnx_qbs.push_back(n);
  // first, bootstrap an ancilla qubit
  circ.add_op<unsigned>(OpType::H, {n});
  Vertex cnx1;
  if (insert_c4xs)
    circ.append_qubits(C4X_normal_decomp(), {cnx_qbs});
  else {
    cnx1 = circ.add_op<unsigned>(
        OpType::CnX,
        {cnx_qbs}); /* the CnXs can be decomposed using lemma 7.3 */
  }
  circ.add_op<unsigned>(OpType::Tdg, {n});
  Vertex cx = circ.add_op<unsigned>(OpType::CX, {n - 1, n});
  if (!insert_c4xs) {
    Edge e1 = circ.get_nth_in_edge(cx, 0);
    lemma73(circ, {e1, cnx1});  // replace first CnX using lemma 7.3
  }
  circ.add_op<unsigned>(OpType::T, {n});
  Vertex cnx2;
  if (insert_c4xs)
    circ.append_qubits(C4X_normal_decomp(), {cnx_qbs});
  else {
    cnx2 = circ.add_op<unsigned>(OpType::CnX, {cnx_qbs});
  }
  circ.add_op<unsigned>(OpType::Tdg, {n});
  cx = circ.add_op<unsigned>(OpType::CX, {n - 1, n});
  Edge e2 = circ.get_nth_in_edge(cx, 0);
  if (!insert_c4xs) lemma73(circ, {e2, cnx2});
  circ.add_op<unsigned>(OpType::T, {n});
  circ.add_op<unsigned>(OpType::H, {n});

  // add incremented shift pattern
  Circuit incrementer = incrementer_borrow_1_qubit(n);
  circ.append(incrementer);

  // z rotation layer #1
  std::vector<Op_ptr> z_rots(n);
  double angle = -0.25;
  for (unsigned i = 0; i < n - 1; ++i) {
    unsigned m = n - 1 - i;
    z_rots[i] = get_op_ptr(OpType::Rz, angle);
    circ.add_op<unsigned>(z_rots[i], {m});
    angle /= 2;
  }

  // decremented shift pattern
  for (unsigned i = 0; i < n; ++i) {
    circ.add_op<unsigned>(OpType::X, {i});
  }
  circ.append(incrementer);
  for (unsigned i = 0; i < n; ++i) {
    circ.add_op<unsigned>(OpType::X, {i});
  }

  // z rotation layer #2
  for (unsigned i = 0; i < n - 1; ++i) {
    unsigned m = n - 1 - i;
    Expr ang = z_rots[i]->get_params()[0];
    circ.add_op<unsigned>(get_op_ptr(OpType::Rz, -ang), {m});
  }
  Expr ang = z_rots[n - 2]->get_params()[0];
  circ.add_op<unsigned>(get_op_ptr(OpType::Rz, -ang), {0});

  const Op_ptr ccx = get_op_ptr(OpType::CCX);
  circ.substitute_all(CCX_normal_decomp(), ccx);

  circ.add_phase(std::pow(0.5, n + 1));
  return circ;
}

// unsigned -> which column the change is in
static unsigned find_first_differing_val(
    const std::deque<bool>& d1, const std::deque<bool>& d2) {
  unsigned N = d1.size();
  if (N != d2.size())
    throw ControlDecompError(
        "Error in `find_first_differing_val`: Deques are of differing "
        "sizes");
  for (unsigned i = 0; i < N; ++i) {
    if (d1[i] != d2[i]) return i;
  }
  throw ControlDecompError(
      "Error in `find_first_differing_val`: No change between deques");
}

// optimal decomposition of CnRy and CnZ for 2 < n < 8 according to 1995
// paper... can do better with ZH calculus?
static Circuit lemma71(
    unsigned arity, const Circuit& v_rep, const Circuit& v_dg_rep) {
  unsigned m_controls = arity - 1;
  if (m_controls < 2)
    throw Unsupported(
        "No point using Lemma 7.1 to decompose a gate with less than 2 "
        "controls");

  GrayCode gc = gen_graycode(m_controls);

  Circuit rep(arity);

  unsigned control_qb = 0;
  unsigned last = 0;
  rep.append_with_map(
      v_rep, {{Qubit(0), Qubit(0)}, {Qubit(1), Qubit(m_controls)}});
  // we ignore the 0...0 term, and the first one is always trivial
  // so start from 2
  for (unsigned i = 2; i < gc.size(); ++i) {
    unsigned change = find_first_differing_val(gc[i], gc[i - 1]);
    for (unsigned j = 1; j < gc[i].size(); ++j) {
      if (gc[i][j] == 1) last = j;
    }
    if (change < control_qb) {
      rep.add_op<unsigned>(OpType::CX, {change, control_qb});
    }
    if (change == control_qb) {
      throw ControlDecompError("Error in graycode iteration");
    }
    if (change > control_qb) {
      rep.add_op<unsigned>(OpType::CX, {control_qb, change});
    }

    if ((i % 2) == 0) {
      rep.append_with_map(
          v_dg_rep, {{Qubit(0), Qubit(last)}, {Qubit(1), Qubit(m_controls)}});
    } else
      rep.append_with_map(
          v_rep, {{Qubit(0), Qubit(last)}, {Qubit(1), Qubit(m_controls)}});
    control_qb = last;
  }
  auto [vit, vend] = boost::vertices(rep.dag);
  VertexSet bin;
  for (auto next = vit; vit != vend; vit = next) {
    ++next;
    Vertex v = *vit;
    if (!bin.contains(v)) {
      Op_ptr op = rep.get_Op_ptr_from_Vertex(v);
      OpType optype = op->get_type();
      if (is_multi_qubit_type(optype) && optype != OpType::CX) {
        Circuit replacement = with_CX(as_gate_ptr(op));
        rep.substitute(replacement, v, Circuit::VertexDeletion::No);
        bin.insert(v);
      }
    }
  }
  rep.remove_vertices(
      bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
  return rep;
}

//`control_m` = number of controls
static Circuit lemma72(unsigned control_m) {
  if (control_m < 3)
    throw Unsupported(
        "Cannot decompose a gate with " + std::to_string(control_m) +
        " controls using Lemma 7.2");
  const unsigned n = 2 * control_m - 1;

  Circuit ccx_circ(n);
  unsigned diff = n - control_m;
  for (unsigned i = control_m - 1; i > 1; --i) {
    ccx_circ.add_op<unsigned>(OpType::CCX, {i, i + diff - 1, i + diff});
  }
  ccx_circ.add_op<unsigned>(OpType::CCX, {0, 1, control_m});
  for (unsigned i = 2; i < control_m; ++i) {
    ccx_circ.add_op<unsigned>(OpType::CCX, {i, i + diff - 1, i + diff});
  }
  for (unsigned i = control_m - 2; i > 1; --i) {
    ccx_circ.add_op<unsigned>(OpType::CCX, {i, i + diff - 1, i + diff});
  }
  ccx_circ.add_op<unsigned>(OpType::CCX, {0, 1, control_m});
  for (unsigned i = 2; i < control_m - 1; ++i) {
    ccx_circ.add_op<unsigned>(OpType::CCX, {i, i + diff - 1, i + diff});
  }
  if (ccx_circ.count_gates(OpType::CCX) != 4 * (control_m - 2))
    throw ControlDecompError("Error in Lemma 7.2: CCX gate count is incorrect");
  return ccx_circ;
}

/* this is specifically for performing corollary 7.4 via lemma 7.3 & lemma 7.2 -
 * optimal decomp of a CnX gate */
/* for corollary 7.4, n >= 7 */
// this is a decomposition of a CnX gate using one dirty ancilla
static void lemma73(Circuit& circ, const std::pair<Edge, Vertex>& pairy) {
  Vertex original_cnx = pairy.second;
  Edge original_spare_edge = pairy.first;
  EdgeVec in_edges = circ.get_in_edges(original_cnx);
  const unsigned N =
      in_edges.size() + 1;  // number of qubits in new circuit to replace
  if (N < 5)
    throw Unsupported(
        "Lemma 7.3 cannot decompose CnX with n = " + std::to_string(N - 1));

  EdgeVec out_edges = circ.get_all_out_edges(original_cnx);

  in_edges.insert(in_edges.begin() + in_edges.size() - 1, original_spare_edge);
  out_edges.insert(
      out_edges.begin() + out_edges.size() - 1, original_spare_edge);

  Subcircuit to_delete{in_edges, out_edges, {original_cnx}};
  bool odd_N = N % 2;
  unsigned m1 = (N + 1) / 2;  // number of controls on the first type of CnX
  unsigned m2 = N - m1 - 1;   // number of controls on the second type of CnX

  // make new circuit to substitute later
  Circuit new_circ(N);
  const Op_ptr cnx_op1 = get_op_ptr(OpType::CnX, std::vector<Expr>{}, m1 + 1);
  const Op_ptr cnx_op2 = get_op_ptr(OpType::CnX, std::vector<Expr>{}, m2 + 1);
  std::vector<unsigned> qbs_m1(m1 + 1);
  std::iota(qbs_m1.begin(), qbs_m1.begin() + m1, 0);
  qbs_m1[m1] = N - 1;
  std::vector<unsigned> qbs_m2(m2 + 1);
  std::iota(qbs_m2.begin(), qbs_m2.end(), N - (m2 + 1));

  // add ladder of CnXs to the correct qubits
  Vertex a = new_circ.add_op<unsigned>(cnx_op1, qbs_m1);
  Vertex b = new_circ.add_op<unsigned>(cnx_op2, qbs_m2);
  Vertex c = new_circ.add_op<unsigned>(cnx_op1, qbs_m1);
  Vertex d = new_circ.add_op<unsigned>(cnx_op2, qbs_m2);

  /* first, replace vertex a, putting circuit at back of new_circ */
  unsigned cutsize = odd_N ? N : N - 1;
  EdgeVec cut1(cutsize);
  VertexVec out_verts = new_circ.q_outputs();
  if (odd_N) {
    for (unsigned i = 0; i < N - 2; ++i)
      cut1[i] = new_circ.get_nth_in_edge(out_verts[i], 0);
    cut1[N - 2] = new_circ.get_nth_in_edge(out_verts[N - 1], 0);
    cut1[N - 1] = new_circ.get_nth_in_edge(out_verts[N - 2], 0);
  } else {
    for (unsigned i = 0; i < cutsize; ++i)
      cut1[i] = new_circ.get_nth_in_edge(out_verts[i], 0);
  }

  Circuit a_replacement;
  if (m1 == 1) {
    a_replacement = CX();
  } else if (m1 == 2) {
    a_replacement = CCX();
  } else
    a_replacement = lemma72(m1);  // returns circuit of size 2*m1 - 1
  new_circ.cut_insert(a_replacement, cut1);
  new_circ.remove_vertex(
      a, Circuit::GraphRewiring::Yes, Circuit::VertexDeletion::Yes);

  VertexSet normal_decomp_vertices;

  Circuit b_replacement;
  if (m2 == 1) {
    b_replacement = CX();
  } else if (m2 == 2) {
    b_replacement = CCX();
  } else
    b_replacement = lemma72(m2);  // returns circuit of size 2*m2 - 1

  unsigned b_qubits = b_replacement.n_qubits();

  // reassign cut to back of circuit
  EdgeVec frontier(N);
  for (unsigned i = 0; i < N; ++i)
    frontier[i] = new_circ.get_nth_in_edge(out_verts[i], 0);

  EdgeVec cut2(b_qubits);
  for (unsigned i = N - m2 - 1; i < N - 1; ++i)
    cut2[i - (N - m2 - 1)] =
        frontier[i];  // N-1 - (N-m2-1) = m2 (all the controls)
  for (unsigned i = 0; i < b_qubits - (m2 + 1); ++i)
    cut2[i + m2] = frontier[i];          // empty wires
  cut2[b_qubits - 1] = frontier[N - 1];  // target

  new_circ.cut_insert(b_replacement, cut2);
  new_circ.remove_vertex(
      b, Circuit::GraphRewiring::Yes, Circuit::VertexDeletion::Yes);
  Vertex last_path = out_verts[N - 1];
  Edge init_back_edge = new_circ.get_nth_in_edge(last_path, 0);
  Vertex init_back_path = new_circ.source(init_back_edge);
  normal_decomp_vertices.insert(init_back_path);
  init_back_edge = new_circ.get_last_edge(init_back_path, init_back_edge);
  init_back_path = new_circ.source(init_back_edge);
  OpType backstop = new_circ.get_OpType_from_Vertex(init_back_path);
  while ((backstop != OpType::CCX) && !is_initial_q_type(backstop)) {
    init_back_edge = new_circ.get_last_edge(init_back_path, init_back_edge);
    init_back_path = new_circ.source(init_back_edge);
    backstop = new_circ.get_OpType_from_Vertex(init_back_path);
  }
  normal_decomp_vertices.insert(init_back_path);
  /* now, replace vertex c */
  EdgeVec cut3(cutsize);
  if (odd_N) {
    for (unsigned i = 0; i < N - 2; ++i)
      cut3[i] = new_circ.get_nth_in_edge(out_verts[i], 0);
    cut3[N - 2] = new_circ.get_nth_in_edge(out_verts[N - 1], 0);
    cut3[N - 1] = new_circ.get_nth_in_edge(out_verts[N - 2], 0);
  } else
    for (unsigned i = 0; i < cutsize; ++i)
      cut3[i] = new_circ.get_nth_in_edge(out_verts[i], 0);

  new_circ.cut_insert(a_replacement, cut3);
  new_circ.remove_vertex(
      c, Circuit::GraphRewiring::Yes, Circuit::VertexDeletion::Yes);

  /* now, replace vertex d */
  for (unsigned i = 0; i < N; ++i)
    frontier[i] = new_circ.get_nth_in_edge(out_verts[i], 0);
  EdgeVec cut4(b_qubits);
  for (unsigned i = N - m2 - 1; i < N - 1; ++i)
    cut4[i - (N - m2 - 1)] =
        frontier[i];  // N-1 - (N-m2-1) = m2 (all the controls)
  for (unsigned i = 0; i < b_qubits - (m2 + 1); ++i)
    cut4[i + m2] = frontier[i];          // empty wires
  cut4[b_qubits - 1] = frontier[N - 1];  // target

  new_circ.cut_insert(b_replacement, cut4);
  new_circ.remove_vertex(
      d, Circuit::GraphRewiring::Yes, Circuit::VertexDeletion::Yes);

  last_path = out_verts[N - 1];
  init_back_edge = new_circ.get_nth_in_edge(last_path, 0);
  init_back_path = new_circ.source(init_back_edge);
  normal_decomp_vertices.insert(init_back_path);
  init_back_edge = new_circ.get_last_edge(init_back_path, init_back_edge);
  init_back_path = new_circ.source(init_back_edge);
  backstop = new_circ.get_OpType_from_Vertex(init_back_path);
  while ((backstop != OpType::CCX) && !is_initial_q_type(backstop)) {
    init_back_edge = new_circ.get_last_edge(init_back_path, init_back_edge);
    init_back_path = new_circ.source(init_back_edge);
    backstop = new_circ.get_OpType_from_Vertex(init_back_path);
  }
  normal_decomp_vertices.insert(init_back_path);

  if (m1 > 2 && m2 > 2 && (new_circ.count_gates(OpType::CCX) != 8 * N - 40))
    throw ControlDecompError("Error in Lemma 7.3: CCX gate count is incorrect");
  /* now, replace CCXs with either CX circuit modulo phase shift or just typical
   * Toffoli decomposition */

  auto [vit, vend] = boost::vertices(new_circ.dag);
  for (auto next = vit; vit != vend; vit = next) {
    ++next;
    if (new_circ.get_OpType_from_Vertex(*vit) == OpType::CCX) {
      Subcircuit sub{
          new_circ.get_in_edges(*vit),
          new_circ.get_all_out_edges(*vit),
          {*vit}};
      if (normal_decomp_vertices.find(*vit) == normal_decomp_vertices.end()) {
        new_circ.substitute(
            CCX_modulo_phase_shift(), sub, Circuit::VertexDeletion::Yes);
      } else {
        new_circ.substitute(
            CCX_normal_decomp(), sub, Circuit::VertexDeletion::Yes);
      }
    }
  }
  if (m1 > 2 && m2 > 2 && new_circ.count_gates(OpType::CX) != 24 * N - 108)
    throw ControlDecompError("Error in Lemma 7.3: CX gate count is incorrect");

  circ.substitute(new_circ, to_delete, Circuit::VertexDeletion::Yes);
}

// N must be >= 3
// Assume the unitary is not identity
// Linear decomposition for multi-controlled special unitaries
static void lemma79(
    Circuit& replacement, unsigned N, const Expr& alpha, const Expr& theta,
    const Expr& beta, candidate_t& CCX_candidates) {
  replacement.add_blank_wires(N);

  // Add Controlled C
  if (!equiv_0(beta - alpha, 8)) {
    replacement.add_op<unsigned>(
        OpType::CRz, {(beta - alpha) / 2.}, {N - 2, N - 1});
  }
  // Add the first CnX
  std::vector<unsigned> cnx_qbs(N - 1);
  std::iota(cnx_qbs.begin(), --cnx_qbs.end(), 0);
  cnx_qbs[N - 2] = N - 1;
  const Op_ptr cnx = get_op_ptr(OpType::CnX, std::vector<Expr>{}, N - 1);
  Vertex firstCnX = replacement.add_op<unsigned>(cnx, cnx_qbs);

  // Add Controlled B
  VertexVec vBs;
  if (!equiv_0(alpha + beta, 8)) {
    Vertex vB1 = replacement.add_op<unsigned>(
        OpType::CRz, {(-alpha - beta) / 2.}, {N - 2, N - 1});
    vBs.push_back(vB1);
  }
  if (!equiv_0(theta, 8)) {
    Vertex vB2 = replacement.add_op<unsigned>(
        OpType::CRy, {-theta / 2.}, {N - 2, N - 1});
    vBs.push_back(vB2);
  }
  // At least one of vB1 and vB2 should be set, otherwise it implies that the
  // SU(2) is identity
  if (vBs.empty()) {
    throw ControlDecompError(
        "Unknown error in Lemma 7.9: identity not rejected");
  }

  // Add the second CnX
  Vertex secondCnX = replacement.add_op<unsigned>(cnx, cnx_qbs);

  // Add Controlled A
  if (!equiv_0(theta, 8)) {
    replacement.add_op<unsigned>(OpType::CRy, {theta / 2.}, {N - 2, N - 1});
  }
  if (!equiv_0(alpha, 4)) {
    replacement.add_op<unsigned>(OpType::CRz, {alpha}, {N - 2, N - 1});
  }
  Edge first_e = replacement.get_nth_in_edge(vBs[0], 0);
  Edge second_e = replacement.get_nth_out_edge(vBs.back(), 0);
  CCX_candidates.push_back({first_e, firstCnX});
  CCX_candidates.push_back({second_e, secondCnX});
}

static Circuit CU_to_CU3(const Eigen::Matrix2cd& u) {
  Circuit c(2);
  std::vector<double> tk1_angles = tk1_angles_from_unitary(u);
  Expr theta = tk1_angles[1];
  Expr phi = tk1_angles[0] - 0.5;
  Expr lambda = tk1_angles[2] + 0.5;
  Expr t = tk1_angles[3] - 0.5 * (tk1_angles[0] + tk1_angles[2]);
  c.add_op<unsigned>(OpType::U1, t, {0});
  c.add_op<unsigned>(OpType::CU3, {theta, phi, lambda}, {0, 1});
  c.remove_noops();
  return c;
}

Circuit CnU_gray_code_decomp(unsigned n, const Eigen::Matrix2cd& u) {
  if (n == 0) {
    // Synthesise U using tk1 and phase
    Circuit cnu_circ(1);
    std::vector<double> tk1_angles = tk1_angles_from_unitary(u);
    cnu_circ.add_op<unsigned>(
        OpType::TK1, {tk1_angles[0], tk1_angles[1], tk1_angles[2]}, {0});
    cnu_circ.add_phase(tk1_angles[3]);
    return cnu_circ;
  }
  if (n == 1) {
    return CU_to_CU3(u);
  }

  Eigen::Matrix2cd v_matrix = nth_root(u, 1ULL << (n - 1));
  Eigen::Matrix2cd v_matrix_dag = v_matrix.adjoint();
  Circuit v_rep = CU_to_CU3(v_matrix);
  Circuit v_dg_rep = CU_to_CU3(v_matrix_dag);
  return lemma71(n + 1, v_rep, v_dg_rep);
}

Circuit CnU_gray_code_decomp(unsigned n, const Gate_ptr& gate) {
  static const std::map<OpType, OpType> cu_type_map = {
      {OpType::Rx, OpType::CRx},
      {OpType::Ry, OpType::CRy},
      {OpType::Rz, OpType::CRz},
      {OpType::U1, OpType::CU1}};

  auto it = cu_type_map.find(gate->get_type());
  if (it == cu_type_map.end()) {
    throw Unsupported(
        "The implementation currently only supports Rx, Ry, Rz, U1");
  }

  const OpType& cu_type = it->second;

  if (n == 0) {
    Circuit cnu_circ(1);
    cnu_circ.add_op<unsigned>(gate, {0});
    return cnu_circ;
  }

  Expr angle = gate->get_params()[0];
  if (n == 1) {
    Circuit cnu_circ(2);
    cnu_circ.add_op<unsigned>(cu_type, angle, {0, 1});
    return cnu_circ;
  }
  Expr param;
  std::optional<double> reduced = eval_expr_mod(angle, 4);
  if (reduced)
    param = reduced.value();
  else
    param = angle;

  param = param / (1ULL << (n - 1));
  Circuit v_rep(2);
  Circuit v_dg_rep(2);
  v_rep.add_op<unsigned>(cu_type, param, {0, 1});
  v_dg_rep.add_op<unsigned>(cu_type, -param, {0, 1});
  return lemma71(n + 1, v_rep, v_dg_rep);
}

Circuit CnRy_normal_decomp(const Op_ptr op, unsigned arity) {
  if (op->get_type() != OpType::CnRy) {
    throw CircuitInvalidity("Operation not CnRy");
  }
  OpDesc desc = op->get_desc();
  Expr angle = op->get_params()[0];
  Circuit rep;
  switch (arity) {
    case 0: {
      throw CircuitInvalidity("Circuit has a CnRy with no in edges!");
    }
    case 1: {
      rep.add_blank_wires(1);
      rep.add_op<unsigned>(OpType::Ry, angle, {0});
      break;
    }
    case 2: {
      rep = CircPool::CRy_using_CX(angle);
      break;
    }
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
    case 8: {
      rep = CnU_gray_code_decomp(
          arity - 1, as_gate_ptr(get_op_ptr(OpType::Ry, angle)));
      break;
    }
    default: {
      rep = CnSU2_linear_decomp(arity - 1, 0., angle, 0.);
      auto [x, xend] = boost::vertices(rep.dag);
      for (auto xnext = x; x != xend; x = xnext) {
        ++xnext;
        OpType type = rep.get_OpType_from_Vertex(*x);
        if (type == OpType::CRy) {
          Expr x_angle = rep.get_Op_ptr_from_Vertex(*x)->get_params()[0];
          Circuit new_circ = CircPool::CRy_using_CX(x_angle);
          Subcircuit sub{rep.get_in_edges(*x), rep.get_all_out_edges(*x), {*x}};
          rep.substitute(new_circ, sub, Circuit::VertexDeletion::Yes);
        } else if (type == OpType::CRz) {
          throw ControlDecompError(
              "Error in Lemma 7.9: Y rotation isn't recognized");
        }
      }
      break;
    }
  }
  return rep;
}

// decompose CnX gate using lemma 7.1
// `n` = no. of controls
Circuit CnX_gray_decomp(unsigned n) {
  switch (n) {
    case 0: {
      return X();
    }
    case 1: {
      return CX();
    }
    case 2: {
      return CCX_normal_decomp();
    }
    case 3: {
      return C3X_normal_decomp();
    }
    case 4: {
      return C4X_normal_decomp();
    }
    default: {
      Circuit circ(n + 1);
      circ.add_op<unsigned>(OpType::H, {n});
      circ.append(
          CnU_gray_code_decomp(n, as_gate_ptr(get_op_ptr(OpType::U1, 1.0))));
      circ.add_op<unsigned>(OpType::H, {n});
      return circ;
    }
  }
}

static void add_cu_using_cu3(
    const unsigned& ctrl, const unsigned& trgt, Circuit& circ,
    const Eigen::Matrix2cd& u) {
  unit_map_t unit_map;
  unit_map.insert({Qubit(0), Qubit(ctrl)});
  unit_map.insert({Qubit(1), Qubit(trgt)});
  Circuit cnu_circ = CU_to_CU3(u);
  circ.append_with_map(cnu_circ, unit_map);
}

// Add pn to qubits {1,...,n}, assume n > 1
static void add_pn(Circuit& circ, unsigned n, bool inverse) {
  TKET_ASSERT(n > 1);
  // pn is self commute
  for (unsigned i = 2; i < n + 1; i++) {
    unsigned long long d = 1ULL << (n - i + 1);
    circ.add_op<unsigned>(OpType::CRx, (inverse ? -1. : 1.) / d, {i - 1, n});
  }
}

// Add pn(u) to qubits {1,...,n}, assume n > 1
static void add_pn_unitary(
    Circuit& circ, const Eigen::Matrix2cd& u, unsigned n, bool inverse) {
  TKET_ASSERT(n > 1);
  // pn_(u) is self commute
  for (unsigned i = 2; i < n + 1; i++) {
    Eigen::Matrix2cd m = nth_root(u, 1ULL << (n - i + 1));
    if (inverse) m.adjointInPlace();
    add_cu_using_cu3(i - 1, n, circ, m);
  }
}

// Add an incrementer without toggling the least significant bit
// Apply to qubits {0,...,n-1}, assume n > 1
static void add_qn(Circuit& circ, unsigned n) {
  TKET_ASSERT(n > 1);
  for (unsigned i = n - 1; i > 1; i--) {
    unsigned long long d = 1ULL << (i - 1);
    add_pn(circ, i, false);
    circ.add_op<unsigned>(OpType::CRx, (double)1 / d, {0, i});
  }
  circ.add_op<unsigned>(OpType::CRx, 1, {0, 1});
  for (unsigned i = 2; i < n; i++) {
    add_pn(circ, i, true);
  }
}

// https://arxiv.org/abs/2203.11882 Equation 5
Circuit incrementer_linear_depth(unsigned n, bool lsb) {
  if (n == 0) {
    return Circuit();
  }
  Circuit circ(n);
  if (n > 1) {
    add_qn(circ, n);
  }
  if (lsb) {
    // Some optimisations might have better handlings for X gates
    // so use X instead of Rx(1)
    circ.add_op<unsigned>(OpType::X, {0});
    circ.add_phase(-0.5);
  }
  return circ;
}

// https://arxiv.org/abs/2203.11882 Equation 3
Circuit CnU_linear_depth_decomp(unsigned n, const Eigen::Matrix2cd& u) {
  if (!is_unitary(u)) {
    throw CircuitInvalidity(
        "Matrix for the controlled operation must be unitary");
  }
  Circuit circ(n + 1);

  if (n == 0) {
    // Synthesis U using tk1 and phase
    std::vector<double> tk1_angles = tk1_angles_from_unitary(u);
    circ.add_op<unsigned>(
        OpType::TK1, {tk1_angles[0], tk1_angles[1], tk1_angles[2]}, {0});
    circ.add_phase(tk1_angles[3]);
    return circ;
  }
  if (n == 1) {
    add_cu_using_cu3(0, 1, circ, u);
    return circ;
  }

  // Add pn(u) to qubits {1,...,n}
  add_pn_unitary(circ, u, n, false);
  // Add CU to {0, n}
  Eigen::Matrix2cd m = nth_root(u, 1ULL << (n - 1));
  add_cu_using_cu3(0, n, circ, m);
  // Add incrementer (without toggling q0) to {0,...,n-1}
  Circuit qn = incrementer_linear_depth(n, false);
  Circuit qn_dag = qn.dagger();
  circ.append(qn);

  // Add pn(u).dagger to qubits {1,...,n}
  add_pn_unitary(circ, u, n, true);

  // Add incrementer inverse (without toggling q0) to {0,...,n-1}
  circ.append(qn_dag);

  return circ;
}

Circuit CnSU2_linear_decomp(
    unsigned n, const Expr& alpha, const Expr& theta, const Expr& beta) {
  // W == I iff one of the following two conditions is met
  // 1. t/2 is even, and (a + b)/2 is even
  // 2. t/2 is odd, and (a + b)/2 is odd
  // check if SU(2) is identity
  if ((equiv_0(theta / 2., 2) && equiv_0((alpha + beta) / 2., 2)) ||
      (equiv_val(theta / 2., 1., 2) && equiv_val((alpha + beta) / 2., 1., 2))) {
    return Circuit(n + 1);
  }

  Circuit circ;

  if (n == 0) {
    // add tk1
    circ.add_blank_wires(1);
    circ.add_op<unsigned>(OpType::TK1, {alpha + 0.5, theta, beta - 0.5}, {0});
    return circ;
  }

  // SU(2) matrix W expressed as Rz(a)Ry(t)Rz(b)
  Expr a = alpha;
  Expr b = beta;
  Expr t = theta;

  // Lemma 4.3: W = A*X*B*X*C
  // By lemma 5.4, C is identity iff W can be expressed as Rz(a')Ry(t')Rz(a')
  // We handle the following two cases
  // if (a-b)/2 is even, a' = (a + b)/2, t' = t
  // if (a-b)/2 is odd, a' = (a + b)/2, t' = - t
  if (equiv_0((a - b) / 2., 2)) {
    a = (a + b) / 2.;
    b = a;
  } else if (equiv_val((a - b) / 2., 1., 2)) {
    a = (a + b) / 2.;
    b = a;
    t = -t;
  }

  // We test whether W can be expressed as a single Ry(t'')
  if (equiv_0((a - b) / 2., 2)) {
    if (equiv_val((a + b) / 2., 1., 2)) {
      // if (a-b)/2 is even, and (a+b)/2 is odd, t'' = 2-t
      a = 0.;
      b = 0.;
      t = 2. - t;
    } else if (equiv_0((a + b) / 2., 2)) {
      // if (a-b)/2 is even, and (a+b)/2 is even, t'' = t
      a = 0.;
      b = 0.;
    }
  } else if (equiv_val((a - b) / 2., 1., 2)) {
    if (equiv_val((a + b) / 2., 1., 2)) {
      // if (a-b)/2 is odd, and (a+b)/2 is odd, t'' = 2-t
      a = 0.;
      b = 0.;
      t = 2. + t;
    } else if (equiv_0((a + b) / 2., 2)) {
      // if (a-b)/2 is odd, and (a+b)/2 is even, t'' = -t
      a = 0.;
      b = 0.;
      t = -t;
    }
  }
  if (n == 1) {
    circ.add_blank_wires(2);
    if (!equiv_0(b - a, 8)) {
      circ.add_op<unsigned>(OpType::Rz, {(b - a) / 2.}, {1});
    }
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    if (!equiv_0(a + b, 8)) {
      circ.add_op<unsigned>(OpType::Rz, {(-a - b) / 2.}, {1});
    }
    if (!equiv_0(t, 8)) {
      circ.add_op<unsigned>(OpType::Ry, {-t / 2.}, {1});
    }
    circ.add_op<unsigned>(OpType::CX, {0, 1});
    if (!equiv_0(t, 8)) {
      circ.add_op<unsigned>(OpType::Ry, {t / 2.}, {1});
    }
    if (!equiv_0(a, 4)) {
      circ.add_op<unsigned>(OpType::Rz, {a}, {1});
    }
    return circ;
  }
  // Using lemma 7.9 for n >= 2
  candidate_t ct;
  lemma79(circ, n + 1, a, t, b, ct);
  for (const std::pair<Edge, Vertex>& pairy : ct) {
    Vertex original_cnx = pairy.second;
    unsigned cnx_arity = circ.n_in_edges(original_cnx);
    switch (cnx_arity) {
      case 2: {
        circ.dag[original_cnx] = {get_op_ptr(OpType::CX)};
        break;
      }
      case 3: {
        circ.dag[original_cnx] = {get_op_ptr(OpType::CCX)};
        break;
      }
      default: {
        lemma73(circ, pairy);
      }
    }
  }
  return circ;
}

}  // namespace CircPool

}  // namespace tket
