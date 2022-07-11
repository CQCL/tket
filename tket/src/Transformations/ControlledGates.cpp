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

#include "ControlledGates.hpp"

#include <math.h>

#include <numeric>
#include <optional>

#include "Circuit/CircPool.hpp"
#include "Circuit/DAGDefs.hpp"
#include "OpType/OpType.hpp"
#include "Transform.hpp"
#include "Utils/EigenConfig.hpp"
#include "Utils/HelperFunctions.hpp"

namespace tket {

namespace Transforms {

/* all of these methods are from https://arxiv.org/pdf/quant-ph/9503016.pdf
or
https://algassert.com/circuits/2015/06/05/Constructing-Large-Controlled-Nots.html
*/

typedef std::vector<std::pair<Edge, Vertex>>
    candidate_t;  // each CnX candidate to decompose needs a spare wire to put
                  // some extra controls on

static Circuit lemma54(const Expr& angle);  // refers to rule lemma 5.4 in paper
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
    if (n > 4)
      circ.append_qubits(CircPool::C4X_normal_decomp(), {0, 1, 2, 3, 4});
    if (n > 3) circ.append_qubits(CircPool::C3X_normal_decomp(), {0, 1, 2, 3});
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
    cnx_top = CircPool::C3X_normal_decomp();
    cnx1_qbs = {0, 1, 2, n};
  } else if (k == 4) {
    cnx_top = CircPool::C4X_normal_decomp();
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
      bottom_incrementer.append_qubits(
          CircPool::C3X_normal_decomp(), {0, 1, 2, 3});
      bottom_incrementer.add_op<unsigned>(OpType::CCX, {0, 1, 2});
      bottom_incrementer.add_op<unsigned>(OpType::CX, {0, 1});
      bottom_incrementer.add_op<unsigned>(OpType::X, {0});
      bot_qbs = {n, n - 3, n - 2, n - 1};
    } else if (j == 5) {
      bottom_incrementer.add_blank_wires(5);
      bottom_incrementer.append_qubits(
          CircPool::C4X_normal_decomp(), {0, 1, 2, 3, 4});
      bottom_incrementer.append_qubits(
          CircPool::C3X_normal_decomp(), {0, 1, 2, 3});
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
    if (n > 4)
      circ.append_qubits(CircPool::C4X_normal_decomp(), {1, 3, 5, 7, 9});
    if (n > 3) circ.append_qubits(CircPool::C3X_normal_decomp(), {1, 3, 5, 7});
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
    circ.append_qubits(CircPool::ladder_down(), ladder_down_qbs);
  }
  circ.add_op<unsigned>(OpType::CX, {N - 2, N - 1});
  for (unsigned i = N - 2; i > 1; --(--i)) {
    std::vector<unsigned> tof_up_qbs = {i - 2, i - 1, i};
    circ.add_op<unsigned>(OpType::CCX, tof_up_qbs);
  }

  for (unsigned i = 2; i < N; ++(++i)) {
    std::vector<unsigned> ladder_down_2_qbs = {i - 2, i - 1, i};
    circ.append_qubits(CircPool::ladder_down_2(), ladder_down_2_qbs);
  }
  circ.add_op<unsigned>(OpType::CX, {N - 2, N - 1});
  for (unsigned i = N - 2; i > 1; --(--i)) {
    std::vector<unsigned> ladder_up_qbs = {i - 2, i - 1, i};
    circ.append_qubits(CircPool::ladder_up(), ladder_up_qbs);
  }
  for (unsigned i = 1; i < N; ++(++i))
    circ.add_op<unsigned>(OpType::CX, {0, i});
  return circ;
}

// decompose CnX gate using
// https://algassert.com/circuits/2015/06/22/Using-Quantum-Gates-instead-of-Ancilla-Bits.html
// `n` = no. of controls
Circuit cnx_normal_decomp(unsigned n) {
  /* handle low qb edge cases */
  bool insert_c4xs;  // dictate whether to add C4Xs or n>4 controlled Xs
                     // when bootstrapping
  switch (n) {
    case 0: {
      return CircPool::X();
    }
    case 1: {
      return CircPool::CX();
    }
    case 2: {
      return CircPool::CCX_normal_decomp();
    }
    case 3: {
      return CircPool::C3X_normal_decomp();
    }
    case 4: {
      return CircPool::C4X_normal_decomp();
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
    circ.append_qubits(CircPool::C4X_normal_decomp(), {cnx_qbs});
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
    circ.append_qubits(CircPool::C4X_normal_decomp(), {cnx_qbs});
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

  decomp_CCX().apply(circ);
  circ.add_phase(std::pow(0.5, n + 1));
  return circ;
}

/* assumes vert is controlled Ry with 1 control */
/* decomposes CRy into 2 CXs and 2 Ry gates */
static Circuit lemma54(const Expr& angle) {
  Circuit new_circ(2);
  std::vector<Expr> A_params = {angle / 2.};
  std::vector<Expr> B_params = {-angle / 2.};
  const Op_ptr A = get_op_ptr(OpType::Ry, A_params);
  const Op_ptr B = get_op_ptr(OpType::Ry, B_params);
  new_circ.add_op<unsigned>(A, {1});
  new_circ.add_op<unsigned>(OpType::CX, {0, 1});
  new_circ.add_op<unsigned>(B, {1});
  new_circ.add_op<unsigned>(OpType::CX, {0, 1});
  return new_circ;
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
    unsigned arity, const Expr& angle, const OpType& cr_type) {
  unsigned m_controls = arity - 1;
  if (m_controls < 2)
    throw Unsupported(
        "No point using Lemma 7.1 to decompose a gate with less than 2 "
        "controls");
  if (m_controls > 7)
    throw Unsupported(
        "Using Lemma 7.1 to decompose a gate with more than 7 controls "
        "is inefficient");
  if (cr_type != OpType::CRy && cr_type != OpType::CU1)
    throw Unsupported(
        "The implementation currently only supports CU1 and CRy ");

  GrayCode gc = gen_graycode(m_controls);
  unsigned n_square_roots = m_controls - 1;

  Circuit rep(arity);
  Expr param;
  std::optional<double> reduced = eval_expr_mod(angle, 4);
  if (reduced)
    param = reduced.value();
  else
    param = angle;

  param = param / (1 << (n_square_roots));

  const Op_ptr V_op = get_op_ptr(cr_type, param);
  const Op_ptr V_dg = get_op_ptr(cr_type, -param);

  unsigned control_qb = 0;
  unsigned last = 0;
  rep.add_op<unsigned>(V_op, {0, m_controls});
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
      rep.add_op<unsigned>(V_dg, {last, m_controls});
    } else
      rep.add_op<unsigned>(V_op, {last, m_controls});
    control_qb = last;
  }
  unsigned correct_gate_count =
      ((1 << m_controls) - 1) + ((1 << m_controls) - 2);
  if (rep.n_gates() != correct_gate_count)
    throw ControlDecompError("Error in Lemma 7.1: Gate count is incorrect");
  auto [vit, vend] = boost::vertices(rep.dag);
  VertexSet bin;
  for (auto next = vit; vit != vend; vit = next) {
    ++next;
    Vertex v = *vit;
    if (!bin.contains(v)) {
      OpType optype = rep.get_OpType_from_Vertex(v);
      if (optype == OpType::CRy || optype == OpType::CU1) {
        Expr v_angle = rep.get_Op_ptr_from_Vertex(v)->get_params()[0];
        Circuit replacement = (optype == OpType::CRy)
                                  ? CircPool::CRy_using_CX(v_angle)
                                  : CircPool::CU1_using_CX(v_angle);
        Subcircuit sub{rep.get_in_edges(v), rep.get_all_out_edges(v), {v}};
        rep.substitute(replacement, sub, Circuit::VertexDeletion::No);
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
    a_replacement = CircPool::CX();
  } else if (m1 == 2) {
    a_replacement = CircPool::CCX();
  } else
    a_replacement = lemma72(m1);  // returns circuit of size 2*m1 - 1
  new_circ.cut_insert(a_replacement, cut1);
  new_circ.remove_vertex(
      a, Circuit::GraphRewiring::Yes, Circuit::VertexDeletion::Yes);

  VertexSet normal_decomp_vertices;

  Circuit b_replacement;
  if (m2 == 1) {
    b_replacement = CircPool::CX();
  } else if (m2 == 2) {
    b_replacement = CircPool::CCX();
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
            CircPool::CCX_modulo_phase_shift(), sub,
            Circuit::VertexDeletion::Yes);
      } else {
        new_circ.substitute(
            CircPool::CCX_normal_decomp(), sub, Circuit::VertexDeletion::Yes);
      }
    }
  }
  if (m1 > 2 && m2 > 2 && new_circ.count_gates(OpType::CX) != 24 * N - 108)
    throw ControlDecompError("Error in Lemma 7.3: CX gate count is incorrect");

  circ.substitute(new_circ, to_delete, Circuit::VertexDeletion::Yes);
}

// N must be >= 3
static void lemma79(
    Circuit& replacement, unsigned N, const Expr& angle,
    candidate_t& CCX_candidates) {
  replacement.add_blank_wires(N);

  std::vector<Expr> A_params = {angle / 2.};
  std::vector<Expr> B_params = {-angle / 2.};
  const Op_ptr A = get_op_ptr(OpType::CnRy, A_params, 2);
  const Op_ptr B = get_op_ptr(OpType::CnRy, B_params, 2);

  Vertex vA = replacement.add_op<unsigned>(A, {N - 2, N - 1});  // A
  std::vector<unsigned> cnx_qbs(N - 1);
  std::iota(cnx_qbs.begin(), --cnx_qbs.end(), 0);
  cnx_qbs[N - 2] = N - 1;
  const Op_ptr cnx = get_op_ptr(OpType::CnX, std::vector<Expr>{}, N - 1);
  Vertex firstCnX = replacement.add_op<unsigned>(cnx, cnx_qbs);
  Vertex vB = replacement.add_op<unsigned>(B, {N - 2, N - 1});  // B
  CCX_candidates.push_back(
      {boost::edge(vA, vB, replacement.dag).first, firstCnX});
  Vertex secondCnX = replacement.add_op<unsigned>(cnx, cnx_qbs);
  Edge out_edge_spare = replacement.get_nth_out_edge(vB, 0);
  CCX_candidates.push_back({out_edge_spare, secondCnX});
}

/* naive decomposition - there are cases we can do better if we can eg. ignore
 * phase */
Transform decomp_CCX() {
  return Transform([](Circuit& circ) {
    const Op_ptr ccx = get_op_ptr(OpType::CCX);
    return circ.substitute_all(CircPool::CCX_normal_decomp(), ccx);
  });
}

Circuit decomposed_CnRy(const Op_ptr op, unsigned arity) {
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
      rep = lemma54(angle);
      break;
    }
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
    case 8: {
      rep = lemma71(arity, angle, OpType::CRy);
      break;
    }
    default: {
      candidate_t ct;
      lemma79(rep, arity, angle, ct);
      if (ct.size() != 2)
        throw ControlDecompError(
            "Unknown error in controlled gate decomposition");
      for (const std::pair<Edge, Vertex>& pairy : ct) lemma73(rep, pairy);
      auto [x, xend] = boost::vertices(rep.dag);
      for (auto xnext = x; x != xend; x = xnext) {
        ++xnext;
        OpType type = rep.get_OpType_from_Vertex(*x);
        if (type == OpType::CnRy) {
          Expr x_angle = rep.get_Op_ptr_from_Vertex(*x)->get_params()[0];
          Circuit new_circ = lemma54(x_angle);
          Subcircuit sub{rep.get_in_edges(*x), rep.get_all_out_edges(*x), {*x}};
          rep.substitute(new_circ, sub, Circuit::VertexDeletion::Yes);
        }
      }
      break;
    }
  }
  return rep;
}

Transform decomp_controlled_Rys() {
  return Transform([](Circuit& circ) {
    bool success = decomp_CCX().apply(circ);
    auto [vit, vend] = boost::vertices(circ.dag);
    for (auto next = vit; vit != vend; vit = next) {
      ++next;
      Vertex v = *vit;
      const Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
      unsigned arity = circ.n_in_edges(v);
      if (op->get_type() == OpType::CnRy) {
        success = true;
        Circuit rep = decomposed_CnRy(op, arity);
        EdgeVec inedges = circ.get_in_edges(v);
        Subcircuit final_sub{inedges, circ.get_all_out_edges(v), {v}};
        circ.substitute(rep, final_sub, Circuit::VertexDeletion::Yes);
      }
    }
    return success;
  });
}

Transform decomp_arbitrary_controlled_gates() {
  return decomp_controlled_Rys() >> decomp_CCX();
}

// decompose CnX gate using lemma 7.1
// `n` = no. of controls
Circuit cnx_gray_decomp(unsigned n) {
  switch (n) {
    case 0: {
      return CircPool::X();
    }
    case 1: {
      return CircPool::CX();
    }
    case 2: {
      return CircPool::CCX_normal_decomp();
    }
    case 3: {
      return CircPool::C3X_normal_decomp();
    }
    case 4: {
      return CircPool::C4X_normal_decomp();
    }
    default: {
      Circuit circ(n + 1);
      circ.add_op<unsigned>(OpType::H, {n});
      circ.append(lemma71(n + 1, 1.0, OpType::CU1));
      circ.add_op<unsigned>(OpType::H, {n});
      return circ;
    }
  }
}

}  // namespace Transforms

}  // namespace tket
