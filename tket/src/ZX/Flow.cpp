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

#include "ZX/Flow.hpp"

#include "Utils/GraphHeaders.hpp"
#include "Utils/MatrixAnalysis.hpp"

namespace tket {

namespace zx {

Flow::Flow(
    const std::map<ZXVert, ZXVertSeqSet>& c,
    const std::map<ZXVert, unsigned>& d)
    : c_(c), d_(d) {}

ZXVertSeqSet Flow::c(const ZXVert& v) const { return c_.at(v); }

ZXVertSeqSet Flow::odd(const ZXVert& v, const ZXDiagram& diag) const {
  sequenced_map_t<ZXVert, unsigned> parities;
  ZXVertSeqSet cv = c(v);
  for (const ZXVert& u : cv.get<TagSeq>()) {
    for (const ZXVert& n : diag.neighbours(u)) {
      if (diag.get_zxtype(n) == ZXType::Output) continue;
      sequenced_map_t<ZXVert, unsigned>::iterator found =
          parities.get<TagKey>().find(n);
      if (found == parities.get<TagKey>().end()) {
        parities.insert({n, 1});
      } else {
        parities.replace(found, {n, found->second + 1});
      }
    }
  }
  ZXVertSeqSet odds;
  for (const std::pair<ZXVert, unsigned>& p : parities.get<TagSeq>()) {
    if (p.second % 2 == 1) {
      odds.insert(p.first);
    }
  }
  return odds;
}

unsigned Flow::d(const ZXVert& v) const { return d_.at(v); }

void Flow::verify(const ZXDiagram& diag) const {
  if (!diag.is_MBQC())
    throw ZXError("Verifying a flow for a diagram that is not in MBQC form");
  std::set<ZXVert> output_set;
  for (const ZXVert& o : diag.get_boundary(ZXType::Output)) {
    output_set.insert(diag.neighbours(o).at(0));
  }
  BGL_FORALL_VERTICES(u, *diag.graph, ZXGraph) {
    ZXType type = diag.get_zxtype(u);
    if (is_boundary_type(type) || output_set.find(u) != output_set.end())
      continue;
    ZXVertSeqSet uc = c(u);
    ZXVertSeqSet uodd = odd(u, diag);
    for (const ZXVert& v : uc.get<TagSeq>()) {
      ZXType vt = diag.get_zxtype(v);
      if (u != v && vt != ZXType::PX && vt != ZXType::PY && d(u) <= d(v))
        throw ZXError("A qubit has an X correction in its past");
      if (u != v && vt == ZXType::PY && d(u) <= d(v) &&
          uodd.find(v) == uodd.end())
        throw ZXError("A past Y vertex receives an X correction");
    }
    for (const ZXVert& v : uodd.get<TagSeq>()) {
      ZXType vt = diag.get_zxtype(v);
      if (u != v && vt != ZXType::PY && vt != ZXType::PZ && d(u) <= d(v))
        throw ZXError("A qubit has a Z correction in its past");
      if (u != v && vt == ZXType::PY && d(u) <= d(v) && uc.find(v) == uc.end())
        throw ZXError("A past Y vertex receives a Z correction");
    }
    bool self_x = (uc.find(u) != uc.end());
    bool self_z = (uodd.find(u) != uodd.end());
    switch (type) {
      case ZXType::XY: {
        if (self_x || !self_z)
          throw ZXError("XY vertex must be corrected with a Z");
        break;
      }
      case ZXType::XZ: {
        if (!self_x || !self_z)
          throw ZXError("XZ vertex must be corrected with a Y");
        break;
      }
      case ZXType::YZ: {
        if (!self_x || self_z)
          throw ZXError("YZ vertex must be corrected with an X");
        break;
      }
      case ZXType::PX: {
        if (!self_z) throw ZXError("PX vertex must be corrected with a Y or Z");
        break;
      }
      case ZXType::PY: {
        if (self_x == self_z)
          throw ZXError("PY vertex must be corrected with an X or Z");
        break;
      }
      case ZXType::PZ: {
        if (!self_x)
          throw ZXError("PZ vertex must be corrected with an X or Y");
        break;
      }
      default:
        throw ZXError("Invalid ZXType for MBQC diagram");
    }
  }
}

void Flow::focus(const ZXDiagram& diag) {
  std::set<ZXVert> output_set;
  for (const ZXVert& o : diag.get_boundary(ZXType::Output)) {
    output_set.insert(diag.neighbours(o).at(0));
  }
  std::map<unsigned, ZXVertVec> order;
  for (const std::pair<const ZXVert, unsigned>& p : d_) {
    auto found = order.find(p.second);
    if (found == order.end())
      order.insert({p.second, {p.first}});
    else
      found->second.push_back(p.first);
  }

  for (const std::pair<const unsigned, ZXVertVec>& p : order) {
    for (const ZXVert& u : p.second) {
      if (output_set.find(u) != output_set.end()) continue;
      ZXVertSeqSet uc = c(u);
      ZXVertSeqSet uodd = odd(u, diag);
      sequenced_map_t<ZXVert, unsigned> parities;
      for (const ZXVert& v : uc.get<TagSeq>()) parities.insert({v, 1});
      for (const ZXVert& v : uc.get<TagSeq>()) {
        if (v == u) continue;
        ZXType vtype = diag.get_zxtype(v);
        if ((vtype != ZXType::XY && vtype != ZXType::PX &&
             vtype != ZXType::PY) ||
            (vtype == ZXType::PY && uodd.find(v) == uodd.end())) {
          ZXVertSeqSet cv = c(v);
          for (const ZXVert& w : cv.get<TagSeq>()) {
            auto found = parities.get<TagKey>().find(w);
            if (found == parities.get<TagKey>().end())
              parities.insert({w, 1});
            else
              parities.replace(found, {w, found->second + 1});
          }
        }
      }
      for (const ZXVert& v : uodd.get<TagSeq>()) {
        if (v == u) continue;
        ZXType vtype = diag.get_zxtype(v);
        if ((output_set.find(v) == output_set.end() && vtype != ZXType::XZ &&
             vtype != ZXType::YZ && vtype != ZXType::PY &&
             vtype != ZXType::PZ) ||
            (vtype == ZXType::PY && uc.find(v) == uc.end())) {
          ZXVertSeqSet cv = c(v);
          for (const ZXVert& w : cv.get<TagSeq>()) {
            auto found = parities.get<TagKey>().find(w);
            if (found == parities.get<TagKey>().end())
              parities.insert({w, 1});
            else
              parities.replace(found, {w, found->second + 1});
          }
        }
      }
      ZXVertSeqSet new_c;
      for (const std::pair<const ZXVert, unsigned> p : parities.get<TagSeq>()) {
        if (p.second % 2 == 1) new_c.insert(p.first);
      }
      c_.at(u) = new_c;
    }
  }
}

Flow Flow::identify_causal_flow(const ZXDiagram& diag) {
  // Check diagram has the expected form for causal flow
  if (!diag.is_MBQC())
    throw ZXError("ZXDiagram must be in MBQC form to identify causal flow");
  std::set<ZXVert> input_set;
  for (const ZXVert& i : diag.get_boundary(ZXType::Input)) {
    input_set.insert(diag.neighbours(i).at(0));
  }
  std::set<ZXVert> output_set;
  for (const ZXVert& o : diag.get_boundary(ZXType::Output)) {
    output_set.insert(diag.neighbours(o).at(0));
  }
  BGL_FORALL_VERTICES(v, *diag.graph, ZXGraph) {
    ZXType vtype = diag.get_zxtype(v);
    if (!is_boundary_type(vtype) && output_set.find(v) == output_set.end() &&
        vtype != ZXType::XY)
      throw ZXError(
          "Causal flow is only defined when all measured vertices are XY");
  }

  // solved contains all vertices for which we have found corrections
  ZXVertSeqSet solved;
  // correctors are those vertices that have been solved but are not yet
  // fl.c(u) for some u
  ZXVertSeqSet correctors;
  // past[v] is undefined if v is not yet solved
  // past[v] is the number of neighbours of v that are still unsolved
  // When past[v] drops to 1, we can correct the unsolved vertex using an X on
  // v and Z on all of its other neighbours
  std::map<ZXVert, unsigned> past;
  Flow fl{{}, {}};

  // Outputs are trivially solved
  for (const ZXVert& o : diag.get_boundary(ZXType::Output)) {
    // ZX Diagrams requires each output to have a unique edge to another vertex
    ZXVert n = diag.neighbours(o).at(0);
    past[n] = diag.degree(n) - 1;
    solved.insert(o);
    solved.insert(n);
    fl.c_.insert({n, {}});
    fl.d_.insert({n, 0});
    // Add output to correctors if it is not an input
    if (input_set.find(n) == input_set.end()) correctors.insert(n);
  }

  unsigned depth = 1;

  do {
    ZXVertSeqSet new_correctors;
    for (const ZXVert& v : correctors.get<TagSeq>()) {
      // Determine whether |N(v) cap unsolved| == 1 to find u
      ZXVert u;
      unsigned n_found = 0;
      for (const ZXVert& vn : diag.neighbours(v)) {
        if (solved.find(vn) == solved.end()) {
          u = vn;
          ++n_found;
        }
      }
      if (n_found != 1) continue;

      // Can correct u by firing stabilizer of v
      fl.c_.insert({u, {v}});
      fl.d_.insert({u, depth});
      solved.insert(u);

      // Determine any new correctors
      n_found = 0;
      bool in = false;
      for (const ZXVert& un : diag.neighbours(u)) {
        if (diag.get_zxtype(un) == ZXType::Input) {
          in = true;
          solved.insert(un);
          continue;
        }
        if (solved.find(un) == solved.end()) {
          ++n_found;
        }
        // Another neighbour of un has been solved, so check if it can now
        // correct something
        auto it = past.find(un);
        if (it != past.end() && it->second > 0) {
          --it->second;
          if (it->second == 1) new_correctors.insert(un);
        }
      }
      // u is a new corrector if u notin I and |N(u) cap unsolved| == 1
      if (!in) {
        past.insert({u, n_found});
        if (n_found == 1) new_correctors.insert(u);
      }
    }
    correctors = new_correctors;
    ++depth;
  } while (!correctors.empty());
  if (solved.size() != diag.n_vertices())
    throw ZXError("ZXDiagram does not have causal flow");
  return fl;
}

std::map<ZXVert, ZXVertSeqSet> Flow::gauss_solve_correctors(
    const ZXDiagram& diag, const boost::bimap<ZXVert, unsigned>& correctors,
    const boost::bimap<ZXVert, unsigned>& preserve, const ZXVertVec& to_solve,
    const boost::bimap<ZXVert, unsigned>& ys) {
  unsigned n_correctors = correctors.size();
  unsigned n_preserve = preserve.size();
  unsigned n_to_solve = to_solve.size();
  unsigned n_ys = ys.size();
  MatrixXb mat = MatrixXb::Zero(n_preserve + n_ys, n_correctors + n_to_solve);
  // Build adjacency matrix
  for (boost::bimap<ZXVert, unsigned>::const_iterator it = correctors.begin(),
                                                      end = correctors.end();
       it != end; ++it) {
    for (const ZXVert& n : diag.neighbours(it->left)) {
      auto in_past = preserve.left.find(n);
      if (in_past != preserve.left.end()) {
        mat(in_past->second, it->right) = true;
      } else {
        auto in_ys = ys.left.find(n);
        if (in_ys != ys.left.end()) {
          mat(n_preserve + in_ys->second, it->right) = true;
        }
      }
    }
  }
  for (boost::bimap<ZXVert, unsigned>::const_iterator it = ys.begin(),
                                                      end = ys.end();
       it != end; ++it) {
    auto found = correctors.left.find(it->left);
    if (found != correctors.left.end())
      mat(n_preserve + it->right, found->second) = true;
  }
  // Add rhs
  for (unsigned i = 0; i < n_to_solve; ++i) {
    ZXVert v = to_solve.at(i);
    switch (diag.get_zxtype(v)) {
      case ZXType::XY:
      case ZXType::PX: {
        mat(preserve.left.at(v), n_correctors + i) = true;
        break;
      }
      case ZXType::XZ: {
        mat(preserve.left.at(v), n_correctors + i) = true;
      }
      // fall through
      case ZXType::YZ:
      case ZXType::PZ: {
        for (const ZXVert& n : diag.neighbours(v)) {
          auto found = preserve.left.find(n);
          if (found != preserve.left.end())
            mat(found->second, n_correctors + i) = true;
          else {
            found = ys.left.find(n);
            if (found != ys.left.end())
              mat(n_preserve + found->second, n_correctors + i) = true;
          }
        }
        break;
      }
      case ZXType::PY: {
        mat(n_preserve + ys.left.at(v), n_correctors + i) = true;
        break;
      }
      default: {
        throw ZXError(
            "Internal error in flow identification: non-MBQC vertex found");
      }
    }
  }

  // Gaussian elimination
  std::vector<std::pair<unsigned, unsigned>> row_ops =
      gaussian_elimination_row_ops(
          mat.block(0, 0, n_preserve + n_ys, n_correctors));
  for (const std::pair<unsigned, unsigned>& op : row_ops) {
    for (unsigned j = 0; j < n_correctors + n_to_solve; ++j) {
      mat(op.second, j) ^= mat(op.first, j);
    }
  }

  // Back substitution
  // For each row i, pick a corrector j for which mat(i,j) == true, else
  // determine that row i has zero lhs
  std::map<unsigned, ZXVert> row_corrector;
  for (unsigned i = 0; i < n_preserve + n_ys; ++i) {
    for (unsigned j = 0; j < n_correctors; ++j) {
      if (mat(i, j)) {
        row_corrector.insert({i, correctors.right.at(j)});
        break;
      }
    }
  }
  // For each past i, scan down column of rhs and for each mat(j,CI+i) == true,
  // add corrector from row j or try next i if row j has zero lhs
  std::map<ZXVert, ZXVertSeqSet> solved_flow;
  for (unsigned i = 0; i < n_to_solve; ++i) {
    bool fail = false;
    ZXVertSeqSet c_i;
    for (unsigned j = 0; j < n_preserve + n_ys; ++j) {
      if (mat(j, n_correctors + i)) {
        auto found = row_corrector.find(j);
        if (found == row_corrector.end()) {
          fail = true;
          break;
        } else {
          c_i.insert(found->second);
        }
      }
    }
    if (!fail) {
      ZXVert v = to_solve.at(i);
      ZXType vt = diag.get_zxtype(v);
      if (vt == ZXType::XZ || vt == ZXType::YZ || vt == ZXType::PZ)
        c_i.insert(v);
      solved_flow.insert({v, c_i});
    }
  }
  return solved_flow;
}

Flow Flow::identify_pauli_flow(const ZXDiagram& diag) {
  // Check diagram has the expected form for pauli flow
  if (!diag.is_MBQC())
    throw ZXError("ZXDiagram must be in MBQC form to identify Pauli flow");

  ZXVertSeqSet solved;
  std::set<ZXVert> inputs;
  Flow fl{{}, {}};

  // Tag input measurements
  for (const ZXVert& i : diag.get_boundary(ZXType::Input)) {
    ZXVert ni = diag.neighbours(i).at(0);
    inputs.insert(ni);
    ZXType nt = diag.get_zxtype(ni);
    if (nt == ZXType::XZ || nt == ZXType::YZ || nt == ZXType::PY)
      throw ZXError(
          "Inputs measured in XZ, YZ, or Y cannot be corrected with Pauli "
          "flow");
  }

  // Indexing of correctors in binary matrix can be preserved between rounds as
  // we will only ever add new correctors
  boost::bimap<ZXVert, unsigned> correctors;
  unsigned corrector_i = 0;

  BGL_FORALL_VERTICES(v, *diag.graph, ZXGraph) {
    switch (diag.get_zxtype(v)) {
      case ZXType::Output: {
        ZXVert n = diag.neighbours(v).at(0);
        // Outputs are trivially solved
        solved.insert(v);
        if (diag.get_zxtype(n) != ZXType::Input) {
          solved.insert(n);
          fl.c_.insert({n, {}});
          fl.d_.insert({n, 0});
        }
        // n is either an Input or PX, in which case it will be added to the
        // correctors in the PX case
        break;
      }
      case ZXType::PX:
      case ZXType::PY: {
        // Can use non-input Xs and Ys to correct
        if (inputs.find(v) == inputs.end()) {
          correctors.insert({v, corrector_i});
          ++corrector_i;
        }
        break;
      }
      default:
        break;
    }
  }

  unsigned depth = 1;

  unsigned n_solved = 0;
  do {
    // Construct Gaussian elimination problem
    boost::bimap<ZXVert, unsigned> preserve;
    boost::bimap<ZXVert, unsigned> unsolved_ys;
    ZXVertVec to_solve;
    BGL_FORALL_VERTICES(v, *diag.graph, ZXGraph) {
      ZXType type = diag.get_zxtype(v);
      if (solved.get<TagKey>().find(v) == solved.get<TagKey>().end() &&
          type != ZXType::Input) {
        to_solve.push_back(v);
        if (type == ZXType::PY)
          unsolved_ys.insert({v, (unsigned)unsolved_ys.size()});
        else if (type != ZXType::PZ)
          preserve.insert({v, (unsigned)preserve.size()});
      }
    }

    std::map<ZXVert, ZXVertSeqSet> new_corrections = gauss_solve_correctors(
        diag, correctors, preserve, to_solve, unsolved_ys);

    n_solved = new_corrections.size();

    for (const std::pair<const ZXVert, ZXVertSeqSet>& nc : new_corrections) {
      fl.c_.insert(nc);
      fl.d_.insert({nc.first, depth});
      solved.insert(nc.first);
      if (inputs.find(nc.first) == inputs.end())
        correctors.insert({nc.first, (unsigned)correctors.size()});
    }

    ++depth;
  } while (n_solved != 0);

  if (solved.size() + inputs.size() != diag.n_vertices())
    throw ZXError("ZXDiagram does not have pauli flow");

  return fl;
}

std::set<ZXVertSeqSet> Flow::identify_focussed_sets(const ZXDiagram& diag) {
  // Check diagram has the expected form for pauli flow
  if (!diag.is_MBQC())
    throw ZXError("ZXDiagram must be in MBQC form to identify gflow");

  std::set<ZXVert> inputs;
  std::set<ZXVert> outputs;

  // Tag input measurements
  for (const ZXVert& i : diag.get_boundary(ZXType::Input)) {
    ZXVert ni = diag.neighbours(i).at(0);
    inputs.insert(ni);
    ZXType nt = diag.get_zxtype(ni);
    if (nt == ZXType::XZ || nt == ZXType::YZ || nt == ZXType::PY)
      throw ZXError(
          "Inputs measured in XZ, YZ, or Y cannot be corrected with Pauli "
          "flow");
  }
  for (const ZXVert& o : diag.get_boundary(ZXType::Output)) {
    ZXVert no = diag.neighbours(o).at(0);
    outputs.insert(no);
  }

  // Build Gaussian elimination problem
  boost::bimap<ZXVert, unsigned> correctors;
  boost::bimap<ZXVert, unsigned> preserve;
  boost::bimap<ZXVert, unsigned> ys;
  unsigned n_correctors = 0;
  unsigned n_preserve = 0;
  unsigned n_ys = 0;

  BGL_FORALL_VERTICES(v, *diag.graph, ZXGraph) {
    switch (diag.get_zxtype(v)) {
      case ZXType::XY: {
        preserve.insert({v, n_preserve});
        ++n_preserve;
        // Can use non-input Xs and Ys to correct
        if (inputs.find(v) == inputs.end()) {
          correctors.insert({v, n_correctors});
          ++n_correctors;
        }
        break;
      }
      case ZXType::PX: {
        // Nonmeasured vertices also covered by PX
        // Only need to preserve measured vertices
        if (outputs.find(v) == outputs.end()) {
          preserve.insert({v, n_preserve});
          ++n_preserve;
        }
        // Can use non-input Xs and Ys to correct
        if (inputs.find(v) == inputs.end()) {
          correctors.insert({v, n_correctors});
          ++n_correctors;
        }
        break;
      }
      case ZXType::PY: {
        ys.insert({v, n_ys});
        ++n_ys;
        // Can use non-input Xs and Ys to correct
        if (inputs.find(v) == inputs.end()) {
          correctors.insert({v, n_correctors});
          ++n_correctors;
        }
        break;
      }
      default:
        break;
    }
  }

  MatrixXb mat = MatrixXb::Zero(n_preserve + n_ys, n_correctors);

  // Build adjacency matrix
  for (boost::bimap<ZXVert, unsigned>::const_iterator it = correctors.begin(),
                                                      end = correctors.end();
       it != end; ++it) {
    for (const ZXVert& n : diag.neighbours(it->left)) {
      auto in_preserve = preserve.left.find(n);
      if (in_preserve != preserve.left.end()) {
        mat(in_preserve->second, it->right) = true;
      } else {
        auto in_ys = ys.left.find(n);
        if (in_ys != ys.left.end()) {
          mat(n_preserve + in_ys->second, it->right) = true;
        }
      }
    }
  }
  for (boost::bimap<ZXVert, unsigned>::const_iterator it = ys.begin(),
                                                      end = ys.end();
       it != end; ++it) {
    auto found = correctors.left.find(it->left);
    if (found != correctors.left.end())
      mat(n_preserve + it->right, found->second) = true;
  }

  // Gaussian elimination
  std::vector<std::pair<unsigned, unsigned>> row_ops =
      gaussian_elimination_row_ops(mat);
  for (const std::pair<unsigned, unsigned>& op : row_ops) {
    for (unsigned j = 0; j < n_correctors; ++j) {
      mat(op.second, j) ^= mat(op.first, j);
    }
  }

  // Back substitution
  // For each column j, it either a leading column (the first column for which
  // mat(i,j) == true for a given i, so set row_corrector[i] = j; by Gaussian
  // Elimination this is the only entry in the column) or it describes the
  // focussed set generator {j} + {row_corrector[i] | mat(i,j) == true}
  std::set<ZXVertSeqSet> focussed;
  std::map<unsigned, ZXVert> row_corrector;
  for (boost::bimap<ZXVert, unsigned>::const_iterator it = correctors.begin(),
                                                      end = correctors.end();
       it != end; ++it) {
    ZXVertSeqSet fset{it->left};
    bool new_row_corrector = false;
    for (unsigned i = 0; i < n_preserve + n_ys; ++i) {
      if (mat(i, it->right)) {
        auto inserted = row_corrector.insert({i, it->left});
        if (inserted.second) {
          // New row_corrector, so move to next column
          new_row_corrector = true;
          break;
        } else {
          // Non-correcting column
          fset.insert(inserted.first->second);
        }
      }
    }
    if (!new_row_corrector) focussed.insert({fset});
  }

  return focussed;
}

}  // namespace zx

}  // namespace tket
