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

#include "Converters.hpp"
#include "ZX/Flow.hpp"

namespace tket {

using namespace zx;

void clean_frontier(
    ZXDiagram& diag, ZXVertVec& frontier, Circuit& circ,
    std::map<ZXVert, unsigned>& qubit_map) {
  std::set<ZXVert> frontier_lookup;
  std::list<std::pair<ZXVert, ZXVert>> czs;
  for (const ZXVert& f : frontier) {
    frontier_lookup.insert(f);
    for (const Wire& w : diag.adj_wires(f)) {
      ZXVert n = diag.other_end(w, f);
      if (diag.get_zxtype(n) == ZXType::Output) {
        unsigned q = qubit_map.at(n);
        if (diag.get_wire_type(w) == ZXWireType::H) {
          diag.set_wire_type(w, ZXWireType::Basic);
          circ.add_op<unsigned>(OpType::H, {q});
        }
        switch (diag.get_zxtype(f)) {
          case ZXType::Input: {
            frontier.erase(std::find(frontier.begin(), frontier.end(), f));
            break;
          }
          case ZXType::XY: {
            const PhasedGen& f_gen = diag.get_vertex_ZXGen<PhasedGen>(f);
            Expr ph = f_gen.get_param();
            if (!equiv_0(ph)) {
              circ.add_op<unsigned>(OpType::U1, ph, {q});
              diag.set_vertex_ZXGen_ptr(
                  f, ZXGen::create_gen(ZXType::PX, *f_gen.get_qtype()));
            }
            break;
          }
          case ZXType::PX: {
            const CliffordGen& f_gen = diag.get_vertex_ZXGen<CliffordGen>(f);
            if (f_gen.get_param()) {
              circ.add_op<unsigned>(OpType::Z, {q});
              diag.set_vertex_ZXGen_ptr(
                  f, ZXGen::create_gen(ZXType::PX, *f_gen.get_qtype()));
            }
            break;
          }
          case ZXType::PY: {
            const CliffordGen& f_gen = diag.get_vertex_ZXGen<CliffordGen>(f);
            circ.add_op<unsigned>(
                f_gen.get_param() ? OpType::S : OpType::Sdg, {q});
            diag.set_vertex_ZXGen_ptr(
                f, ZXGen::create_gen(ZXType::PX, *f_gen.get_qtype()));
            break;
          }
          default:
            throw ZXError(
                "Error during extraction from ZX diagram: unexpected ZXType in "
                "frontier");
        }
      } else if (frontier_lookup.find(n) != frontier_lookup.end()) {
        czs.push_back({f, n});
        diag.remove_wire(w);
      }
    }
  }
  for (const std::pair<ZXVert, ZXVert>& pair : czs) {
    circ.add_op<unsigned>(
        OpType::CZ, {qubit_map.at(pair.first), qubit_map.at(pair.second)});
  }
}

ZXVertSeqSet neighbours_of_frontier(
    const ZXDiagram& diag, ZXVertVec& frontier) {
  ZXVertSeqSet n_set;
  for (const ZXVert& f : frontier) {
    for (const Wire& w : diag.adj_wires(f)) {
      ZXVert n = diag.other_end(w, f);
      ZXType n_type = diag.get_zxtype(n);
      if (n_type == ZXType::Output || n_type == ZXType::Input) continue;
      n_set.insert(n);
    }
  }
  return n_set;
}

static void bipartite_complementation(
    ZXDiagram& diag, const ZXVertSeqSet& sa, const ZXVertSeqSet& sb) {
  for (const ZXVert& a : sa.get<TagSeq>()) {
    for (const ZXVert& b : sb.get<TagSeq>()) {
      std::optional<Wire> wire = diag.wire_between(a, b);
      if (wire)
        diag.remove_wire(*wire);
      else
        diag.add_wire(a, b, ZXWireType::H);
    }
  }
}

bool remove_all_gadgets(ZXDiagram& diag, const ZXVertVec& frontier) {
  bool removed_gadget = false;
  for (const ZXVert& f : frontier) {
    ZXVertVec f_ns = diag.neighbours(f);
    ZXVert o;
    for (const ZXVert& n : f_ns) {
      // Each frontier vertex is connected to a unique output, find it
      if (diag.get_zxtype(n) == ZXType::Output) {
        o = n;
        break;
      }
    }
    for (const ZXVert& n : f_ns) {
      if (diag.get_zxtype(n) == ZXType::YZ) {
        // Pivot
        // Identify three subsets of neighbours
        ZXVertSeqSet excl_f{f_ns.begin(), f_ns.end()};
        excl_f.erase(n);
        excl_f.erase(o);
        ZXVertSeqSet excl_n, joint;
        auto& lookup_f = excl_f.get<TagKey>();
        for (const ZXVert& nn : diag.neighbours(n)) {
          if (lookup_f.find(nn) != lookup_f.end())
            joint.insert(nn);
          else
            excl_n.insert(nn);
        }
        excl_n.erase(f);
        excl_f.erase(joint.begin(), joint.end());
        // The is_MBQC check in zx_to_circuit guarantees QuantumType::Quantum
        bipartite_complementation(diag, joint, excl_n);
        bipartite_complementation(diag, joint, excl_f);
        bipartite_complementation(diag, excl_n, excl_f);
        Wire ow = *diag.wire_between(f, o);
        diag.set_wire_type(
            ow, (diag.get_wire_type(ow) == ZXWireType::Basic)
                    ? ZXWireType::H
                    : ZXWireType::Basic);
        const PhasedGen& yz_gen = diag.get_vertex_ZXGen<PhasedGen>(n);
        diag.set_vertex_ZXGen_ptr(
            n, ZXGen::create_gen(
                   ZXType::XY, -yz_gen.get_param(), QuantumType::Quantum));
        for (const ZXVert& nn : joint.get<TagSeq>()) {
          ZXGen_ptr new_gen;
          switch (diag.get_zxtype(nn)) {
            case ZXType::XY: {
              const PhasedGen& ph_gen = diag.get_vertex_ZXGen<PhasedGen>(nn);
              new_gen = ZXGen::create_gen(
                  ZXType::XY, ph_gen.get_param() + 1., QuantumType::Quantum);
              break;
            }
            case ZXType::PX:
            case ZXType::PY: {
              const CliffordGen& cl_gen =
                  diag.get_vertex_ZXGen<CliffordGen>(nn);
              new_gen = ZXGen::create_gen(
                  cl_gen.get_type(), !cl_gen.get_param(), QuantumType::Quantum);
              break;
            }
            case ZXType::XZ:
            case ZXType::YZ: {
              const PhasedGen& ph_gen = diag.get_vertex_ZXGen<PhasedGen>(nn);
              new_gen = ZXGen::create_gen(
                  ph_gen.get_type(), -ph_gen.get_param(), QuantumType::Quantum);
              break;
            }
            case ZXType::PZ: {
              new_gen = diag.get_vertex_ZXGen_ptr(nn);
              break;
            }
            default:
              throw ZXError(
                  "Error during extraction from ZX diagram: unexpected ZXType "
                  "during local complementation");
          }
          diag.set_vertex_ZXGen_ptr(nn, new_gen);
        }
        removed_gadget = true;
      }
    }
  }
  return removed_gadget;
}

Circuit zx_to_circuit(const ZXDiagram& d) {
  ZXDiagram diag = d;
  if (!diag.is_MBQC())
    throw ZXError("Can only extract a circuit from a ZX diagram in MBQC form");
  ZXVertVec ins = diag.get_boundary(ZXType::Input);
  ZXVertVec outs = diag.get_boundary(ZXType::Output);
  if (ins.size() != outs.size())
    throw ZXError("Can only extract a circuit from a unitary ZX diagram");

  Circuit circ(ins.size());

  ZXVertVec frontier;
  std::map<ZXVert, unsigned> qubit_map;
  for (const ZXVert& o : outs) {
    ZXVert f_i = diag.neighbours(o).at(0);
    frontier.push_back(f_i);
    qubit_map.insert({o, qubit_map.size()});
    qubit_map.insert({f_i, qubit_map.size()});
  }

  clean_frontier(diag, frontier, circ, qubit_map);
  while (!frontier.empty()) {
    if (remove_all_gadgets(diag, frontier)) continue;

    ZXVertSeqSet neighbours = neighbours_of_frontier(diag, frontier);
    boost::bimap<ZXVert, unsigned> correctors, preserve, ys;
    ZXVertVec to_solve;
    for (const ZXVert& f : frontier)
      correctors.insert({f, (unsigned)correctors.size()});
    for (const ZXVert& n : neighbours.get<TagSeq>()) {
      ZXType n_type = diag.get_zxtype(n);
      if (n_type == ZXType::XY || n_type == ZXType::PX ||
          n_type == ZXType::PY) {
        preserve.insert({n, (unsigned)preserve.size()});
        to_solve.push_back(n);
      }
    }
    std::map<ZXVert, ZXVertSeqSet> candidates =
        Flow::gauss_solve_correctors(diag, correctors, preserve, to_solve, ys);

    if (candidates.empty())
      throw ZXError(
          "Error during extraction from ZX diagram: diagram does not have "
          "gflow");

    unsigned min = UINT_MAX;
    ZXVert best;
    for (const std::pair<const ZXVert, ZXVertSeqSet>& p : candidates) {
      if (p.second.size() < min) {
        min = p.second.size();
        best = p.first;
      }
    }
    ZXVertSeqSet g_best = candidates.at(best);

    ZXVert f_to_isolate = g_best.get<TagSeq>().front();
    unsigned f_q = qubit_map.at(f_to_isolate);
    for (const ZXVert& f : g_best) {
      if (f != f_to_isolate) {
        circ.add_op<unsigned>(OpType::CX, {f_q, qubit_map.at(f)});
      }
    }
    ZXVert out;
    Wire w_out;
    for (const Wire& w : diag.adj_wires(f_to_isolate)) {
      ZXVert n = diag.other_end(w, f_to_isolate);
      if (diag.get_zxtype(n) == ZXType::Output) {
        out = n;
        w_out = w;
        break;
      }
    }
    diag.add_wire(
        best, out,
        (diag.get_wire_type(w_out) == ZXWireType::Basic) ? ZXWireType::H
                                                         : ZXWireType::Basic);
    diag.remove_vertex(f_to_isolate);
    qubit_map.erase(qubit_map.find(f_to_isolate));
    qubit_map.insert({best, f_q});
    for (ZXVertVec::iterator it = frontier.begin(); it != frontier.end();
         ++it) {
      if (*it == f_to_isolate) {
        *it = best;
        break;
      }
    }

    clean_frontier(diag, frontier, circ, qubit_map);
    // Remove inputs from frontier
    neighbours_of_frontier(diag, frontier);
  }

  qubit_map_t qm;
  for (unsigned i = 0; i < ins.size(); ++i) {
    ZXVert in = ins.at(i);
    ZXVert out = diag.neighbours(in).at(0);
    if (diag.get_zxtype(out) != ZXType::Output)
      throw ZXError(
          "Error during extraction from ZX diagram: input not adjacent to "
          "output after extracting");
    qm.insert({Qubit(qubit_map.at(out)), Qubit(i)});
  }
  circ.permute_boundary_output(qm);

  // Reverse gates in circuit (all gates added are self-transpose)
  return circ.transpose();
}

}  // namespace tket