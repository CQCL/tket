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

#include "ZX/ZXDiagram.hpp"

namespace tket {

namespace zx {

/**
 * Class for compositional rewrites.
 * The broad structure is similar to the Transform class for Circuits.
 * Each rewrite is encapsulated by a single method `apply` which searches for
 * and performs all instances of a rewrite on a given diagram. The `apply`
 * method has an optional `Subdiagram` argument to operate locally on some
 * region of the diagram. This will also update the `Subdiagram` to use the new
 * vertices and edges. The `apply` method returns whether or not a rewrite was
 * performed. The main way to compose rewrites is in sequence or loops. Each
 * rewrite is expected to preserve global semantics of the diagram (or
 * restricted subdiagram) including global scalar.
 *
 * Previous designs used a separate method to identify independent matches from
 * performing the rewrite on a single match. This added extra complexity to the
 * structure, more reliance on users to be sensible, and reduced efficiency of
 * the rewrites. We decided that interactive rewriting was something we were
 * happy to go without as tket is not intended to be a proof assistant.
 */
class Rewrite {
 public:
  typedef std::function<bool(ZXDiagram&)> RewriteFun;
  typedef std::function<unsigned(const ZXDiagram&)> Metric;

  /**
   * The actual rewrite to be applied.
   * Performs the rewrite in place (optionally restricted to some subdiagram)
   * and returns true iff some change is made.
   */
  const RewriteFun apply;

  ///////////////
  // Combinators//
  ///////////////

  static Rewrite sequence(const std::vector<Rewrite>& rvec);
  static Rewrite repeat(const Rewrite& rw);
  static Rewrite repeat_with_metric(const Rewrite& rw, const Metric& eval);
  static Rewrite repeat_while(const Rewrite& cond, const Rewrite& body);

  //////////////////
  // Decompositions//
  //////////////////

  /**
   * Replaces every ZXBox by its internal diagram recursively until no ZXBoxes
   * remain.
   */
  static Rewrite decompose_boxes();

  /**
   * Replaces every Hadamard wire by an explicit Hbox node.
   */
  static Rewrite basic_wires();

  /**
   * Expands every generator into ZSpiders, XSpiders, and a combination of Basic
   * and Hadamard edges.
   */
  static Rewrite rebase_to_zx();

  /**
   * Expands every generator into MBQC vertices.
   */
  static Rewrite rebase_to_mbqc();

  //////////
  // Axioms//
  //////////

  /**
   * Converts all red spiders (XSpider) to green (ZSpider) with Hadamards around
   * it. The Hadamards are applied by flipping the wire type of incident edges
   * between Basic and H.
   */
  static Rewrite red_to_green();

  /**
   * Merges two adjacent ZX spiders (XSpider, ZSpider) of the same colour
   * connected by a Basic wire into a single spider. Also merges two adjacent
   * spiders of different colour connected by a H edge.
   */
  static Rewrite spider_fusion();

  /**
   * Removes both H and Basic self loop edges around ZX spiders.
   * Basic edges can simply be removed.
   * Removing H loops introduces an extra pi phase on the spider.
   */
  static Rewrite self_loop_removal();

  /**
   * Remove parallel edges between ZX spiders (Hopf rule).
   * Matches either pairs of H edges between spiders of the same colour or Basic
   * edges between spiders of different colour. This applies to Quantum edges
   * between a pair of Classical spiders.
   */
  static Rewrite parallel_h_removal();

  /////////////////
  // GraphLikeForm//
  /////////////////

  /**
   * Guarantees that each boundary vertex is adjacent to a unique ZSpider.
   * This adds identity chains when two boundaries are either directly connected
   * or are adjacent to the same spider.
   */
  static Rewrite separate_boundaries();

  /**
   * Guarantees that the edge on each boundary vertex is Basic.
   * If a boundary has a Hadamard, then we add a ZSpider identity as in I/O
   * extensions in MBQC.
   */
  static Rewrite io_extension();

  ///////////////////////////
  // GraphLikeSimplification//
  ///////////////////////////

  /**
   * Removes interior proper Cliffords (spiders where the phase is an odd
   * multiple of pi/2). Performs local complementation about the vertex and
   * removes it.
   */
  static Rewrite remove_interior_cliffords();

  /**
   * Removes adjacent interior Paulis (spiders where the phase is an integer
   * multiple of pi). Pivots about the edge connecting the vertices and removes
   * them.
   */
  static Rewrite remove_interior_paulis();

  /**
   * Identifies adjacent Pauli spiders where one is adjacent to a boundary.
   * This rule applies I/O extensions to push the match into the interior from
   * which it can be handled by `remove_interior_paulis`.
   */
  static Rewrite extend_at_boundary_paulis();

 private:
  Rewrite(const RewriteFun& fun);

  static bool decompose_boxes_fun(ZXDiagram& diag);
  static bool basic_wires_fun(ZXDiagram& diag);
  static bool rebase_to_zx_fun(ZXDiagram& diag);
  static bool rebase_to_mbqc_fun(ZXDiagram& diag);
  static bool red_to_green_fun(ZXDiagram& diag);
  static bool spider_fusion_fun(ZXDiagram& diag);
  static bool self_loop_removal_fun(ZXDiagram& diag);
  static bool parallel_h_removal_fun(ZXDiagram& diag);
  static bool separate_boundaries_fun(ZXDiagram& diag);
  static bool io_extension_fun(ZXDiagram& diag);
  static bool remove_interior_cliffords_fun(ZXDiagram& diag);
  static bool remove_interior_paulis_fun(ZXDiagram& diag);
  static bool extend_at_boundary_paulis_fun(ZXDiagram& diag);
};

}  // namespace zx

}  // namespace tket
