#include <algorithm>
#include <array>
#include <iostream>
#include <iterator>
#include <memory>
#include <numeric>
#include <set>
#include <string>
#include <functional>

#include "Circuit/CircUtils.hpp"
#include "Circuit/Circuit.hpp"
#include "Circuit/DAGDefs.hpp"
#include "Circuit/ThreeQubitConversion.hpp"
#include "OpType/EdgeType.hpp"
#include "OpType/OpType.hpp"
#include "Utils/Assert.hpp"
#include "Utils/GraphHeaders.hpp"

namespace tket {

// Helper class defining a pure-quantum subcircuit of up to n qubits.
class QInteraction {
 public:
  // Construct an empty subcircuit with a single edge that is both input and
  // output, and no vertices.
  QInteraction(const Circuit &circ, const Edge &e);

  // Combine with another subcircuit disjoint from this. Disjointness, and the
  // fact that the combined subcircuit has at most three wires, is assumed.
  void combine(const QInteraction &other);

  EdgeVec out_edges() const;

  VertexSet vertices() const;

  unsigned n_wires() const;

  unsigned n_vertices() const;

  Subcircuit subcircuit() const;

  // Append a vertex following the subcircuit. It is assumed that every input
  // edge of the vertex matches exactly one output edge of the existing
  // subcircuit.
  void append(const Vertex &v);

 private:
  const Circuit &circ_;
  EdgeVec in_edges_;
  EdgeVec out_edges_;
  unsigned n_wires_;    // number of in/out edge pairs
  VertexSet vertices_;  // all internal vertices
};

typedef std::unique_ptr<QInteraction> iptr;

// Helper class representing a system of disjoint interactions, each with at
// most three qubits. The interactions are represented by integer labels.
class QISystem {
 public:
  // Construct an empty system.
  explicit QISystem(Circuit &circ, std::function<Circuit(Circuit)> replacement_func);

  // Add a new interaction to the system consisting of a single edge, and
  // return its index.
  int create_new_interaction_from_edge(const Edge &e);

  // Return the set of (indices of) interactions in the system that have v as a
  // direct successor.
  std::vector<int> interactions_feeding_vertex(const Vertex &v) const;

  // The total width (number of wires) of a subset of the interactions.
  unsigned total_n_wires(const std::vector<int> &S) const;

  // From a set of indices, choose the one indexing the largest interaction,
  // in terms of vertex count.
  int largest_interaction(const std::vector<int> &S) const;

  // Combine a set of existing interactions into one and append the vertex v.
  // It is assumed that the interactions are combinable and the vertex
  // appendable.
  void combine_and_append(const std::vector<int> &S, const Vertex &v);

  // Close an interaction, squashing it if possible, and erasing it from the
  // set. Return true iff any substitution was made, and the (possibly new)
  // vector of outgoing edges from the region of the interaction.
  std::pair<bool, EdgeVec> close_interaction(int i);

  // Close an interaction and spawn new ones on its outgoing edges. Return true
  // iff any substitution is made.
  bool close_interaction_and_spawn(int i);

  // Close all interactions that have v as a direct successor, and start new
  // ones following them. Return true iff any substitution is made.
  bool close_interactions_feeding_vertex(const Vertex &v);

  // Close all interactions. Return true iff any substitution is made.
  bool close_all_interactions();

  // Delete all vertices marked for deletion.
  void destroy_bin();

 private:
  Circuit &circ_;
  VertexList bin_;
  std::map<int, iptr> interactions_;
  int idx_;
  std::function<Circuit(Circuit)> replacement_func_;
};

}  // namespace tket