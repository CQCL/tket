#ifndef _TKET_TokenSwapping_NeighboursFromArchitecture_H_
#define _TKET_TokenSwapping_NeighboursFromArchitecture_H_

#include "ArchitectureMapping.hpp"
#include "NeighboursInterface.hpp"

namespace tket {
namespace tsa_internal {

/** Stores and returns upon request the adjacent vertices to a given vertex
 *  on a graph, using an underlying Architecture object.
 */
class NeighboursFromArchitecture : public NeighboursInterface {
 public:
  /** The objects must remain valid AND unchanged
   *  for the lifetime of this object.
   *  @param arch_mapping An object which contains a reference to an
   *    Architecture object internally, and handles Node -> vertex size_t
   * conversions.
   */
  explicit NeighboursFromArchitecture(const ArchitectureMapping& arch_mapping);

  /** For extra convenience, the list of neighbours is always sorted
   *  in increasing order (so you can do binary search, etc.)
   *  @param vertex A vertex.
   *  @return A sorted list of all adjacent vertices, stored internally.
   */
  virtual const std::vector<size_t>& operator()(size_t vertex) override;

 private:
  const ArchitectureMapping& m_arch_mapping;

  /** The key is the vertex, the value is the list of neighbours. */
  std::map<size_t, std::vector<size_t>> m_cached_neighbours;
};

}  // namespace tsa_internal
}  // namespace tket
#endif
