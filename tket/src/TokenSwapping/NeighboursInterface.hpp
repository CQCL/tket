#ifndef _TKET_TokenSwapping_NeighboursInterface_H_
#define _TKET_TokenSwapping_NeighboursInterface_H_

#include <cstdint>
#include <vector>

namespace tket {
namespace tsa_internal {

/** What are the adjacent vertices to a given vertex on a graph?
 *  For larger, sparse graphs, it might
 *  calculate and store neighbours only when required.
 */
class NeighboursInterface {
 public:
  /** Returns the neighbours of the given vertex.
   *  The vector of neighbours is required to be stored internally.
   *  However, no guarantee that the reference will remain valid
   *  once another function call occurs.
   *  By default, throws (not implemented).
   *  (It's useful to be able to create a "null" object like this,
   *  because some algorithms don't actually need a neighbours object,
   *  but others do).
   *  @param vertex A vertex.
   *  @return A sorted list of all adjacent vertices, stored internally.
   */
  virtual const std::vector<size_t>& operator()(size_t vertex);

  virtual ~NeighboursInterface();
};

}  // namespace tsa_internal
}  // namespace tket
#endif
