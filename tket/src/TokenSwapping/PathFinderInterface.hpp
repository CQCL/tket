#ifndef _TKET_TokenSwapping_PathFinderInterface_H_
#define _TKET_TokenSwapping_PathFinderInterface_H_

#include <string>
#include <vector>

namespace tket {
namespace tsa_internal {

/** What is SOME shortest path between vertices?
 *  This might involve an arbitrary choice,
 *  because some paths will not be unique if the graph is not a tree.
 *  For algorithms, we might need to choose in a vaguely consistent way,
 *  and use a random number generator.
 */
class PathFinderInterface {
 public:
  PathFinderInterface();

  /** By default, simply throws (not implemented).
   *  Returns a shortest path from v1 to v2, including v1 at the start
   *  and v2 at the end. This should usually return the same result for
   *  (v1, v2) each time it is called, but may change slightly over time.
   *  Although the path is stored internally, there's no guarantee
   *  that the reference will remain valid once another function call occurs.
   *  There's no guarantee that the path for (v1, v2) will be the reverse of
   *  the path for (v2, v1).
   *  Could take time O(length of path), if it is built up anew each time.
   *  @param vertex1 First vertex v1.
   *  @param vertex2 Second vertex v2.
   *  @return A list of vertices, starting with v1 and ending with v2,
   *    giving a shortest path from v1 to v2 (not unique, maybe not constant
   *    over time, and maybe not a valid reference after any other call).
   */
  virtual const std::vector<size_t>& operator()(size_t vertex1, size_t vertex2);

  virtual ~PathFinderInterface();

  /** Some path finders use randomness; if so, override this to reset
   *  the source of randomness to some default seed
   *  to ensure reproducibility. By default, does nothing.
   */
  virtual void reset();

  /** If some other algorithm has made use of an edge v1-v2,
   * without going through this path finder object,
   * call this function to inform this object.
   * (E.g., some classes remember which previous operator() calls were made,
   * and use them to decide future paths when there is a nonunique choice).
   * By default, does nothing.
   * @param vertex1 First vertex v1.
   * @param vertex2 Second vertex v2.
   */
  virtual void register_edge(size_t vertex1, size_t vertex2);

  /** For convenience, if "register_edge" does nothing, return false so that the
   * caller knows and doesn't waste time repeatedly calling "register_edge".
   * @return True if the function "register_edge" has been overridden to do
   * something, false if the function does nothing
   */
  virtual bool edge_registration_has_effect() const;

  /** For debugging purposes, every object has a name.
   *  @return The name of the object.
   */
  const std::string& name() const;

 protected:
  std::string m_name;
};

}  // namespace tsa_internal
}  // namespace tket
#endif
