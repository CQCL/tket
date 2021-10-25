
#ifndef _TKET_TokenSwapping_TableLookup_PartialMappingLookup_H_
#define _TKET_TokenSwapping_TableLookup_PartialMappingLookup_H_

#include <cstdint>
#include <set>
#include <vector>

#include "ExactMappingLookup.hpp"

namespace tket {
namespace tsa_internal {

/** This is the same as ExactMappingLookup, except that we allow vertices not to
 * have tokens. It works simply by going through possible permutations of empty
 * vertices and doing an exact permutation lookup (limiting the number of
 * permutations to avoid excessive slowdown).
 */
class PartialMappingLookup {
 public:
  /** Parameters controlling the partial mapping lookup. Sensible defaults,
   * found by experimentation. */
  struct Parameters {
    /** To speed up, don't try all permutations if there are many empty
     * vertices; limit them to this number. */
    unsigned max_number_of_empty_vertex_permutations;

    Parameters();
  };

  /** If desired, change some internal parameters.
   * @return Internal parameters object, to be changed if desired.
   */
  Parameters& get_parameters();

  /** The result is stored internally. The same format as ExactMappingLookup.
   * @param desired_mapping A (source vertex) -> (target vertex) permutation.
   * @param edges Edges which exist between the vertices (equivalently, the
   * swaps which we are permitted to use). Edges with vertices not appearing in
   * desired_mapping will simply be ignored.
   * @param vertices_with_tokens_at_start Every vertex mentioned within
   * desired_mapping which has a token, just BEFORE the swaps are performed to
   * enact the desired_mapping, must be mentioned here. Other vertices not
   * mentioned in the mapping are allowed; they will simply be ignored.
   * @param max_number_of_swaps Stop looking if every sequence of swaps in the
   * table which enacts the desired mapping exceeds this length (or doesn't
   * exist at all).
   */
  const ExactMappingLookup::Result& operator()(
      const VertexMapping& desired_mapping, const std::vector<Swap>& edges,
      const std::set<size_t>& vertices_with_tokens_at_start,
      unsigned max_number_of_swaps = 16);

 private:
  Parameters m_parameters;
  ExactMappingLookup m_exact_mapping_lookup;
  std::vector<size_t> m_empty_source_vertices;
  std::vector<size_t> m_empty_target_vertices;
  VertexMapping m_altered_mapping;
};

}  // namespace tsa_internal
}  // namespace tket
#endif
