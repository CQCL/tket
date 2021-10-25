#ifndef _TKET_TokenSwapping_TSAUtils_DebugFunctions_H_
#define _TKET_TokenSwapping_TSAUtils_DebugFunctions_H_

#include <string>

#include "VertexMappingFunctions.hpp"

namespace tket {
namespace tsa_internal {

/** Get a string representation.
 *  @param vertex_mapping A mapping, usually representing a desired
 * source->target mapping for a Token Swapping problem.
 *  @return A string representation.
 */
std::string str(const VertexMapping& vertex_mapping);

/** Get a string representation.
 *  @param swaps An ordered list of swaps, usually the solution to a Token
 * Swapping problem.
 *  @return A string representation.
 */
std::string str(const SwapList& swaps);

/** Get a string representation.
 *  @param swaps An ordered list of swaps, usually the solution to a Token
 * Swapping problem.
 *  @return A string representation.
 */
std::string str(const std::vector<Swap>& swaps);

}  // namespace tsa_internal
}  // namespace tket
#endif
