#ifndef _TKET_TESTS_TokenSwapping_TestUtils_PartialTsaTesting_H_
#define _TKET_TESTS_TokenSwapping_TestUtils_PartialTsaTesting_H_

#include "Architecture/Architectures.hpp"
#include "TokenSwapping/PartialTsaInterface.hpp"
#include "TokenSwapping/RNG.hpp"

namespace tket {
namespace tsa_internal {
namespace tests {

enum class RequiredTsaProgress { NONE, FULL, NONZERO };
enum class TokenOption {
  ALLOW_EMPTY_TOKEN_SWAP,
  DO_NOT_ALLOW_EMPTY_TOKEN_SWAP
};

/// Returns a summary string of the results, as well as doing the checks.
std::string run_tests(
    const Architecture& arch, const std::vector<VertexMapping>& problems,
    PathFinderInterface& path_finder, PartialTsaInterface& partial_tsa,
    RequiredTsaProgress progress,
    TokenOption token_option = TokenOption::DO_NOT_ALLOW_EMPTY_TOKEN_SWAP);

/// If no path finder is specified, will use the RiverFlowPathFinder
/// (which needs an RNG).
std::string run_tests(
    const Architecture& arch, const std::vector<VertexMapping>& problems,
    RNG& rng, PartialTsaInterface& partial_tsa, RequiredTsaProgress progress,
    TokenOption token_option = TokenOption::DO_NOT_ALLOW_EMPTY_TOKEN_SWAP);

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket
#endif
