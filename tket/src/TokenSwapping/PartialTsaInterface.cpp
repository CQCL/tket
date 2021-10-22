#include "PartialTsaInterface.hpp"

#include "Utils/Exceptions.hpp"

namespace tket {
namespace tsa_internal {

const std::string& PartialTsaInterface::name() const { return m_name; }

}  // namespace tsa_internal
}  // namespace tket
