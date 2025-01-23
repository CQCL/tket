// Copyright Quantinuum
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

#include "tkwsm/Common/DyadicFraction.hpp"

#include <cmath>
#include <sstream>
#include <tkassert/Assert.hpp>

#include "tkwsm/Common/BitFunctions.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

DyadicFraction::DyadicFraction(UInt x) : m_value(x), m_exponent(0) {}

// We have 0 < x <= y, but x*y may overflow,
// and x,y have already had all their trailing zeros removed.
// If x*y DOES overflow, right shift x,y as little as possible,
// until x*y can safely be calculated,
// and return the total right shifts.
static unsigned remove_bits_lossy(
    DyadicFraction::UInt& x, unsigned x_bit_length, DyadicFraction::UInt& y,
    unsigned y_bit_length) {
  unsigned total_bitlength = x_bit_length + y_bit_length;

  // OK, some slight inaccuracy here;
  //    x*y would overflow  <==>  BitLength(xy) > 64,
  // but of course we're using BitLength(x) + BitLength(y)
  // as an estimate (upper bound) for BitLength(xy).
  // Might be out by 1, so we might remove 1 extra bit;
  // but don't worry about it, insignificant.
  if (total_bitlength <= 64) {
    return 0;
  }
  unsigned remaining_shift = total_bitlength - 64;
  TKET_ASSERT(remaining_shift >= 1);
  // We're definitely going to shift by exactly this much,
  // so might as well store it now.
  const auto shift_to_return = remaining_shift;

  // First, y is larger, so put all the reduction onto y, until
  // x,y are equal in length.
  const unsigned ly_minus_lx = y_bit_length - x_bit_length;
  if (ly_minus_lx >= remaining_shift) {
    // We're lucky! We need only shift y.
    y >>= remaining_shift;
    return shift_to_return;
  }
  // Shift y to get an equal bitlength with x.
  y >>= ly_minus_lx;
  remaining_shift -= ly_minus_lx;
  // Now, split the shifts equally between x,y,
  // since they now have equal bitlength.
  const unsigned common_shift = remaining_shift / 2;
  x >>= common_shift;
  y >>= common_shift;
  remaining_shift -= common_shift * 2;
  TKET_ASSERT(remaining_shift <= 1);
  if (remaining_shift == 0) {
    return shift_to_return;
  }
  // Must be exactly 1 bit remaining. Are we lucky, did we create a zero?
  if ((x & 1) == 0) {
    x >>= 1;
  } else {
    y >>= 1;
  }
  return shift_to_return;
}

DyadicFraction& DyadicFraction::mult(UInt x) {
  if (m_value == 0 || x == 1) {
    return *this;
  }
  if (x == 0) {
    m_value = 0;
    m_exponent = 0;
    return *this;
  }
  auto current_bit_length = BitFunctions::get_bit_length(m_value);
  auto x_bit_length = BitFunctions::get_bit_length(x);

  // Of course, the bit length almost has the log property, but not quite:
  //    Bitlength(xy) <= Bitlength(x) + Bitlength(y),
  // but inequality can hold, i.e. 8,15, 8*15=120 < 128
  // have bitlengths 4,4,7 respectively.
  // We're just using the total bitlength to get an upper bound
  // on the size of the product.
  auto total_bit_length = current_bit_length + x_bit_length;
  if (total_bit_length <= 64) {
    m_value *= x;
    return *this;
  }
  // First, losslessly remove bits.
  const auto current_trailing_zeros =
      BitFunctions::get_number_of_rightmost_zero_bits(m_value);

  const auto x_trailing_zeros =
      BitFunctions::get_number_of_rightmost_zero_bits(x);

  current_bit_length -= current_trailing_zeros;
  m_value >>= current_trailing_zeros;
  x_bit_length -= x_trailing_zeros;
  x >>= x_trailing_zeros;

  // We've divided nx by 2^a.2^b, so we must compensate.
  // (But we haven't actually multiplied by x yet!)
  m_exponent += current_trailing_zeros + x_trailing_zeros;

  // We've losslessly removed as many bits as we can,
  // but maybe still not enough.

  if (m_value < x) {
    m_exponent +=
        remove_bits_lossy(m_value, current_bit_length, x, x_bit_length);
  } else {
    m_exponent +=
        remove_bits_lossy(x, x_bit_length, m_value, current_bit_length);
  }
  m_value *= x;
  return *this;
}

DyadicFraction& DyadicFraction::mult(const DyadicFraction& other) {
  mult(other.m_value);
  m_exponent += other.m_exponent;
  return *this;
}

DyadicFraction& DyadicFraction::mult_n_over_k(UInt n) {
  mult(n);
  // K = 1024 = 2^10.
  m_exponent -= 10;
  return *this;
}

double DyadicFraction::get_double() const {
  const double result =
      std::pow(2.0, m_exponent) * static_cast<double>(m_value);
  return result;
}

double DyadicFraction::get_log() const {
  return std::log(static_cast<double>(m_value)) + m_exponent * std::log(2.0);
}

std::string DyadicFraction::str() const {
  std::stringstream ss;
  ss << "val=" << m_value << " exp=" << m_exponent;
  return ss.str();
}

bool DyadicFraction::operator<(const DyadicFraction& other) const {
  // We're asking if   2^n.x < 2^m.y,   i.e.   x < 2^d.y ?
  const auto exponent_diff = other.m_exponent - m_exponent;
  if (exponent_diff >= 0) {
    // d >= 0.
    auto y = other.m_value;
    const unsigned y_bitlength = BitFunctions::get_bit_length(y);
    if (y_bitlength + (unsigned)exponent_diff > 64) {
      // 2^d.y would overflow, so definitely bigger than x.
      return true;
    }
    // We can shift without overflow.
    y <<= (unsigned)exponent_diff;
    return m_value < y;
  }
  // Now we are asking if 2^d.x < y for d>0.
  auto x = m_value;
  const unsigned x_shift = -exponent_diff;
  const unsigned x_bitlength = BitFunctions::get_bit_length(m_value);
  if (x_shift + x_bitlength > 64) {
    // 2^d.x overflows.
    return false;
  }
  x <<= x_shift;
  return x < other.m_value;
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket
