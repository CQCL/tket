// Copyright 2019-2023 Cambridge Quantum Computing
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

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <numeric>
#include <set>
#include <sstream>
#include <tkrng/RNG.hpp>

using std::stringstream;
using std::vector;

namespace tket {

// Check that the RNG really is identical across all platforms.

SCENARIO("RNG: test get_size_t") {
  RNG rng;
  stringstream ss;
  {
    std::array<size_t, 5> counts;
    counts.fill(0);
    const size_t max_v = counts.size() - 1;

    // Crudely check that we don't have too much bias.
    for (int nn = 0; nn < 100000; ++nn) {
      ++counts[rng.get_size_t(max_v)];
    }
    ss << "[ Counts for v=" << max_v << " : ";
    for (auto cc : counts) {
      ss << cc << " ";
    }
  }
  {
    const size_t max_v = 99;
    ss << "Values for v=" << max_v << " : ";
    for (int nn = 0; nn < 30; ++nn) {
      ss << rng.get_size_t(max_v) << " ";
    }
  }
  {
    const size_t min_v = 100;
    const size_t max_v = 105;
    ss << "Values for min_v=" << min_v << ", max_v=" << max_v << " : ";
    for (int nn = 0; nn < 20; ++nn) {
      ss << rng.get_size_t(min_v, max_v) << " ";
    }
    ss << "]";
  }
  CHECK(
      ss.str() ==
      "[ Counts for v=4 : 19878 19996 19936 20230 19960 "
      "Values for v=99 : 41 49 24 48 92 58 47 58 15 94 25 53 30 28 81 80"
      " 54 19 75 1 60 88 20 90 21 33 36 48 84 30 Values for min_v=100,"
      " max_v=105 : 104 100 104 103 104 105 102 104 100 100 105 101 101"
      " 104 104 104 104 103 102 103 ]");
}

SCENARIO("RNG: check_percentage bool sequence") {
  RNG rng;
  rng.set_seed(11111);
  stringstream ss;
  size_t number_of_true = 0;
  for (int nn = 0; nn < 100; ++nn) {
    if (rng.check_percentage(30)) {
      ++number_of_true;
      ss << "1";
    } else {
      ss << "0";
    }
  }
  CHECK(
      ss.str() ==
      "10100000000001100110100000101100100010000100000011100"
      "10000100011110000100000100000100000111000001000");
  CHECK(number_of_true == 29);
}

SCENARIO("RNG: vector operations") {
  RNG rng;
  rng.set_seed(22222);
  stringstream ss;
  vector<char> letters(26);
  std::iota(letters.begin(), letters.end(), 'a');
  for (int nn = 0; nn < 10; ++nn) {
    ss << rng.get_element(letters);
  }
  ss << "@";

  auto print_vector = [&ss, &letters]() {
    for (const char& ch : letters) {
      ss << ch;
    }
    ss << "#";
  };
  print_vector();
  rng.do_shuffle(letters);
  print_vector();
  while (!letters.empty()) {
    ss << rng.get_and_remove_element(letters);
  }
  CHECK(
      ss.str() ==
      "csifqddrqs@abcdefghijklmnopqrstuvwxyz#"
      "ejhdckamvzpbfsuirxonlgwtqy#tnljsfdvbgyhimpwxcourzkqea");
}

SCENARIO("RNG: permutations") {
  RNG rng;
  const size_t size = 100;
  const auto numbers = rng.get_permutation(size);
  REQUIRE(numbers.size() == size);

  const std::set<size_t> numbers_again(numbers.cbegin(), numbers.cend());
  REQUIRE(numbers_again.size() == size);
  REQUIRE(*numbers_again.crbegin() == size - 1);
  stringstream ss;
  ss << "[ ";
  for (auto x : numbers) {
    ss << x << " ";
  }
  ss << "]";
  REQUIRE(
      ss.str() ==
      "[ 37 15 43 14 44 5 21 29 55 35 28 9 13 32 56 3 30 40"
      " 16 22 20 33 2 79 58 38 17 11 47 73 46 51 1 61 74 0 49 6 75 39 53 19"
      " 62 18 60 93 87 8 59 89 78 96 25 10 27 76 70 72 80 99 86 77 91 63 50"
      " 66 82 88 83 7 12 31 41 90 54 67 57 85 42 48 4 94 81 92 23 34 95 26"
      " 69 24 68 71 64 84 36 65 97 98 52 45 ]");
}

SCENARIO("RNG: default seed") {
  vector<size_t> numbers;
  {
    RNG rng;
    for (unsigned ii = 0; ii < 10; ++ii) {
      numbers.emplace_back(rng.get_size_t(100));
    }
  }
  vector<size_t> numbers_again;
  {
    RNG rng;
    rng.set_seed();
    for (unsigned ii = 0; ii < 10; ++ii) {
      numbers_again.emplace_back(rng.get_size_t(100));
    }
  }
  CHECK(numbers == numbers_again);
  CHECK(numbers == vector<size_t>{79, 25, 71, 95, 1, 40, 25, 2, 52, 34});
}

SCENARIO("RNG: 64-bit uints") {
  RNG rng;
  // Check some fixed 4-bit blocks for approximate uniformity.
  // Count how often each of the 16 possible bits occurs.
  // OK, this is not a very precise statistical test,
  // but OK for a quick sanity check.
  std::array<std::array<unsigned, 16>, 5> counts;
  for (auto& list : counts) {
    std::fill(list.begin(), list.end(), 0);
  }
  vector<std::uint64_t> numbers;
  const std::array<unsigned, 5> shift_amounts{0, 10, 25, 41, 60};

  for (unsigned ii = 0; ii < 10000; ++ii) {
    const std::uint64_t number = rng();
    if (ii < 10) {
      numbers.push_back(number);
    }
    for (unsigned jj = 0; jj < shift_amounts.size(); ++jj) {
      const auto block_value = (number >> shift_amounts[jj]) & 0xf;
      counts.at(jj).at(block_value) += 1;
    }
  }
  CHECK(
      numbers == vector<std::uint64_t>{
                     0xc96d191cf6f6aea6, 0x401f7ac78bc80f1c, 0xb5ee8cb6abe457f8,
                     0xf258d22d4db91392, 0x4eef2b4b5d860cc, 0x67a7aabe10d172d6,
                     0x40565d50e72b4021, 0x5d07b7d1e8de386, 0x8548dea130821acc,
                     0x583c502c832e0a3a});

  // The counts should be roughly equal. Looks OK at a glance.
  const std::array<std::string, 5> expected_count_strings{
      "[ 610 686 615 602 598 613 607 610 632 589 599 630 626 647 674 662 ]",
      "[ 603 577 623 674 619 638 630 682 628 614 653 597 644 596 644 578 ]",
      "[ 636 600 621 604 660 636 608 601 662 623 622 635 587 607 640 658 ]",
      "[ 584 619 614 622 640 607 594 639 599 647 660 659 661 619 616 620 ]",
      "[ 619 682 617 616 619 574 573 632 655 640 616 637 643 631 628 618 ]"};
  for (unsigned ii = 0; ii < counts.size(); ++ii) {
    std::stringstream ss;
    ss << "[ ";
    for (auto x : counts[ii]) {
      ss << x << " ";
    }
    ss << "]";
    CHECK(expected_count_strings[ii] == ss.str());
  }
}

}  // namespace tket
