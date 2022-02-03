// Copyright 2019-2022 Cambridge Quantum Computing
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

#include <catch2/catch.hpp>
#include <set>
#include <sstream>

#include "Utils/Assert.hpp"

using Catch::Matchers::Contains;

// An assert function with abort obviously cannot be tested here;
// but we CAN test assert functions which only throw.
namespace tket {
namespace {
// Just ensure that we have checked every message.
class MessageChecker {
 public:
  explicit MessageChecker(const std::vector<std::string>& calc_messages)
      : m_ii_count(0), m_calc_messages(calc_messages) {}

  const std::string& get_message(int ii) {
    ++m_ii_count;
    m_values_of_ii_checked.insert(ii);
    return m_calc_messages.at(ii);
  }

  void final_checks() const {
    CHECK(m_values_of_ii_checked.size() == m_calc_messages.size());
    // the ii should be [0,1,2,...,m].
    CHECK(m_values_of_ii_checked.size() == m_ii_count);
    CHECK(*m_values_of_ii_checked.cbegin() == 0);
    CHECK(
        *m_values_of_ii_checked.crbegin() == m_values_of_ii_checked.size() - 1);
  }

 private:
  unsigned m_ii_count;
  const std::vector<std::string>& m_calc_messages;
  std::set<int> m_values_of_ii_checked;
};

}  // namespace

static void check_filename_is_included(
    const std::vector<std::string>& messages) {
  for (const auto& message : messages) {
    CHECK_THAT(message, Contains("test_TketAssertWithThrow.cpp"));
  }
}

static int get_number(int nn) {
  if (nn > 15) {
    throw std::runtime_error("Error!!");
  }
  return nn - 10;
}

SCENARIO("Simple asserts with throws") {
  std::vector<std::string> calc_messages;
  std::vector<int> values_of_nn_with_error;

  for (int nn = 0; nn <= 20; ++nn) {
    try {
      // Should throw for nn in [3,5]
      TKET_ASSERT_WITH_THROW((nn - 3) * (nn - 5) > 0);

      // Should throw for nn in [8,10]
      TKET_ASSERT_WITH_THROW(
          (nn - 8) * (nn - 10) > 0 || AssertMessage() << "N=" << nn);

      // Should throw for [16,20] (the function throws).
      TKET_ASSERT_WITH_THROW(get_number(nn) < 20);
    } catch (const std::exception& e) {
      values_of_nn_with_error.push_back(nn);
      std::stringstream ss;
      ss << "CHECK: nn=" << nn << " ; " << e.what();
      calc_messages.emplace_back(ss.str());
    }
  }

  CHECK(calc_messages.size() == 11);
  check_filename_is_included(calc_messages);

  MessageChecker checker(calc_messages);

  for (int ii = 0; ii <= 2; ++ii) {
    const auto& message = checker.get_message(ii);
    CHECK_THAT(
        message,
        Contains(std::string("CHECK: nn=") + std::to_string(ii + 3) + " ; "));
    CHECK_THAT(message, Contains("Assertion '(nn - 3) * (nn - 5) > 0'"));
  }
  for (int ii = 3; ii <= 5; ++ii) {
    const auto& message = checker.get_message(ii);
    const std::string n_value = std::to_string(ii + 5);
    CHECK_THAT(message, Contains(std::string("CHECK: nn=") + n_value + " ; "));
    CHECK_THAT(message, Contains("Assertion"));
    CHECK_THAT(message, Contains("failed:"));
    CHECK_THAT(message, Contains(std::string("'N=") + n_value + "'"));
  }
  for (int ii = 6; ii <= 10; ++ii) {
    const auto& message = checker.get_message(ii);
    const std::string n_value = std::to_string(ii + 10);
    CHECK_THAT(message, Contains(std::string("CHECK: nn=") + n_value + " ; "));
    CHECK_THAT(
        message,
        Contains("Evaluating assertion condition 'get_number(nn) < 20'"));
    CHECK_THAT(message, Contains("threw unexpected exception: 'Error!!'"));
  }
  CHECK(
      values_of_nn_with_error ==
      std::vector<int>{3, 4, 5, 8, 9, 10, 16, 17, 18, 19, 20});
  checker.final_checks();
}

// Throws for nn in [2,5] or [8,10] with message.
static int get_number_with_asserts(int nn) {
  TKET_ASSERT_WITH_THROW((nn - 2) * (nn - 5) > 0);

  TKET_ASSERT_WITH_THROW(
      (nn - 8) * (nn - 10) > 0 || AssertMessage() << "N=" << nn << ": second");

  return nn + 5;
}

SCENARIO("Asserts with throws within calls") {
  std::vector<std::string> calc_messages;
  std::vector<int> values_of_nn_with_error;
  for (int nn = 0; nn <= 30; ++nn) {
    try {
      // Throws for [2,5] or [8,10].
      const int mm = get_number_with_asserts(nn);

      // Throws for mm=15,16, so nn=10,11,
      // but NOT for 10 because of the above! So only for nn=11.
      TKET_ASSERT_WITH_THROW(!(mm >= 15 && mm <= 16));

      // Throws for [26,30], since mm=n+5.
      TKET_ASSERT_WITH_THROW(
          mm <= 30 || AssertMessage() << "N=" << nn << ", M=" << mm);

      // Should throw from nn-10, so [12,15] or [18,20] (with message).
      TKET_ASSERT_WITH_THROW(get_number_with_asserts(nn - 10) >= nn - 5);

      // Should throw from nn-15, so [17,20]
      // (except that [18,20] are covered above, so nn=17 only)
      // or [23,25].
      TKET_ASSERT_WITH_THROW(
          get_number_with_asserts(nn - 15) >= nn - 10 ||
          AssertMessage() << "assert with N=" << nn);
    } catch (const std::exception& e) {
      values_of_nn_with_error.push_back(nn);
      std::stringstream ss;
      ss << "CHECK: nn=" << nn << " ; " << e.what();
      calc_messages.emplace_back(ss.str());
    }
  }
  CHECK(calc_messages.size() == 24);
  check_filename_is_included(calc_messages);

  MessageChecker checker(calc_messages);

  for (int ii = 0; ii <= 3; ++ii) {
    const auto& message = checker.get_message(ii);
    CHECK_THAT(
        message,
        Contains(std::string("CHECK: nn=") + std::to_string(ii + 2) + " ; "));
    // comes from "get_number_with_asserts"
    CHECK_THAT(message, Contains("Assertion '(nn - 2) * (nn - 5) > 0'"));
    // the function name
    CHECK_THAT(message, Contains("get_number_with_asserts"));
  }
  for (int ii = 4; ii <= 6; ++ii) {
    const auto& message = checker.get_message(ii);
    const auto n_value = std::to_string(ii + 4);
    CHECK_THAT(message, Contains(std::string("CHECK: nn=") + n_value + " ; "));
    // comes from "get_number_with_asserts"
    CHECK_THAT(message, Contains("Assertion"));
    // the function name
    CHECK_THAT(message, Contains("get_number_with_asserts"));
    CHECK_THAT(message, Contains(std::string("'N=") + n_value + ": second'"));

    // comes from the second assert in the function, without a message.
    CHECK_THAT(message, !Contains("(nn - 2) * (nn - 5)"));
  }
  {
    const auto& message = checker.get_message(7);
    CHECK_THAT(message, Contains("CHECK: nn=11 ; "));
    CHECK_THAT(message, Contains("Assertion '!(mm >= 15 && mm <= 16)'"));

    CHECK_THAT(message, !Contains("get_number_with_asserts"));
  }
  for (int ii = 8; ii <= 11; ++ii) {
    const auto& message = checker.get_message(ii);
    const auto n_value = std::to_string(ii + 4);
    CHECK_THAT(message, Contains(std::string("CHECK: nn=") + n_value + " ; "));
    CHECK_THAT(
        message, Contains("Evaluating assertion condition "
                          "'get_number_with_asserts(nn - 10) >= nn - 5'"));
    CHECK_THAT(message, Contains("threw unexpected exception"));
    CHECK_THAT(message, Contains("Assertion '(nn - 2) * (nn - 5) > 0'"));

    CHECK_THAT(message, !Contains("AssertMessage()"));
  }
  {
    const auto& message = checker.get_message(12);
    CHECK_THAT(message, Contains(std::string("CHECK: nn=17 ; ")));
    CHECK_THAT(
        message, Contains("Evaluating assertion condition "
                          "'get_number_with_asserts(nn - 15) >= nn - 10 || "
                          "AssertMessage() << "));
    CHECK_THAT(message, Contains("threw unexpected exception"));
    CHECK_THAT(message, Contains("Assertion '(nn - 2) * (nn - 5) > 0'"));
  }
  for (int ii = 13; ii <= 15; ++ii) {
    const auto& message = checker.get_message(ii);
    CHECK_THAT(
        message,
        Contains(std::string("CHECK: nn=") + std::to_string(ii + 5) + " ; "));
    CHECK_THAT(
        message, Contains("Evaluating assertion condition "
                          "'get_number_with_asserts(nn - 10) >= nn - 5'"));
    CHECK_THAT(message, Contains("threw unexpected exception"));
    CHECK_THAT(message, Contains("Assertion"));
    CHECK_THAT(
        message,
        Contains(std::string("'N=") + std::to_string(ii - 5) + ": second"));

    CHECK_THAT(message, !Contains("(nn - 2) * (nn - 5)"));
    CHECK_THAT(message, !Contains("AssertMessage()"));
  }
  for (int ii = 16; ii <= 18; ++ii) {
    const auto& message = checker.get_message(ii);
    CHECK_THAT(
        message,
        Contains(std::string("CHECK: nn=") + std::to_string(ii + 7) + " ; "));
    CHECK_THAT(
        message,
        Contains(
            "Evaluating assertion condition "
            "'get_number_with_asserts(nn - 15) >= nn - 10 || AssertMessage()"));
    CHECK_THAT(message, Contains("threw unexpected exception"));
    CHECK_THAT(message, Contains("Assertion"));
    CHECK_THAT(
        message,
        Contains(std::string("'N=") + std::to_string(ii - 8) + ": second"));

    CHECK_THAT(message, !Contains("(nn - 2) * (nn - 5)"));
  }
  for (int ii = 19; ii <= 23; ++ii) {
    const auto& message = checker.get_message(ii);
    const auto n_value = std::to_string(ii + 7);
    CHECK_THAT(message, Contains(std::string("CHECK: nn=") + n_value + " ; "));
    CHECK_THAT(message, Contains("Assertion "));
    CHECK_THAT(message, Contains("failed: "));
    CHECK_THAT(
        message, Contains("'N=" + n_value + ", M=" + std::to_string(ii + 12)));

    CHECK_THAT(message, !Contains("Evaluating assertion condition"));
    CHECK_THAT(message, !Contains("get_number_with_asserts"));
    CHECK_THAT(message, !Contains("threw unexpected exception"));
    CHECK_THAT(message, !Contains("Assertion()"));
    CHECK_THAT(message, !Contains("(nn - 2) * (nn - 5)"));
  }
  CHECK(values_of_nn_with_error == std::vector<int>{2,  3,  4,  5,  8,  9,
                                                    10, 11, 12, 13, 14, 15,
                                                    17, 18, 19, 20, 23, 24,
                                                    25, 26, 27, 28, 29, 30});
  checker.final_checks();
}

SCENARIO("Asserts with various bool conversions") {
  // First, list things which do throw.
  bool throws = true;
  try {
    TKET_ASSERT_WITH_THROW(!"");
    throws = false;
  } catch (const std::exception&) {
  }
  CHECK(throws);

  throws = true;
  try {
    TKET_ASSERT_WITH_THROW(0);
    throws = false;
  } catch (const std::exception&) {
  }
  CHECK(throws);

  int xx = 1;
  try {
    // Now, list non-throwing things first.
    TKET_ASSERT_WITH_THROW("");
    ++xx;
    TKET_ASSERT_WITH_THROW("aaaaa");
    ++xx;
    TKET_ASSERT_WITH_THROW(xx);
    ++xx;
    TKET_ASSERT_WITH_THROW(true);
    ++xx;
    TKET_ASSERT_WITH_THROW(-1);
    ++xx;
    TKET_ASSERT_WITH_THROW(xx > 0);
    ++xx;
    // Throws
    TKET_ASSERT_WITH_THROW(!"bbbbb");
    xx *= 1000;
  } catch (const std::exception&) {
    xx *= 100;
  }
  CHECK(xx == 700);
}

}  // namespace tket
