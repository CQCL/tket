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

#pragma once

#include <variant>
#include <type_traits>

// Helper function to convert a value of one type to another
// Conversion must be possible to compile
template <typename TargetType, typename StartType>
TargetType convertTo(const StartType& value) {
    return static_cast<TargetType>(value);
}

// Function to convert a variant with value of any of its types
// to a value of its first type
// All other types must be convertible to first
template <typename First, typename... Other>
First convertVariantToFirstType(const std::variant<First, Other...>& var) {
    return std::visit([](const auto& value) {
        return convertTo<First>(value);
    }, var);
}

// Function to convert a vector of a variant with values of any of its types
// to a vector of values of its first type
// All other types must be convertible to first
template <typename First, typename... Other>
std::vector<First> convertVariantVectorToFirstTypeVector(const std::vector<std::variant<First, Other...>>& variant_vec){
    auto first_vec = std::vector<First>();
    first_vec.reserve(variant_vec.size());
    for (const auto& variant_scalar: variant_vec){
        first_vec.push_back(convertVariantToFirstType(variant_scalar));
    }
    return first_vec;
}
