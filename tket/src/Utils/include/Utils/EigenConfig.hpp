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

#pragma once

/**
 * @file
 * @brief Include this file rather than including the Eigen headers directly
 */

#include "Utils/Json.hpp"

#if !defined(_MSC_VER)
#pragma GCC diagnostic push

#if defined(__clang__)
#if __has_warning("-Wdeprecated-copy")
#pragma GCC diagnostic ignored "-Wdeprecated-copy"
#endif
#elif __GNUC__ >= 11
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#endif

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/MatrixFunctions>

#if !defined(_MSC_VER)
#pragma GCC diagnostic pop
#endif

namespace Eigen {

template <typename _Scalar, int _Rows, int _Cols>
void to_json(nlohmann::json& j, const Matrix<_Scalar, _Rows, _Cols>& matrix) {
  for (Index i = 0; i < matrix.rows(); ++i) {
    nlohmann::json row = nlohmann::json::array();
    for (Index j = 0; j < matrix.cols(); ++j) {
      row.push_back(matrix(i, j));
    }
    j.push_back(row);
  }
}

template <typename _Scalar, int _Rows, int _Cols>
void from_json(const nlohmann::json& j, Matrix<_Scalar, _Rows, _Cols>& matrix) {
  for (size_t i = 0; i < j.size(); ++i) {
    const auto& j_row = j.at(i);
    for (size_t j = 0; j < j_row.size(); ++j) {
      matrix(i, j) =
          j_row.at(j).get<typename Matrix<_Scalar, _Rows, _Cols>::Scalar>();
    }
  }
}

}  // namespace Eigen
