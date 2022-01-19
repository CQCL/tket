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

#include "Circuit/Circuit.hpp"
#include "Converters.hpp"

namespace tket {

class CXMaker {
 public:
  explicit CXMaker(unsigned qubits, bool reverse_cx_dirs = false)
      : _circ(qubits), _reverse_cx_dirs(reverse_cx_dirs) {}
  void row_add(unsigned r0, unsigned r1);
  Circuit _circ;
  bool _reverse_cx_dirs;
};

class DiagMatrix {
 public:
  DiagMatrix() {}
  explicit DiagMatrix(const MatrixXb& matrix) : _matrix(matrix) {}
  void row_add(unsigned r0, unsigned r1);
  void col_add(unsigned c0, unsigned c1);
  void gauss(CXMaker& cxmaker, unsigned blocksize = 6);
  friend std::ostream& operator<<(std::ostream& out, const DiagMatrix& diam);
  bool is_id() const;
  bool is_id_until_columns(unsigned limit) const;
  unsigned n_rows() const;
  unsigned n_cols() const;

  MatrixXb _matrix;
};

}  // namespace tket
