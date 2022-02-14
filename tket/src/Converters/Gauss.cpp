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

#include "Gauss.hpp"

namespace tket {

void CXMaker::row_add(unsigned r0, unsigned r1) {
  if (_reverse_cx_dirs)
    _circ.add_op<unsigned>(OpType::CX, {r1, r0});
  else
    _circ.add_op<unsigned>(OpType::CX, {r0, r1});
}

void DiagMatrix::row_add(unsigned r0, unsigned r1) {
  for (unsigned i = 0; i < _matrix.row(r0).size(); ++i) {
    _matrix.row(r1)[i] ^= _matrix.row(r0)[i];
  }
}

void DiagMatrix::col_add(unsigned c0, unsigned c1) {
  for (unsigned i = 0; i < _matrix.col(c0).size(); ++i) {
    _matrix.col(c1)[i] ^= _matrix.col(c0)[i];
  }
}

void DiagMatrix::gauss(CXMaker& cxmaker, unsigned blocksize) {
  std::vector<std::pair<unsigned, unsigned>> row_ops =
      gaussian_elimination_row_ops(_matrix, blocksize);
  for (std::pair<unsigned, unsigned> op : row_ops) {
    row_add(op.first, op.second);
    cxmaker.row_add(op.first, op.second);
  }
}

bool DiagMatrix::is_id() const { return _matrix.isIdentity(); }

bool DiagMatrix::is_id_until_columns(unsigned limit) const {
  TKET_ASSERT(limit <= n_rows());

  for (unsigned i = 0; i < n_rows(); ++i) {
    if (_matrix(i, i) == 0) return false;
  }

  for (unsigned i = 0; i < n_rows(); ++i) {
    for (unsigned j = 0; j < n_cols(); ++j) {
      if ((i > j) && (_matrix(i, j) == 1)) return false;
    }
  }

  for (unsigned i = 0; i < n_rows(); ++i) {
    for (unsigned j = 0; j < n_cols(); ++j) {
      if (j > limit) {
        if ((i < j) && (_matrix(i, j) == 1)) return false;
      }
    }
  }

  return true;
}

unsigned DiagMatrix::n_rows() const { return _matrix.rows(); }

unsigned DiagMatrix::n_cols() const { return _matrix.cols(); }

std::ostream& operator<<(std::ostream& out, const DiagMatrix& diam) {
  out << "give the DiagMatrix: " << std::endl;
  for (unsigned i = 0; i < diam._matrix.row(0).size(); ++i) {
    out << " ";
    for (unsigned j = 0; j < diam._matrix.row(0).size(); ++j) {
      out << diam._matrix.row(i)[j] << ", ";
    }
    out << std::endl;
  }
  out << std::endl;

  return out;
}

}  // namespace tket
