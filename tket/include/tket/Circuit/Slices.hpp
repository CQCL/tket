// Copyright 2019-2024 Cambridge Quantum Computing
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

#include <memory>

#include "tket/Circuit/DAGDefs.hpp"
#include "tket/Utils/SequencedContainers.hpp"

namespace tket {

typedef VertexVec Slice;

typedef sequenced_map_t<UnitID, Edge> unit_frontier_t;
// typedef sequenced_map_t<Qubit, Edge> q_frontier_t;
typedef sequenced_map_t<Bit, EdgeVec> b_frontier_t;
// typedef sequenced_map_t<Bit, Edge> c_frontier_t;

struct CutFrontier {
  std::shared_ptr<Slice> slice;
  std::shared_ptr<unit_frontier_t> u_frontier;
  std::shared_ptr<b_frontier_t> b_frontier;
  void init() {
    slice = std::make_shared<Slice>();
    u_frontier = std::make_shared<unit_frontier_t>();
    b_frontier = std::make_shared<b_frontier_t>();
  }
};

class Circuit;

/*SliceIterator class is used for lazy evaluation of slices */
class SliceIterator {
 public:  // these are currently public to allow skip_func slicing.
  CutFrontier cut_;
  std::shared_ptr<b_frontier_t> prev_b_frontier_;
  const Circuit *circ_;

  class Sliceholder {
   private:
    Slice current_slice_;

   public:
    explicit Sliceholder(Slice slice) : current_slice_(slice) {}
    Slice operator*() const { return current_slice_; }
  };

  // take in an unsigned 'n' and a circuit and give the 'n'th slice
  // note: n=0 gives an empty SliceIterator
  // n=1 gives the first slice

  SliceIterator(
      const Circuit &circ, const std::function<bool(Op_ptr)> &skip_func);
  explicit SliceIterator(const Circuit &circ);
  SliceIterator() : cut_(), circ_() { cut_.init(); }
  Slice operator*() const { return *cut_.slice; }
  bool operator==(const SliceIterator &other) const {
    return *cut_.slice == *other.cut_.slice;
  }
  bool operator!=(const SliceIterator &other) const {
    return !(*this == other);
  }
  std::shared_ptr<const unit_frontier_t> get_u_frontier() const {
    return cut_.u_frontier;
  }
  std::shared_ptr<const b_frontier_t> get_b_frontier() const {
    return cut_.b_frontier;
  }
  std::shared_ptr<const b_frontier_t> get_prev_b_frontier() const {
    return prev_b_frontier_;
  }
  // A postfix increment operator overload
  Sliceholder operator++(int);
  // A prefix increment operator overload
  SliceIterator &operator++();
  bool finished() const;
};

}  // namespace tket
