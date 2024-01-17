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

#include "tket/ZX/Rewrite.hpp"

namespace tket {

namespace zx {

Rewrite Rewrite::to_graphlike_form() {
  return Rewrite::sequence(
      {Rewrite::rebase_to_zx(), Rewrite::red_to_green(),
       Rewrite::spider_fusion(), Rewrite::parallel_h_removal(),
       Rewrite::io_extension(), Rewrite::separate_boundaries()});
}

Rewrite Rewrite::reduce_graphlike_form() {
  Rewrite reduce = Rewrite::sequence(
      {Rewrite::repeat(Rewrite::remove_interior_cliffords()),
       Rewrite::extend_at_boundary_paulis(),
       Rewrite::repeat(Rewrite::remove_interior_paulis()),
       Rewrite::gadgetise_interior_paulis()});
  return Rewrite::sequence(
      {reduce, Rewrite::repeat_while(Rewrite::merge_gadgets(), reduce)});
}

Rewrite Rewrite::to_MBQC_diag() {
  return Rewrite::sequence(
      {Rewrite::rebase_to_mbqc(), Rewrite::extend_for_PX_outputs(),
       Rewrite::repeat(Rewrite::internalise_gadgets())});
}

}  // namespace zx

}  // namespace tket
