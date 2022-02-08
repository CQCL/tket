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

#include "GateNode.hpp"

#include "BitOperations.hpp"
#include "Utils/Assert.hpp"

/*
The following is intended to be helpful, but there are no guarantees
(that it IS helpful)...

Suppose that a circuit has n qubits, and a gate acts on k qubits [q0, q1, ...],
where n >= k. Let M be the (2^n)*(2^n) unitary matrix representing the gate.

Let U be the (2^k)*(2^k) unitary matrix which would represent the same gate
acting on qubits [0,1,2,...,k-1] of a k-qubit circuit.

How do we obtain M from U, q0, q1, ..., k, n? Use ILO-BE convention.

For a binary string 001010010... of length n, the vector
|001010010...> corresponds to an element of C^(2^n),
a 2^n-dimensional vector space, as follows:

Replace each 0 by e(0) = (1 0)^T, and each 1 by e(1) = (0 1)^T, in C^2.

Between these 2D vectors, insert tensor product symbols
(still in the same order, from left to right).

The resulting vector in C^(2^n) has exactly one row, say p, containing 1,
and 0 everywhere else; and therefore equals e(p) in C^(2^n),
for some 0 <= p < 2^n.

p is just the number of zeros above the entry 1. But, how to find p?

p is exactly the usual integer value of the original binary string 001010010...
used to create the vector with tensor products.
(This is easy to prove, by induction on n).

We need to find M|abcdefgh...z> for each length-n binary string abc...z.

Now, "a" corresponds to qubit 0, "b" to qubit 1, ... Mark the bits
in positions [q0, q1,...] to get abc@e@@h...z, using "@" to denote a marked bit.

M|abc@e@@h...z> comes from calculating U|###> on the bits @@@ alone,
where ### are the bits @@@, but permuted if necessary (because the gate acts on
qubits [q0, q1, q2, ...] which is not necessarily [0,1,2,...]).

We "tear up" the values of U|###> and copy/move them around
to fill up the (2^n)*(2^n) matrix.

In the final matrix M, there are 2^(n-k) copies of each nonzero entry u_{i,j} of
U, placed in various positions. The positions depend only on q0, q1, ...  and
i,j,k,n, not on the value u_{i,j}. There are no collisions. Thus, M is given by
a linear function of the entries u_{i,j}.

For each i,j with 0 <= i,j < k, define the (2^k)*(2^k) matrix D = D(i,j) by

D_{r,s} = 1 if (r,s)=(i,j), and 0 otherwise.

Thus, U = sum_{i,j} u_{i,j} D(i,j).

So now, if we fix (i,j), replace U by D(i,j), and allow all bits to vary,
it defines a linear transformation represented by a matrix M(i,j),
and then the final matrix M we desire is just M = sum_{i,j} u_{i,j} M(i,j).

[Note that the D matrices are not unitary any more, but that's OK;
the formula still makes sense, it is all just linear equations].

So, what happens with D? An alternative characterisation of
D = D(i,j) is:  D(e(j)) = e(i),  and D(e(t)) = 0 for all other t.

So now, fix (i,j). The only way for M acting on |***@**@*@@***>
to be nonzero is for the bits @@@ to be (a suitable permutation of)
the bits ### which represent j in binary (because, then |###> = e(j)).

Now D|###> = e(i). Hence, when M acts on the vector |***@**@*@@***>,
it changes the bits @@@ to (a suitable permutation of) the bits representing i.

This, then, is the algorithm:

-   For any binary string T_k = b0.b1.b2.... of length k,
    and any length n binary string S_n,
    let f(T_k, S_n) be the length n binary string
    obtained by overwriting the bit in position q(r) with b(r).
    (Counting positions left to right).

-   Let g be the reverse operation: given any length n binary string S_n,
    let g(S_n) = T_k be the unique length k binary string such that
    f(T_k, S_n) = S_n.
    (Of course, there are 2^{n-k} different S_n giving the same T_k).

-   Pick any nonzero element u_{i,j} of U. Thus 0 <= i,j < 2^k.
    Consider any fixed    x = S_n    for which g(S_n) = T_k
    represents the integer j.
    Let   y = f(T', x)   where T' is the binary string
    representing the integer i.
    The linear transformation x -> y which maps all other length n strings
    to zero (regarding them as e(p) in C^{2^n}, where p is exactly the integer
    value of the binary string) is a (2^n * 2^n) matrix W with 1 in exactly one
    position, in fact (y,x) exactly (converting the binary strings x,y to
    integers), because W(e(x)) = e(y).

-   To build up M, let the value u_{i,j} be copied to position (y,x).

-   However, we can vary those bits of S_n which are NOT in the positions
    [q0, q1, ...] to get different (y,x), but still with the same i,j
    (because the i,j depend only on the bits at [q0, q1, ...]).
    The x values are all distinct from each other (and so, too, are the y),
    so we get 2^{n-k} distinct positions (y,x) where we place the value u_{i,j}.

SPARSITY: U is (2^k)*(2^k). Let 0 < s <= 1 be the sparsity of U,
i.e. U has s.4^k nonzero entries. The final matrix M then has
2^{n-k}.s.4^k = s.2^{n+k} nonzero entries, so dividing by 4^n,
the sparsity of M is s.2^{k-n}. Thus, even if U is dense, in most applications
M becomes very sparse as n increases since k <= 2 for almost all gates.
(In the few gates where k > 2, usually U is sparse, so that s << 1.
e.g., CCX has an 8x8 matrix U, with 8 nonzero entries, so s=1/8).

This sparsity is essential. On an ordinary 2020 laptop,
for n=8, with 560 gates, the full unitary M is found in <1 second.
A single statevector for n=16, with 1120 gates, is found in ~3 seconds.
This is without any really fancy, clever optimisation.

*/

namespace tket {
namespace tket_sim {
namespace internal {

namespace {
/** There are k distinct qubits [q0, q1, q2, ...] chosen from the set
 * {0,1,2,...,n-1}. We want each binary string x of length k to have its bits
 * extracted and moved to the positions [q0, q1, q2, ...] within a larger binary
 * string of length n.
 *
 *  Thus 0 <= x < 2^k if we interpret x as an unsigned integer in the usual way.
 *
 *  Element[x] of the returned vector gives exactly the bits corresponding to x.
 */
struct LiftedBitsResult {
  // Element[x] tells us exactly the result of translating the bits of x.
  std::vector<SimUInt> translated_bits;

  // Simply the OR of all translated bits; has a "1" everywhere there is a slot
  // for a translated bit to appear.
  SimUInt translated_bits_mask;

  /** Set the data members, given the list of relevant qubits, and the full
   * number of qubits from which they are chosen.
   *  @param qubits The distinct qubit indices chosen from {0,1,2,...,n-1}.
   *  @param full_number_of_qubits The value n, the number of qubits in the full
   * set, from which "qubits" is extracted.
   */
  void set(const std::vector<unsigned>& qubits, unsigned full_number_of_qubits);
};

void LiftedBitsResult::set(
    const std::vector<unsigned>& qubits, unsigned full_number_of_qubits) {
  TKET_ASSERT(full_number_of_qubits >= qubits.size());
  TKET_ASSERT(full_number_of_qubits < 32);
  TKET_ASSERT(!qubits.empty());

  translated_bits.assign(get_matrix_size(qubits.size()), 0);

  // Build up the elements by ORing with shifted bits.
  translated_bits_mask = 0;
  SimUInt k_string_bit = 1;

  for (unsigned count = 0; count < qubits.size(); ++count) {
    TKET_ASSERT(full_number_of_qubits >= qubits[count] + 1);

    // This will be a bit within the length n string.
    SimUInt long_string_bit = 1;
    long_string_bit <<=
        (full_number_of_qubits - (qubits[qubits.size() - 1 - count] + 1));

    translated_bits_mask |= long_string_bit;

    // Go through the entries and add the new long string bit
    // when appropriate
    // (i.e., if "k_string_bit" occurs in the binary expansion of i).
    for (SimUInt ii = 0; ii < translated_bits.size(); ++ii) {
      // If "ii" has a 1 in the appropriate position,
      // insert a "1" in the length n string for this ii.

      // Trick: convert 0 to 0, but convert any nonzero number to 1.
      // This hopefully is faster by avoiding branching.
      const SimUInt has_bit = std::min(ii & k_string_bit, SimUInt(1));
      translated_bits[ii] |= has_bit * long_string_bit;
    }
    k_string_bit <<= 1;
  }
}
}  // namespace

static void set_lifted_triplets(
    const std::vector<TripletCd>& triplets, LiftedBitsResult& lifted_bits,
    std::vector<TripletCd>& lifted_triplets, ExpansionData& expansion_data,
    const std::vector<unsigned>& qubits, unsigned full_number_of_qubits) {
  // translated_bits[j] gives the bits within the length n binary string,
  // which correspond to the length k binary representation of j,
  // but permuted and moved around so as to fit in the "slots" for qubits
  // [q0, q1, q2,...].
  // Since U is unitary of size 2^k, no column or row can be all zeros,
  // so every index j will be needed at least once.
  // Hence, it is not a waste of time to compute all,
  // even if U is very sparse.
  lifted_bits.set(qubits, full_number_of_qubits);
  lifted_triplets.clear();

  expansion_data = get_expansion_data(
      lifted_bits.translated_bits_mask, full_number_of_qubits);

  const SimUInt free_bits_limit =
      get_matrix_size(full_number_of_qubits - qubits.size());
  TKET_ASSERT(free_bits_limit != 0 || !"Too many bits");

  for (SimUInt free_bits = 0; free_bits < free_bits_limit; ++free_bits) {
    const SimUInt expanded_free_bits =
        get_expanded_bits(expansion_data, free_bits);

    // Now go down the U matrix entries and copy them to the appropriate
    // positions.
    for (const auto& triplet : triplets) {
      lifted_triplets.emplace_back(
          lifted_bits.translated_bits.at(triplet.row()) | expanded_free_bits,
          lifted_bits.translated_bits.at(triplet.col()) | expanded_free_bits,
          triplet.value());
    }
  }
}

namespace {
// Contains data potentially of size roughly 2^n or larger,
// to avoid expensive memeory reallocation.
struct LargeWorkData {
  LiftedBitsResult lifted_bits;
  std::vector<TripletCd> lifted_triplets;
  SparseMatrixXcd sparse_matrix;
  ExpansionData expansion_data;

  static LargeWorkData& get_work_data();
};

LargeWorkData& LargeWorkData::get_work_data() {
  static LargeWorkData data;
  return data;
}
}  // namespace

void GateNode::apply_full_unitary(
    Eigen::MatrixXcd& matr, unsigned full_number_of_qubits) const {
  auto& work_data = LargeWorkData::get_work_data();

  set_lifted_triplets(
      triplets, work_data.lifted_bits, work_data.lifted_triplets,
      work_data.expansion_data, qubit_indices, full_number_of_qubits);

  work_data.sparse_matrix =
      get_sparse_square_matrix(work_data.lifted_triplets, matr.rows());

  matr = work_data.sparse_matrix * matr;
}

}  // namespace internal
}  // namespace tket_sim
}  // namespace tket
