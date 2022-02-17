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

Some brief explanation of the Token Swapping algorithms here may be helpful.

PROBLEM: let G be a graph (undirected, no loops or multiple edges), with labelled distinct tokens (counters) on some or all of the vertices.

- An allowed move is to choose two adjacent vertices v1, v2 and swap whatever tokens T1, T2 are currently on those vertices. (Or, if one of the vertices, say v2, does not have a token, simply move the token T1 from v1 to v2).

- We are given a desired final rearrangement of the tokens.

- The problem is to compute a swap sequence, of shortest length, which will transform the initial token configuration into the desired final configuration.

- Thus, if every vertex contains a token, we are trying to perform a given permutation on the vertices of the graph.

- It is not hard to show: that if the graph is connected, then EVERY rearrangement is possible.

- More generally, a solution exists if and only if: for every token, the initial and final vertex are in the same connected component.

The 2016 paper "Approximation and Hardness of Token Swapping" by Tillmann Miltzow, Lothar Narins, Yoshio Okamoto, Gunter Rote, Antonis Thomas, Takeaki Uno is very useful.

One of the main results is an algorithm, based upon finding cycles in the graph, which (in the full case, where every vertex has a token) is guaranteed to use no more than 4x the optimal number of swaps, or 2x for trees.

The 2019 paper "Circuit Transformations for Quantum Architectures" by Andrew M. Childs, Eddie Schoute and Cem M. Unsal generalises to the partial case (where some vertices might not contain a token).


KNOWN RESTRICTIONS:

If the graph G is not connected, our routines may fail.

- We plan to fix this; until then, a workaround is to split the problem into connected components.

- Of course, in real problems, architectures are connected so this doesn't arise.


OUR ALGORITHMS:

Let L be the sum of distances of each token T from its final destination vertex.

Thus L >= 0 always, and the problem is finished if and only if L=0.

Thus, the goal of a token swapping algorithm (TSA) is to reduce L to zero, using as few swaps as possible.

(1) Cycle finding

We make some observations:

- The algorithms in the 2016 and 2019 papers try to find a cycle [v(0), v(1), ..., v(n), v(0)] in G. Initially, let vertex v(j) have token T(j). The swaps [v(0), v(1)], [v(1), v(2)], ..., [v(n-1), v(n)] then have the effect of performing a cyclic shift: T(0) -> v(n), T(1) -> v(0), T(2) -> v(1), ..., T(n) -> v(n-1).

- However, we do not need the swap [v(n), v(0)]. Thus, we can perform cyclic shifts along paths [v(0), v(1), ..., v(n)] instead of cycles.

- The paper algorithms search for cycles where EVERY token move T->v is beneficial, i.e. moves the token T one step closer to its final destination. Thus, the cycle on n vertices reduces L by n, at the cost of n-1 swaps.

- Cyclic shifts on n vertices may exist which contain swaps which don't decrease L (they might even increase L). However, as long as the overall effect is a decrease in L, the cyclic shift is not bad.

Thus, our algorithm works by searching for good PATHS instead of cycles; it allows some individual swaps to be bad, as long as the overall effect is still good.

(2) Alternative moves when cycles fail

When no good cycles exist, a different move must be performed. The performed move might not decrease L (although L never increases, it may leave L unchanged), and so we must worry whether the algorithm will ever terminate.

- In the papers, the performed moves are carefully chosen so that theorems guarantee that the algorithm terminates in a reasonable number of moves.

- We have no such theoretical guarantee. Instead, we perform a "trivial TSA", i.e., a sequence of swaps guaranteed to reduce L to zero eventually, even if the number of swaps is not so good.

However, notice that we need not perform the "trivial TSA" swaps until the end; we can break off as soon as L decreases, and switch back to the normal cyclic-shift mode, which is expected to reduce L more quickly.

(3) Additional reductions

We perform two more reductions, which are new:

- General swap list optimisation: given a sequence S of swaps, various simple substitutions may reduce the length of S, whilst keeping the overall effect of S unchanged; provided that we only use the same swaps that were present in S, the new sequence is still valid (does not use any further edges of the graph).

- General table lookup reduction: we have a large precomputed table which contains optimal swap sequences on graphs with <= 6 vertices. Thus, given our computed swap sequence S, we find the vertex mapping between two times, look up an optimal swap sequence for the mapping in the table (using only edges in our given graph, i.e. valid swaps), and replace the swap segment if the new sequence is shorter.



THE MAIN ALGORITHMIC CLASSES:


BestFullTsa:

The main end-to-end mapping. Uses HybridTsa to compute a solution; then general reduction rules (SwapListOptimiser) and table lookups (SwapListSegmentOptimiser, SwapListTableOptimiser) to reduce the swap sequence length, whilst preserving the mapping.


HybridTsa:

Combines CyclesPartialTsa, TrivialTSA to get an overall full TSA (token swapping algorithm). This works by running CyclesPartialTsa whenever possible; switching to TrivialTSA when it gets stuck; and breaking off TrivialTSA and switching back to CyclesPartialTsa when progress is made (i.e., L is decreased).

Note that HybridTsa also uses simple heuristics about which abstract cycles to perform first, and which abstract swap to omit.

(I.e., to perform the abstract cyclic shift v0->v1->v2->...->v(n)->v0, we often write [v(n), v(n-1)], [v(n-1), v(n-2)], ..., [v2, v1], [v1, v0]. Thus we have OMITTED [v0, v(n)] from the set of swaps [v(i), v(i+1)] for i=0,1,2,...,n, where v(n+1) = v(0) by definition. However, we could have omitted any of the n+1 swaps).

HybridTsa tries to estimate which ordering is likely to reduce L the quickest; although TrivialTSA works, it is expected to be worse than CyclesPartialTsa when that class works, so we want to break off as soon as possible.


CyclesPartialTsa:

The main "cycle finding" class, corresponding to the cycles in the papers. "Partial" because it is not guaranteed to find a swap sequence. Constructs "concrete" paths, i.e. actual paths in the graph, which give rise to "abstract cycles" (i.e. a cyclic shift, performed by swapping along the path).


TrivialTSA:

Performs any desired vertex permutation, as follows:

(i) split the permutation into disjoint cycles. (Called "abstract" cycles, because they are unrelated to the cycles in the actual graph. They are TOTALLY UNRELATED to the cycles of "CyclesPartialTsa"! In fact they know nothing of the underlying graph).

(ii) decompose the abstract cycles into "abstract swaps", i.e. without knowing the edges of the graph, the cyclic shift v0->v1->v2->...->vn->v0 can be rewritten as the abstract swaps [vn, v(n-1)], [v(n-1), v(n-2)], ..., [v2, v1], [v1, v0], which might not be possible in the graph.

(iii) decompose the abstract swaps into concrete swaps. I.e., choose a path [u, v0, v1, ..., v(n), v] between given (u,v), so that the abstract swap(u,v) can be performed by swapping along the path.


RiverFlowPathFinder:

Actually computes the path required by TrivialTSA, for part (iii). We don't just choose a path at random; we deliberately make the paths overlap as much as possible, for better optimisation later (and tests showed that this really is significant).


