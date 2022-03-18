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


BASIC PROBLEM:

We are interested in the standard subgraph monomorphism problem, also known as the non-induced subgraph isomorphism problem:

Given graphs P, T (the pattern graph and target graph), find a subgraph monomorphism f from P to T. Thus, f is an injective function

    f : V(P) -> V(T)

such that

    e in E(P) => f(e) in E(T).

In other words, f maps the vertices of P into the vertices of T, without clashes, such that if (v1,v2) are adjacent in P, then (f(v1), f(v2)) are adjacent in T.

(It is "non-induced" because we allow (f(v1), f(v2)) to be adjacent in T even if (v1,v2) are not adjacent in P).

Our graphs are always undirected, without loops or multiple edges.

In applications to qubit routing and placement problems, P is usually the interaction graph between the logical qubits (the qubits of a quantum circuit), and T is the interaction graph of an actual quantum computer with physical qubits.

We also consider a new problem, called Weighted Subgraph Monomorphism (WSM): each pattern edge and target edge additionally has a weight w >= 0 (which is an integer). For each subgraph monomorphism f, let

    S(f) = sum_(e in E(P)) w(e).w(f(e))

We call S(f) the "total weight" or "scalar product" of a monomorphism f.

Then, the WSM problem is to find f which minimises S(f).

(It seems unlikely that this is a new problem, but we have not found this in existing papers).

The point is that we can adjust the P-edge weights to account for less "important" qubit interactions: i.e., less often used, or occurring further in the future of the circuit.

We can also adjust the T-edge weights to account for different 2-qubit gate fidelities on the hardware, and also add extra target edges to account for the fact that it is unlikely that an f exists mapping EVERY P-edge to a T-edge on the hardware device: thus we expect that we will have to apply extra SWAP gates to move the qubits around; thus we still want adjacent P-vertices to be mapped to closeby T-vertices, even if not exactly adjacent.

(So, as time progresses, the locations of the logical qubits, i.e. the physical qubits to which they have been assigned, are expected to "drift" away from the initial positions, and so we adjust the P-edge weights with some kind of "time decay", to take account of the fact that far future assignments are not really worth taking into account very much, because the qubit locations may be quite different by then).


MAIN FEATURES:

- We can pass in a timeout (in milliseconds), and also an exact max number of iterations (useful for reproducible tests).

- We can pass in "suggestions", i.e. a list of assignments PV->TV which it will try to satisfy first.

- If it terminates due to timeout (or max iterations), we can still continue solving later.

- We can set it to terminate immediately upon finding a solution (especially useful for unweighted problems, i.e. standard subgraph monomorphism problems).

- If it terminates without timeout (or hitting max iterations), it is GUARANTEED to have found an optimal solution, OR proved that none exists (unless it was set to terminate early upon finding any solution, and this did occur).


NOT YET IMPLEMENTED:

- (Trivial): pass in reduced domains, i.e. BEGIN by constraining Dom(x1) = {v1,v2,v3,...}, ... for some x1, x2, ...

- (Easy): implementing various bespoke constraints, like  Cost(f) <= K, where Cost(f) is an expression like   Cost(f) = sum_i m_i.Distance(y_i, f(v_i))  or similar, with given weights m_i > 0 and vertices y_i.

The exact form of Cost(F) is not very relevant; since the search tree goes through all possibilities, basically ANY cost function is possible. But it would work much better for functions like the above with a monotonic property, i.e. if F is a partial mapping and F' extends F, then Cost(F') >= Cost(F). This would enable pruning of the search tree.


MORE COMMENTS ABOUT "SUGGESTIONS":

- The solver will try to satisfy all of them, and then proceed as if the variable and value choosers had chosen those at random.

- They can be bad or invalid, which can cause poor runtimes. Suggestions are not very useful for the full weighted case running to completion. But, this mainly appears to be because they don't reduce the time taken to prove optimality; they can often enable the optimal solution to be found more quickly, BUT then the solver still takes about the same time to prove that no better solution exists.


ALGORITHMIC COMMENTS:

Call this version 0.2.

- Some of the basic ideas are similar to the Glasgow Subgraph Solver, described for example in the paper "A Parallel, Backjumping Subgraph Isomorphism Algorithm using Supplemental Graphs" by Ciaran McCreesh and Patrick Prosser.

- We do a backtracking search, using some simple graph theory to decide if a p-vertex could be mapped to a t-vertex, and simple heuristics to guide the search.

- We tried parallel searching, restarts, nogood recording as recommended in the Glasgow Subgraph Solver paper, but it did not appear to be useful for us.

(However, that paper only considers the unweighted problem, and the case where we stop at the FIRST solution. Thus, maybe it is still worth trying when we only need one solution).

- The data structures (composed of std::maps, std::sets, etc.) are far from optimised at present; fancier data structures are planned once the basic algorithms and heuristics are nailed down.

