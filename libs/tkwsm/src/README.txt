// Copyright Quantinuum
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

Version 0.41.

BASIC PROBLEM:

We are interested in the standard subgraph monomorphism problem, also known as the non-induced subgraph isomorphism problem:

Given graphs P, T (the pattern graph and target graph), find a subgraph monomorphism f from P to T. Thus, by definition f is an injective function

    f : V(P) -> V(T)

(i.e., f maps the vertices of P into the vertices of T, without clashes), such that if (v1,v2) are adjacent in P, then (f(v1), f(v2)) are adjacent in T.

Thus we can define (i.e., f "induces") a new function  F : E(P) -> E(T)  (i.e., mapping pattern edges to target edges) via  F((v1,v2)) = (f(v1), f(v2)), which is uniquely defined, well-defined, and injective, precisely because of the monomorphism property of f. Note that:

- f is "non-induced" because we allow (f(v1), f(v2)) to be adjacent in T even if (v1,v2) are not adjacent in P.

- We usually abuse notation and write f instead of F.

- Our graphs are always undirected, without loops or multiple edges.

In applications to qubit routing and placement problems, P is usually the interaction graph between the logical qubits (the qubits of a quantum circuit), and T is the interaction graph of an actual quantum computer with physical qubits.

We solve a new problem, which we call Weighted Subgraph Monomorphism (WSM): each pattern edge and target edge additionally has a weight w >= 0 (which is an integer). For each subgraph monomorphism f, let

    S(f) = sum_(e in E(P)) w(e).w(f(e))

We call S(f) the "scalar product" of a monomorphism f. Then, the WSM problem is to find f which minimises S(f).

(It seems unlikely that this is really new, but we have not seen this in existing papers).

The point is that we can adjust the P-edge weights to account for less "important" qubit interactions: i.e., less often used, or occurring further in the future of the circuit.

We can also adjust the T-edge weights to account for different 2-qubit gate fidelities on the hardware, and also add extra target edges to account for the fact that it is unlikely that an f exists mapping EVERY P-edge to a T-edge on the hardware device: thus we expect that we will have to apply extra SWAP gates to move the qubits around; thus we still want adjacent P-vertices to be mapped to nearby T-vertices, even if not exactly adjacent.

(So, as time progresses, the locations of the logical qubits, i.e. the physical qubits to which they have been assigned, are expected to "drift" away from the initial positions, and so we might adjust the P-edge weights with some kind of "time decay", to take account of the fact that far future assignments are not really worth taking into account very much, because the qubit locations may be quite different by then).


MAIN FEATURES:

- We can pass in a timeout (in milliseconds), and also an exact max number of iterations (useful for reproducible tests).

- If it terminates due to timeout (or max iterations) and we retain the solver object, we can resume solving later.

- We can set it to terminate immediately upon finding a solution (especially useful for unweighted problems, i.e. standard subgraph monomorphism problems).

- If it terminates without timeout (or hitting max iterations), it is GUARANTEED to have found an optimal solution, OR proved that none exists (unless it was set to terminate early upon finding any solution, and this did occur).


NOT YET IMPLEMENTED:

- (Trivial): pass in reduced domains, i.e. BEGIN by constraining Dom(x1) = {v1,v2,v3,...}, ... for some x1, x2, ...

- (Easy): implementing various bespoke constraints, like  Cost(f) <= K, where Cost(f) is an expression like   Cost(f) = sum_i m_i.Distance(y_i, f(v_i))  or similar, with given weights m_i > 0 and vertices y_i.

The exact form of Cost(F) is not very relevant; since the search tree goes through all possibilities, basically ANY cost function is possible. But it would work much better for functions like the above with a monotonic property, i.e. if F is a partial mapping and F' extends F, then Cost(F') >= Cost(F). This would enable pruning of the search tree.

- (Easy): if PV->TV, and we know Dom(x) for all neighbours x of PV, we can use maximum weighted bipartite matching to get a guaranteed lower bound on the contribution to the scalar product of all edges containing PV. In some cases, the optimal matching can prove that f(PV)=TV is impossible. This would work well with the WeightNogoodDetector.

  We tried a quick crude version of this, but it was too slow (even though finding a matching is fast, there were just too many of them). Something like the WeightNogoodDetectorManager will probably work well (i.e., we observe how successful the matchings are at detecting impossible assignments, or finding large lower bounds, and dynamically adjust the probabilities up or down as searching progresses).

- (Medium): if the target graph is complete and the pattern graph is a cycle, then WSM is exactly the Travelling Salesman Problem. Thus, for general pattern graphs but complete target graphs, we should be able to get reasonable solutions using simulated annealing and similar algorithms. Notice that the general WSM routines are NOT well suited to complete target graphs (because, they are basically brute force: no graph theoretic reductions are possible. Weight nogood detection is the only possible help here, but it seems a little too crude to work well).

  However, for the application to qubit placement, we don't need a really good solution; we propose finding an OK solution very quickly using trivial greedy heuristics, jumping about a bit with simulated annealing, and then ERASING some higher weight target edges and applying the normal WSM.


- (Hard): for approximate results (or maybe to consider restarts intelligently), it would be good to get a good estimator function for when the current solution is likely to be close to optimal. (E.g., there seem to be some problems where the solver finds a solution which is close to optimal or even optimal quite quickly, but takes much longer to attain or PROVE optimality - which is probably not needed for our applications). Maybe something fancy with machine learning could be used? There are several thousand test problems (and it is easy to generate many more). We could run the WSM solver through to completion on many classes of problems of varying sizes, and record the progress of the search tree (the shape, depth, weights, etc. of each branch we generate, in order, together with the total and partial scalar products). Maybe some useful, general numerical patterns will emerge?


FURTHER COMMENTS:

- Some of the basic ideas are similar to the Glasgow Subgraph Solver, described for example in the paper "A Parallel, Backjumping Subgraph Isomorphism Algorithm using Supplemental Graphs" by Ciaran McCreesh and Patrick Prosser. This paper definitely should be read.

- We do a backtracking search, using some simple graph theory to decide if a p-vertex could be mapped to a t-vertex, and simple heuristics to guide the search.

- We tried parallel searching, restarts, nogood recording as recommended in the Glasgow Subgraph Solver paper, but it did not appear to be useful for us.

(However, that paper only considers the unweighted problem, and the case where we stop at the FIRST solution. Thus, maybe it is still worth trying when we only need one solution. Also, the test problems in the paper are probably a lot larger than ours).

- Runtime seems to vary enormously from problem to problem, and merely counting vertices and edges is not enough to get accurate estimates of runtimes.

- It is often much faster if you only want a reasonable solution, or if you don't care about proving optimality. The solver can be stopped and restarted.

- Don't think of the list of nodes as a historical record of previous states (although it looks that way). Instead, think of them as future candidate states to collapse further (by making more assignments). E.g., when we decide to make a new assignment PV->TV, we remove TV from Domain(PV)={TV,a,b,c,...} before creating the next node, to leave Domain(PV) = {TV} in the new node, but Domain(PV)={a,b,c,...} in the old node. Thus the OLD node, now changed, is NOT any previously considered state, since Domain(PV) is different from both the original node AND the new node.


NOTABLE CHANGES FROM PREVIOUS VERSIONS:

Version 0.4 -> 0.41:
- A refactor has taken place: now, we use boost::dynamic_bitset to represent domains, rather than std::set<VertexWSM>. There are a few places where the old std::set<VertexWSM> code is still there (with ugly added temporary conversion code), which should be finished off; although, not as important, as none in critical loops. It seems like practically every solve is faster; depending on the problems, the speed improvement factor is 2-10 approximately (often around 3).

Version 0.31 -> 0.4:
- The InitPlacement code has been added. This is like an initial proof-of-concept showing how to use a WSM formulation to get reasonable initial qubit placements in a real scenario. Preliminary results are very encouraging but far more testing and experimentation should be done.

Version 0.3 -> 0.31:
- Refactor: internally relabels vertices if necessary (and converts them back again at the end - no visible difference for the user), to make them lie in {0,1,2,...,N}. This allows std::vector instead of std::map in some places. This simple change increases speed by >10% !

Version 0.2 -> 0.3:

- We used to have a std::vector of Node objects, with each Node containing a   std::map<vertex, std::set<vertex>>   object. Thus, for each PV, we have a set Dom(PV) which is the collection of all target vertices which PV could take (if not yet assigned). This corresponds exactly to the obvious mathematical formulation of the algorithm. HOWEVER, this requires a lot of copying of data. Instead, we now split into  std::map<vertex, history>,  where a "history" object records instead all versions of Domain(PV), for a SINGLE PV. We now only make a new domain object when it CHANGES. [REMARK: in version 0.31, these std::maps have become std::vectors!]

- When deciding on a new candidate p-vertex to assign (variable ordering), version 0.2 gave priority to those adjacent to newly assigned PV. However, removing this and simply treating all unassigned vertices equally was equally fast or significantly faster in almost every test case (even ignoring the extra slight speedup from erasing the extra bookkeeping code). Maybe this is not so surprising; initially we think of "nearby" vertices as being "more constrained", and so more natural to assign first. But we have many graph-theoretic reduction methods, some of which are more "long range". So, viewed purely as an abstract constrained variables problem, there is no compelling reason why we should pay attention to adjacent vertices.

- Version 0.2 separated the domains from assignments, i.e. we stored the main x -> Dom(x)  data in std::maps, but if  Dom(x)={y}  then x was not included, instead occurring in a separate  x->y map.  Now we've simplified and just have one main std::map to represent  x -> Dom(x).

- When searching and pruning the domains etc., we distinguish between a CHECK (which takes a single assignment PV->TV and sees if it is obviously impossible, independent of other domains) and a REDUCER (which reduces domains of OTHER nearby pattern vertices).

  In most cases, a reducer gives a corresponding check. E.g., the most obvious one is the NeighbourReducer: if  f(x)=y  has been set, let v be a neighbour of x (in the pattern graph). Then obviously  f(v)  must be a neighbour of  y  (in the target graph). Hence we can INTERSECT Dom(v) with {t in V(T) : t is adjacent to y}.

  Notice that the corresponding check is simply: #{neighbours of x} <= #{neighbours of y}.

  Some reducers don't quite fit into this framework, e.g. Hall set reduction: if, e.g., Dom(x1) = {y,z}, Dom(x2) = {y,z}, Dom(x3) = {s,t,y,z}, it is clear that only x1,x2 can take the values y,z. Thus we can immediately erase them from all other domains, so that Dom(x3) = {s,t}. This, of course, is a generalisation of the alldiff constraint  f(x)=y  =>  y is not in Dom(v) for all v!=x.

- Derived graphs (similar, but sometimes different, to "supplemental graphs" in the Glasgow Subgraph Solver paper): there are various transformations which transform a graph  (V,E)  into a new graph  (V,E'),  i.e. we keep the vertices but change the edges, sometimes quite drastically; and sometimes adding extra data, e.g. edge weights to E'.

  For some transformations, it has the property that if   f : V(P) -> V(T)  is a subgraph monomorphism from (V(P), E(P)) to (V(T), E(T)), then the same f is a subgraph monomorphism from (V(P), E'(P)) to (V(T), E'(T)).

  The most obvious example is the D2 derived graph (this is just my notation, I don't know of any standard notation):  for u,v in V, let  (u,v) in E'  <==>  there is a path [u,x,v] of length 2 (i.e., u,x,v are all distinct, and (u,x), (x,v) are in E).

  Notice that this is similar but different to dist(u,v)=2. Furthermore, one can give (u,v) in E' a weight (completely unrelated to the edge weights in the original graph): simply let it be the number of such distinct paths; then  f : E'(P) -> E'(T)  must increase edge weights (in addition to being well-defined).

  Similar things can be done with D3, D4, etc. but after D3 the computation becomes quite heavy. Furthermore, this raises the possibility of recursive deriving: one could take D2(D2), D2(D3), etc. etc. and do simple checks on f for a whole sequence of them.

  It sounds like some interesting theory could be worked out. However in practice these higher order derived graphs don't seem to be worth the computation: they seemed to take longer to compute than the time saved by using them for reducing.


POSSIBLE FUTURE PLANS/REFACTORINGS:

- Obviously, encoding a domain with a dynamic bitset instead of a std::set could be worth doing.

But it's not necessarily better unless we have hundreds or thousands of vertices, and possibly quite large domains. How long does it take to copy a dynamic bitset with N bits, vs. a std::set? If N is large, say N~1000, but the std::set has size << N, say size~5, the std::set may actually be faster. In many problems, as the search nears an end, many domains have size ~2 or so. Thus it's not so clear cut. Would have to implement it to see what difference it makes.

Of course, if we could GUARANTEE that we had <= 64 target vertices, we could use instead a std::uint64_t. This surely WOULD be faster.

Thus, a highly optimised version would abstract away all the interfaces, and have 3 possible specialisations: std::set (as now), dynamic bitset, and std::uint64_t, and choose one at runtime based upon sizes. Quite a lot of work, though, for maybe not much benefit.

