# Code documentation

SCIP-based branch-and-price solver for the k-vertex-cut problem. 

## Quick start
- Run with a DIMACS format graph: `./bin/kvertexcut -f instance.dimacs` (loads `scip.set` by default).
- Load custom SCIP settings: `./bin/kvertexcut -f instance.dimacs -s params/params_k10.set`.

## Inputs
- **DIMACS graph file** (`.dimacs`): problem line `p edge nnodes nedges`, edges `e u v` (1-indexed nodes).
- **Optional weights**: if `options/weighted = 1`, a companion file `instance.dimacs.w` with one vertex weight per line (node order) is required. Defaults to weight 1 otherwise.
- **Problem parameter**: `params/k` (int, k >= 1), the target minimum number of vertex clusters in the residual graph.

## Parameters

### Required paths
- `output/results` (string): base directory for `.res` results file.
- `output/solution` (string): base directory for `.sol` solution file.
- `output/plot` (string): base directory for solution plots.

### Problem configuration
- `params/k` (int, default: 2): value of k for the k-vertex-cut problem.
- `options/weighted` (bool, default: 0): enable weighted variant (requires `.w` file).

### Algorithm options
- `options/symmetrymethod` (int, default: 2): symmetry handling method.
  - `0`: none
  - `1`: symresacks
  - `2`: lexicographic/orbital reduction
- `options/cliquegeneration` (int, default: 2): clique generation strategy.
  - `0`: edges only
  - `1`: edge-clique partitioning
  - `2`: edge-clique covering
- `options/relaxconstraints` (bool, default: 1): relax equality constraints to `>=` (see the paper for further details).
- `options/connectivity/cut` (bool, default: 1): enforce minimum connectivity cut.

### Warm start
- `options/connectivity/warmstart` (bool, default: 1): connectivity-based warm start.
- `options/ILPwarmstart` (bool, default: 0): ILP warm start.

### Solve mode
- `options/solveLP` (bool, default: 0): solve LP relaxation only (sets `limits/nodes = 1` and relaxes the integrality on binary variables).
- `options/plot` (bool, default: 0): enable PNG plot of the solution.

## Outputs

### Results file (`.res`)
Tab-separated file `<output/results>/<instance>_k<k>.res` (appends each run). Columns:

1. instance name
2. number of vertices in original graph
3. number of edges in original graph
4. number of vertices in residual graph (after preprocessing)
5. number of edges in residual graph
6. value of parameter k
7. symmetry method (0/1/2)
8. connectivity warmstart flag (0/1)
9. ILP warmstart flag (0/1)
10. connectivity cut flag (0/1)
11. clique generation method (0/1/2)
12. relax constraints flag (0/1)
13. LP-only mode flag (0/1)
14. weighted variant flag (0/1)
15. number of vertices fixed by preprocessing
16. total cost of prefixed vertices
17. solution status (Optimal, TimeLimit, Trivial, etc.)
18. best primal bound found
19. best dual bound found
20. total cost of the computed cut (including the prefixed vertices)
21. total solving time (seconds)
22. time spent in pricing (seconds)
23. number of variables (columns) generated
24. total branch-and-bound nodes explored
25. number of original variables
26. number of original constraints
27. columns generated at root node
28. columns generated via Farkas pricing
29. maximum branch-and-bound tree depth
30. dual bound at root node

### Solution file (`.sol`)
File `<output/solution>/<instance>_k<k>.sol` with solution details:
- **Line 1**: total cut cost (sum of costs of vertices in the cut).
- **Line 2**: space-separated list of vertex indices (1-indexed) in the cut.

Example:
```
42
3 7 12 15 20
```

### Plot file (`.png`)
If `options/plot = 1`, generates `<output/plot>/<instance>_k<k>.png` visualizing the solution:
- **Green nodes**: vertices not in the cut, grouped by connected component.
- **Red nodes**: vertices in the cut (separator).
- **Layout**: largest component at center, smaller components orbiting around, cut vertices in a grid at the top.
- **Edges**: solid black (between non-cut nodes), dashed gray (incident to cut nodes).

Requires Graphviz. A `.dot` file is generated in any case if `options/plot = 0`.

## Notes
- The directories in `output/*` must exist before running; otherwise execution aborts. If they are provided,by default the result and solution files are saved in the current directory.
- With `options/solveLP = 1`, the node limit is forced to 1 to retrieve only the LP solution.
- For weighted runs, ensure the `.w` file is in the same directory of the `.dimacs` instance.
