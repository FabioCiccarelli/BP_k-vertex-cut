# Solution Files

This directory contains the best incumbent solutions found by our branch-and-price algorithm for the k-vertex cut problem.

## Directory Structure

### `unweighted/`
Contains solution files for unweighted graph instances.
### `weighted/`
Contains solution files for weighted graph instances.

## Solution File Format

Each solution file follows the naming convention: `<instance_name>_k<value>.sol`

Where:
- `<instance_name>` is the name of the graph instance
- `<value>` is the k parameter (target number of connected components)

### File Content Structure

Each solution file contains exactly **two lines**:

1. **Line 1**: The cost of the vertex separator (cost of vertices removed)
2. **Line 2**: Space-separated list of vertex indices that form the separator

### Example

For a file named `karate_k5.sol`:
```
2
32 33
```

This means:
- The optimal solution removes **2 vertices** from the graph (cost=2)
- The vertices removed are **vertex 32** and **vertex 33**
- After removing these vertices, the graph has at least 5 connected components

### Vertex Indexing

- Vertices are **0-indexed** (starting from 0)
- Vertex indices do not correspond to the node numbering in the original DIMACS graph file (which is 1-indexed)


## Notes

- All solutions represent the best incumbent found during our computational experiments, on our hardware and with a time limit of one hour
- Solutions may be optimal or best-known feasible solutions depending on whether the instance was solved to optimality
- For optimality status and detailed computational statistics, refer to `BPresults.xlsx` in the parent `results/` directory
