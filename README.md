#  [Julia](http://julialang.org) interface to the Metis graph-partitioning algorithms

## Installation

```julia
Pkg.add("Metis")
```

## Obtaining a fill-reducing permutation for a sparse matrix

```julia
using Graphs, Metis
## if A is a sparse Hermitian matrix
perm, iperm = NodeND(sparse2adjacencylist(A)) 
## returns a fill-reducing permutation and its inverse
```
