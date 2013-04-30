#  [Julia](http://julialang.org) interface to the Metis graph-partitioning algorithms

## Installation

```julia
Pkg.add("Metis")
```

Adding the package will download, configure and install metis-5.1.0 in the directory
```julia
joinpath(Pkg.dir(), "Metis", "deps")
```
Configuration requires Cmake version 2.8 or later.  At present the build is only available on Linux/Unix.  Contributions of build stanzas for Windows or OS X are welcome.

## Functions

The names of the available Julia functions are those from the Metis API

`nodeND(al)`
: create a fill-reducing permutation from the adjacency list, `al`, of a symmetric sparse matrix
`partGraphKWay(al, nparts::Integer)`
: partition a graph given as an adjacency list, `al`, into nparts 
`partGraphRecursive(al, nparts::Integer)`
: partition a graph given as an adjacency list, `al`, into nparts 

## Examples

The function `Metis.testgraph` can be used to read one of the sample graphs available in the `graphs` directory of `metis-5.1.0`.  These graphs correspond to 2D and 3D finite element meshes.

4elt
: a smaller sample graph (15606 vertices, 45878 edges)
copter2 
: a medium size sample graph (55476 vertices, 352238 edges)
mdual
: a larger sample graph (258569 vertices, 513132 edges)

```julia
using Graphs, Metis
copter2 = Metis.testgraph("copter2");
perm, iperm = nodeND(copter2)
objval, part = partGraphKWay(copter2, 6)
counts = zeros(Int, 6);
for p in part counts[p] += 1 end
println(counts)
```

## Obtaining a fill-reducing permutation for a sparse matrix

The function sparse2adjacencylist in the Graphs package creates and adjacency list representation of the pattern of a `SparseMatrixCSC` object.  To be used with nodeND the original sparse matrix must be Hermitian.
```julia
## if A is a sparse Hermitian matrix
perm, iperm = nodeND(sparse2adjacencylist(A)) 
## returns a fill-reducing permutation and its inverse
```
