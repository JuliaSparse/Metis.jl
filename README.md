#  [Julia](http://julialang.org) interface to the Metis graph-partitioning algorithms

[![Build Status](https://travis-ci.org/JuliaSparse/Metis.jl.svg?branch=master)](https://travis-ci.org/JuliaSparse/Metis.jl)
[![Coverage Status](https://coveralls.io/repos/dmbates/Metis.jl/badge.png?branch=master)](https://coveralls.io/r/dmbates/Metis.jl?branch=master)
[![Metis](http://pkg.julialang.org/badges/Metis_release.svg)](http://pkg.julialang.org/?pkg=Metis&ver=release)

## Installation

```julia
Pkg.add("Metis")
```

On systems without a pre-packaged Metis library, adding this package will download, configure and install metis-5.1.0 in the directory
```julia
Pkg.dir("Metis", "deps")
```
Configuration requires Cmake version 2.8 or later.  At present the build is only available on Linux/Unix.  Contributions of build stanzas for Windows or OS X are welcome.

## Functions

The names of the available Julia functions are those from the Metis API

`nodeND(al)`
: create a fill-reducing permutation from the adjacency list, `al`, of a symmetric sparse matrix
`nodeND(m)`
: create a fill-reducing permutation from a symmetric sparse matrix, `m`.
`nodeND!(m)`
: mutating version of `nodeND`
`partGraphKway(al, nparts::Integer)`
: partition a graph given as an adjacency list, `al`, into nparts 
`partGraphRecursive(al, nparts::Integer)`
: partition a graph given as an adjacency list, `al`, into nparts 

## Examples

The function `Metis.testgraph` can be used to read one of the sample graphs available in the `graphs` directory of `metis-5.1.0`.  These graphs correspond to 2D and 3D finite element meshes.

`4elt`
	: a smaller sample graph (15606 vertices, 45878 edges)
`copter2` 
	: a medium size sample graph (55476 vertices, 352238 edges)
`mdual`
	: a larger sample graph (258569 vertices, 513132 edges)

```julia
using Graphs, Metis
copter2 = Metis.testgraph("copter2");
perm, iperm = nodeND(copter2)
objval, part = partGraphKway(copter2, 6)
counts = zeros(Int, 6);
for p in part counts[p] += 1 end
println(counts)
```
