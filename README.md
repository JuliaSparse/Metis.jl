#  [Julia](http://julialang.org) interface to the Metis graph-partitioning algorithms

[![Build Status](https://travis-ci.org/JuliaSparse/Metis.jl.svg?branch=master)](https://travis-ci.org/JuliaSparse/Metis.jl)
[![Coverage Status](https://coveralls.io/repos/dmbates/Metis.jl/badge.png?branch=master)](https://coveralls.io/r/dmbates/Metis.jl?branch=master)
[![Metis](http://pkg.julialang.org/badges/Metis_release.svg)](http://pkg.julialang.org/?pkg=Metis&ver=release)

## Installation

```julia
Pkg.add("Metis")
```

Adding the package will install the Metis library itself on OS-X, Windows and Linux systems using `apt` or `yum`.

On other operating systems this package will download, configure and install metis-5.1.0 in the directory
```julia
Pkg.dir("Metis", "deps")
```
Configuration requires Cmake version 2.8 or later.

## Functions

The names of the available Julia functions are those from the Metis API

* `nodeND(al)` : recursively bisect the undirected graph implied by the adjacency list, 
  `al`, and return the fill-reducing permutation that results
* `nodeND(m)` : recursively bisect the undirected graph of the nonzero structure of the 
  symmetric sparse matrix, `m`, and return the fill-reducing permutation
* `vertexSep(al)`: compute a vertex separator for the adjacency list, `al`, of an
  undirected graph
* `vertexSep(m)`: compute a vertex separator for the undirected graph of the nonzero 
  structure of the symmetric sparse matrix, `m`
* `partGraphKway(al, nparts::Integer)`: partition a graph given as an adjacency list, `al`, into nparts 
* `partGraphRecursive(al, nparts::Integer)`: partition a graph given as an adjacency list, `al`, into nparts 

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
sizes, part = vertexSep(copter2)
objval, part = partGraphKway(copter2, 6)
counts = zeros(Int, 6);
for p in part counts[p] += 1 end
println(counts)
```
