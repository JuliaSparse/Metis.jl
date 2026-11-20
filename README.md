# Metis

| **Build Status**                                                      |
|:--------------------------------------------------------------------- |
| [![][gh-actions-img]][gh-actions-url] [![][codecov-img]][codecov-url] |

*Metis.jl* is a Julia wrapper to the [Metis library][metis-url] which is a
library for partitioning unstructured graphs, partitioning meshes, and
computing fill-reducing orderings of sparse matrices.

## Graph partitioning
`Metis.partition` calculates graph partitions. As an example, here we partition
a small graph into two, three and four parts, and visualize the result:

| ![][partition2-url]     | ![][partition3-url]     | ![][partition4-url]     |
|:----------------------- |:----------------------- |:----------------------- |
| `Metis.partition(g, 2)` | `Metis.partition(g, 3)` | `Metis.partition(g, 4)` |

`Metis.partition` calls `METIS_PartGraphKway` or `METIS_PartGraphRecursive` from the Metis
C API, depending on the optional keyword argument `alg`:
 - `alg = :KWAY`:  multilevel k-way partitioning (`METIS_PartGraphKway`).
 - `alg = :RECURSIVE`:  multilevel recursive bisection (`METIS_PartGraphRecursive`).

### Edge weights
`Metis.partition` uses edge weights when the graph has them:
 - Graphs from SimpleWeightedGraphs.jl use their weights by default. The weights must be
   positive integers, otherwise `Metis.graph` throws an error. This requires Julia 1.9 or
   later, since the fallback for Julia versions without package extensions doesn't include
   SimpleWeightedGraphs support.
 - The values of a sparse matrix `S` are ignored by default. To use them as edge weights,
   call `Metis.partition(Metis.graph(S; weights = true), n)`. The values must be
   integer-valued, and `weights = true` isn't supported for `Hermitian` or `Symmetric`
   wrapped matrices.

`Metis.separator` and `Metis.permutation` don't use edge weights.

## Vertex separator
`Metis.separator` calculates a [vertex separator](https://en.wikipedia.org/wiki/Vertex_separator)
of a graph. `Metis.separator` calls `METIS_ComputeVertexSeparator` from the Metis C API.
As an example, here we calculate a vertex separator (green) of a small graph:

| ![][separator-url]   |
|:-------------------- |
| `Metis.separator(g)` |

`Metis.separator` returns a vector with one label per vertex: `1` and `2` for the two parts
and `3` for the separator. No edge joins part `1` and part `2`. A part or the separator can
be empty, in which case its label is missing from the result.

## Fill reducing permutation
`Metis.permutation` calculates the fill reducing permutation
for a sparse matrix. `Metis.permutation` calls `METIS_NodeND` from the Metis
C API. As an example, we calculate the fill reducing permutation
for a sparse matrix `S` originating from a typical (small) FEM problem, and
visualize the sparsity pattern for the original matrix and the permuted matrix:

```julia
perm, iperm = Metis.permutation(S)
```

| <pre>⠛⣤⢠⠄⠀⣌⠃⢠⠀⠐⠈⠀⠀⠀⠀⠉⠃⠀⠀⠀⠀⠀⠀⠀⠀⠘⠀⠂⠔⠀<br>⠀⠖⠻⣦⡅⠘⡁⠀⠀⠀⠀⠐⠀⠁⠀⢂⠀⠀⠠⠀⠀⠀⠁⢀⠀⢀⠀⠀⠄⢣<br>⡀⢤⣁⠉⠛⣤⡡⢀⠀⠂⠂⠀⠂⠃⢰⣀⠀⠔⠀⠀⠀⠀⠀⠀⠀⠀⠀⠄⠄⠀<br>⠉⣀⠁⠈⠁⢊⠱⢆⡰⠀⠈⠀⠀⠀⠀⢈⠉⡂⠀⠐⢀⡞⠐⠂⠀⠄⡀⠠⠂⠀<br>⢀⠀⠀⠀⠠⠀⠐⠊⠛⣤⡔⠘⠰⠒⠠⠀⡈⠀⠀⠀⠉⠉⠘⠂⠀⠀⠀⡐⢈⠀<br>⠂⠀⢀⠀⠈⠀⠂⠀⣐⠉⢑⣴⡉⡈⠁⡂⠒⠀⠁⢠⡄⠀⠐⠀⠠⠄⠀⠁⢀⡀<br>⠀⠀⠄⠀⠬⠀⠀⠀⢰⠂⡃⠨⣿⣿⡕⠂⠀⠨⠌⠈⠆⠀⠄⡀⠑⠀⠀⠘⠀⠀<br>⡄⠀⠠⢀⠐⢲⡀⢀⠀⠂⠡⠠⠱⠉⢱⢖⡀⠀⡈⠃⠀⠀⠀⢁⠄⢀⣐⠢⠀⠀<br>⠉⠀⠀⠀⢀⠄⠣⠠⠂⠈⠘⠀⡀⡀⠀⠈⠱⢆⣰⠠⠰⠐⠐⢀⠀⢀⢀⠀⠌⠀<br>⠀⠀⠀⠂⠀⠀⢀⠀⠀⠀⠁⣀⡂⠁⠦⠈⠐⡚⠱⢆⢀⢀⠡⠌⡀⡈⠸⠁⠂⠀<br>⠀⠀⠀⠀⠀⠀⣠⠴⡇⠀⠀⠉⠈⠁⠀⠀⢐⠂⠀⢐⣻⣾⠡⠀⠈⠀⠄⠀⡉⠄<br>⠀⠀⠁⢀⠀⠀⠰⠀⠲⠀⠐⠀⠀⠡⠄⢀⠐⢀⡁⠆⠁⠂⠱⢆⡀⣀⠠⠁⠉⠇<br>⣀⠀⠀⢀⠀⠀⠀⠄⠀⠀⠀⠆⠑⠀⠀⢁⠀⢀⡀⠨⠂⠀⠀⢨⠿⢇⠀⡸⠀⢀<br>⠠⠀⠀⠀⠀⠄⠀⡈⢀⠠⠄⠀⣀⠀⠰⡘⠀⠐⠖⠂⠀⠁⠄⠂⣀⡠⠻⢆⠄⠃<br>⠐⠁⠤⣁⠀⠁⠈⠀⠂⠐⠀⠰⠀⠀⠀⠀⠂⠁⠈⠀⠃⠌⠧⠄⠀⢀⠤⠁⠱⢆</pre> | <pre>⣕⢝⠀⠀⢸⠔⡵⢊⡀⠂⠀⠀⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣑⠑<br>⠀⠀⠑⢄⠀⠳⠡⢡⣒⣃⢣⠯⠆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠌<br>⢒⠖⢤⡀⠑⢄⢶⡈⣂⠎⢎⠉⠩⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀<br>⡱⢋⠅⣂⡘⠳⠻⢆⡥⣈⠆⡨⡩⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠀<br>⠠⠈⠼⢸⡨⠜⡁⢫⣻⢞⢔⠀⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠠⠠<br>⠀⠀⡭⡖⡎⠑⡈⡡⠐⠑⠵⣧⣜⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀<br>⠀⠁⠈⠁⠃⠂⠃⠊⠀⠘⠒⠙⠛⢄⠀⠀⢄⠀⠤⢠⠀⢄⢀⢀⠀⡀⠀⠀⢄⢄<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⠊⠀⣂⠅⢓⣤⡄⠢⠠⠀⠌⠉⢀⢁<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⠊⠀⠑⢄⠁⣋⠀⢀⢰⢄⢔⢠⡖⢥⠀⠁<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣃⠌⠜⡥⢠⠛⣤⠐⣂⡀⠀⡀⡁⠍⠤⠒⠀<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢄⠙⣴⠀⢀⠰⢠⠿⣧⡅⠁⠂⢂⠂⠋⢃⢀<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢐⠠⡉⠐⢖⠀⠈⠅⠉⢕⢕⠝⠘⡒⠠⠀⠀<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠠⠀⠂⠐⣑⠄⠨⠨⢀⣓⠁⣕⢝⡥⢉⠁⠠<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡆⠁⠜⣍⠃⡅⡬⠀⠘⡈⡅⢋⠛⣤⡅⠒<br>⢕⠘⡂⠄⠀⠀⠁⠀⠀⡂⠀⢠⠀⢕⠄⢐⠄⠀⠘⠀⠉⢐⠀⠀⠁⡀⢡⠉⢟⣵</pre> |
|:---------------------- |:--------------------------------- |
| `S` (5% stored values) | `S[perm,perm]` (5% stored values) |

We can also visualize the sparsity pattern of the Cholesky factorization of
the same matrix. It is here clear that using the fill reducing permutation
results in a sparser factorization:

|<pre>⠙⢤⢠⡄⠀⣜⠃⢠⠀⠐⠘⠀⠀⠀⠀⠛⠃⠀⠀⠀⠀⠀⠀⠀⠀⠘⠀⠂⡔⠀<br>⠀⠀⠙⢦⡇⠾⡃⠰⠀⠀⠀⠐⠀⠃⠀⢂⠀⠀⠠⠀⠀⠀⠃⢀⠀⢀⠀⠀⠆⢣<br>⠀⠀⠀⠀⠙⢼⣣⢠⠀⣂⣂⢘⡂⡃⢰⣋⡀⣔⢠⠀⠀⠀⡃⠈⠀⢈⠀⡄⣄⡋<br>⠀⠀⠀⠀⠀⠀⠑⢖⡰⠉⠉⠈⠁⠁⢘⢙⠉⡊⢐⢐⢀⣞⠱⠎⠀⠌⡀⡣⡊⠉<br>⠀⠀⠀⠀⠀⠀⠀⠀⠙⢤⣴⢸⣴⡖⢠⣤⡜⢣⠀⠀⠛⠛⡜⠂⠀⢢⠀⡔⢸⡄<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢼⣛⣛⣛⣛⣓⣚⡃⢠⣖⣒⣓⢐⢠⣜⠀⡃⢘⣓<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿⣿⣿⣿⣿⣿⢸⣿⣿⣿⣾⣿⣿⠀⣿⢸⣿<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣒⣿⣺⣿<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿⣿⣿⣿⣿⣿⣿⣿⣿⣤⣿⣼⣿<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿⣿⣿⣿⣿⣿⣿⣿⣿<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿⣿⣿⣿⣿⣿⣿<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿⣿⣿⣿⣿<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿⣿⣿<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿</pre> | <pre>⠑⢝⠀⠀⢸⠔⡵⢊⡀⡂⠀⠀⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣕⢕<br>⠀⠀⠑⢄⠀⠳⠡⢡⣒⣃⢣⠯⠆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠌<br>⠀⠀⠀⠀⠑⢄⢶⡘⣂⡎⢎⡭⠯⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠶⠴<br>⠀⠀⠀⠀⠀⠀⠙⢎⣷⣏⢷⣯⡫⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠛⡛<br>⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿⢼⣧⣧⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠤⡤<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢿⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣭⣯<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢄⠀⠀⢄⠀⠤⢠⠀⢄⢀⢀⠀⡀⠀⠀⢟⢟<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⠊⠀⣂⠅⢓⣤⡄⠢⠠⠀⠌⠉⢀⢁<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⠉⣋⠀⢁⢰⢔⢔⢠⡖⢥⠁⠃<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢤⠘⣶⡂⠠⡀⣡⠭⣤⢓⢗<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢷⡇⡇⣢⣢⠂⣯⣷⣶<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢕⢟⢝⣒⠭⠭⡭<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢝⣿⣿⡭⡯<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿⣿⣿<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿</pre> |
|:----------------------------- |:--------------------------------------- |
| `chol(S)` (16% stored values) | `chol(S[perm,perm])` (6% stored values) |

## Direct access to the Metis C API
For more fine tuned usage of Metis consider calling the C API directly.
The wrappers aren't exported, so call them with the module prefix, e.g. `Metis.METIS_NodeND`.
The following functions are currently exposed:
- `METIS_PartGraphRecursive`
- `METIS_PartGraphKway`
- `METIS_PartMeshNodal`
- `METIS_PartMeshDual`
- `METIS_MeshToNodal`
- `METIS_MeshToDual`
- `METIS_NodeND`
- `METIS_SetDefaultOptions`
- `METIS_Free`

all with the same arguments and argument order as described in the
[Metis manual][metis-manual-url].

The following functions are also exposed. The manual doesn't describe them, so see
`metis.h` (included in METIS_jll) for their arguments, which the wrappers follow:
- `METIS_ComputeVertexSeparator`
- `METIS_NodeNDP`
- `METIS_NodeRefine`

[gh-actions-img]: https://github.com/JuliaSparse/Metis.jl/workflows/CI/badge.svg
[gh-actions-url]: https://github.com/JuliaSparse/Metis.jl/actions?query=workflow%3ACI

[codecov-img]: http://codecov.io/github/JuliaSparse/Metis.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/JuliaSparse/Metis.jl?branch=master

[metis-url]: https://karypis.github.io/glaros/software/metis/overview.html
[metis-manual-url]: https://karypis.github.io/glaros/files/sw/metis/manual.pdf

[partition2-url]: https://user-images.githubusercontent.com/11698744/38196819-65950f1e-3685-11e8-8db4-6aa9563bbd62.png
[partition3-url]: https://user-images.githubusercontent.com/11698744/38196820-65b11c9a-3685-11e8-95a0-b3b280359b31.png
[partition4-url]: https://user-images.githubusercontent.com/11698744/38196821-65ddc1dc-3685-11e8-8eb1-ce44ef1646f3.png
[separator-url]: https://user-images.githubusercontent.com/11698744/38196822-65fffc34-3685-11e8-9575-4dba41faec41.png
