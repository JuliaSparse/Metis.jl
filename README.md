# Metis

| **Build Status**                                                                                |
|:----------------------------------------------------------------------------------------------- |
| [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] [![][codecov-img]][codecov-url] |

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

## Vertex separator
`Metis.separator` calculates a [vertex separator](https://en.wikipedia.org/wiki/Vertex_separator)
of a graph. `Metis.separator` calls `METIS_ComputeVertexSeparator` from the Metis C API.
As an example, here we calculate a vertex separator (green) of a small graph:

| ![][separator-url]   |
|:-------------------- |
| `Metis.separator(g)` |

## Fill reducing permutation
`Metis.permutation` calculates the fill reducing permutation
for a sparse matrices. `Metis.permutation` calls `METIS_NodeND` from the Metis
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
The following functions are currently exposed:
- `METIS_PartGraphRecursive`
- `METIS_PartGraphKway`
- `METIS_ComputeVertexSeparator`
- `METIS_NodeND`

all with the same arguments and argument order as described in the
[Metis manual][metis-manual-url].


[travis-img]: https://travis-ci.org/JuliaSparse/Metis.jl.svg?branch=master
[travis-url]: https://travis-ci.org/JuliaSparse/Metis.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/76i2pnyq8495sgfa/branch/master?svg=true
[appveyor-url]: https://ci.appveyor.com/project/JuliaSparse/metis-jl/branch/master

[codecov-img]: http://codecov.io/github/JuliaSparse/Metis.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/JuliaSparse/Metis.jl?branch=master

[metis-url]: http://glaros.dtc.umn.edu/gkhome/metis/metis/overview
[metis-manual-url]: http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/manual.pdf

[S-url]: https://user-images.githubusercontent.com/11698744/38196722-dd9877c2-3684-11e8-8c02-a767604824d1.png
[Spp-url]: https://user-images.githubusercontent.com/11698744/38196723-ddb62fba-3684-11e8-89ff-181128644294.png
[C-url]: https://user-images.githubusercontent.com/11698744/38196720-dd5dd748-3684-11e8-8413-a52d336abe49.png
[Cpp-url]: https://user-images.githubusercontent.com/11698744/38196721-dd7ac16e-3684-11e8-8a35-761e97d11235.png
[partition2-url]: https://user-images.githubusercontent.com/11698744/38196819-65950f1e-3685-11e8-8db4-6aa9563bbd62.png
[partition3-url]: https://user-images.githubusercontent.com/11698744/38196820-65b11c9a-3685-11e8-95a0-b3b280359b31.png
[partition4-url]: https://user-images.githubusercontent.com/11698744/38196821-65ddc1dc-3685-11e8-8eb1-ce44ef1646f3.png
[separator-url]: https://user-images.githubusercontent.com/11698744/38196822-65fffc34-3685-11e8-9575-4dba41faec41.png
[vertex-separator-url]: https://en.wikipedia.org/wiki/Vertex_separator
