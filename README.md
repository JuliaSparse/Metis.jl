# Metis.jl

[![][travis-img]][travis-url] [![][codecov-img]][codecov-url]

*Metis.jl* is a Julia wrapper to the [Metis library][metis-url] which is a
library for partitioning unstructured graphs, partitioning meshes, and
computing fill-reducing orderings of sparse matrices.

## Fill reducing permutation
`Metis.permutation` can be used to calculate the fill reducing permutation
for sparse matrices. `Metis.permutation` calls `METIS_NodeND` from the Metis
C API. As an example, we calculate the fill reducing permutation
for a sparse matrix `S` originating from a typical FEM problem, and visualize
the sparsity pattern for the original matrix and the permuted matrix.

```julia
perm, iperm = Metis.permutation(S)
```

| `spy(S)` (5% stored values) | `spy(S[perm,perm])` (5% stored values) |
|:--------------------------- |:-------------------------------------- |
| <pre>⠛⣤⢠⠄⠀⣌⠃⢠⠀⠐⠈⠀⠀⠀⠀⠉⠃⠀⠀⠀⠀⠀⠀⠀⠀⠘⠀⠂⠔⠀<br>⠀⠖⠻⣦⡅⠘⡁⠀⠀⠀⠀⠐⠀⠁⠀⢂⠀⠀⠠⠀⠀⠀⠁⢀⠀⢀⠀⠀⠄⢣<br>⡀⢤⣁⠉⠛⣤⡡⢀⠀⠂⠂⠀⠂⠃⢰⣀⠀⠔⠀⠀⠀⠀⠀⠀⠀⠀⠀⠄⠄⠀<br>⠉⣀⠁⠈⠁⢊⠱⢆⡰⠀⠈⠀⠀⠀⠀⢈⠉⡂⠀⠐⢀⡞⠐⠂⠀⠄⡀⠠⠂⠀<br>⢀⠀⠀⠀⠠⠀⠐⠊⠛⣤⡔⠘⠰⠒⠠⠀⡈⠀⠀⠀⠉⠉⠘⠂⠀⠀⠀⡐⢈⠀<br>⠂⠀⢀⠀⠈⠀⠂⠀⣐⠉⢑⣴⡉⡈⠁⡂⠒⠀⠁⢠⡄⠀⠐⠀⠠⠄⠀⠁⢀⡀<br>⠀⠀⠄⠀⠬⠀⠀⠀⢰⠂⡃⠨⣿⣿⡕⠂⠀⠨⠌⠈⠆⠀⠄⡀⠑⠀⠀⠘⠀⠀<br>⡄⠀⠠⢀⠐⢲⡀⢀⠀⠂⠡⠠⠱⠉⢱⢖⡀⠀⡈⠃⠀⠀⠀⢁⠄⢀⣐⠢⠀⠀<br>⠉⠀⠀⠀⢀⠄⠣⠠⠂⠈⠘⠀⡀⡀⠀⠈⠱⢆⣰⠠⠰⠐⠐⢀⠀⢀⢀⠀⠌⠀<br>⠀⠀⠀⠂⠀⠀⢀⠀⠀⠀⠁⣀⡂⠁⠦⠈⠐⡚⠱⢆⢀⢀⠡⠌⡀⡈⠸⠁⠂⠀<br>⠀⠀⠀⠀⠀⠀⣠⠴⡇⠀⠀⠉⠈⠁⠀⠀⢐⠂⠀⢐⣻⣾⠡⠀⠈⠀⠄⠀⡉⠄<br>⠀⠀⠁⢀⠀⠀⠰⠀⠲⠀⠐⠀⠀⠡⠄⢀⠐⢀⡁⠆⠁⠂⠱⢆⡀⣀⠠⠁⠉⠇<br>⣀⠀⠀⢀⠀⠀⠀⠄⠀⠀⠀⠆⠑⠀⠀⢁⠀⢀⡀⠨⠂⠀⠀⢨⠿⢇⠀⡸⠀⢀<br>⠠⠀⠀⠀⠀⠄⠀⡈⢀⠠⠄⠀⣀⠀⠰⡘⠀⠐⠖⠂⠀⠁⠄⠂⣀⡠⠻⢆⠄⠃<br>⠐⠁⠤⣁⠀⠁⠈⠀⠂⠐⠀⠰⠀⠀⠀⠀⠂⠁⠈⠀⠃⠌⠧⠄⠀⢀⠤⠁⠱⢆</pre> | <pre>⣕⢝⠀⠀⢸⠔⡵⢊⡀⠂⠀⠀⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣑⠑<br>⠀⠀⠑⢄⠀⠳⠡⢡⣒⣃⢣⠯⠆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠌<br>⢒⠖⢤⡀⠑⢄⢶⡈⣂⠎⢎⠉⠩⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀<br>⡱⢋⠅⣂⡘⠳⠻⢆⡥⣈⠆⡨⡩⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠀<br>⠠⠈⠼⢸⡨⠜⡁⢫⣻⢞⢔⠀⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠠⠠<br>⠀⠀⡭⡖⡎⠑⡈⡡⠐⠑⠵⣧⣜⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣀<br>⠀⠁⠈⠁⠃⠂⠃⠊⠀⠘⠒⠙⠛⢄⠀⠀⢄⠀⠤⢠⠀⢄⢀⢀⠀⡀⠀⠀⢄⢄<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⠊⠀⣂⠅⢓⣤⡄⠢⠠⠀⠌⠉⢀⢁<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⠊⠀⠑⢄⠁⣋⠀⢀⢰⢄⢔⢠⡖⢥⠀⠁<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣃⠌⠜⡥⢠⠛⣤⠐⣂⡀⠀⡀⡁⠍⠤⠒⠀<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢄⠙⣴⠀⢀⠰⢠⠿⣧⡅⠁⠂⢂⠂⠋⢃⢀<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢐⠠⡉⠐⢖⠀⠈⠅⠉⢕⢕⠝⠘⡒⠠⠀⠀<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠠⠀⠂⠐⣑⠄⠨⠨⢀⣓⠁⣕⢝⡥⢉⠁⠠<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡆⠁⠜⣍⠃⡅⡬⠀⠘⡈⡅⢋⠛⣤⡅⠒<br>⢕⠘⡂⠄⠀⠀⠁⠀⠀⡂⠀⢠⠀⢕⠄⢐⠄⠀⠘⠀⠉⢐⠀⠀⠁⡀⢡⠉⢟⣵</pre> |

We may also visualize the sparsity pattern of the Cholesky factorization of
the same matrix. It is here clear that using the fill reducing permutation
results in a sparser factorization.

| `spy(chol(S))` (16% stored values) | `spy(chol(S[perm,perm]))` (6% stored values) |
|:---------------------------------- |:-------------------------------------------- |
|<pre>⠙⢤⢠⡄⠀⣜⠃⢠⠀⠐⠘⠀⠀⠀⠀⠛⠃⠀⠀⠀⠀⠀⠀⠀⠀⠘⠀⠂⡔⠀<br>⠀⠀⠙⢦⡇⠾⡃⠰⠀⠀⠀⠐⠀⠃⠀⢂⠀⠀⠠⠀⠀⠀⠃⢀⠀⢀⠀⠀⠆⢣<br>⠀⠀⠀⠀⠙⢼⣣⢠⠀⣂⣂⢘⡂⡃⢰⣋⡀⣔⢠⠀⠀⠀⡃⠈⠀⢈⠀⡄⣄⡋<br>⠀⠀⠀⠀⠀⠀⠑⢖⡰⠉⠉⠈⠁⠁⢘⢙⠉⡊⢐⢐⢀⣞⠱⠎⠀⠌⡀⡣⡊⠉<br>⠀⠀⠀⠀⠀⠀⠀⠀⠙⢤⣴⢸⣴⡖⢠⣤⡜⢣⠀⠀⠛⠛⡜⠂⠀⢢⠀⡔⢸⡄<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢼⣛⣛⣛⣛⣓⣚⡃⢠⣖⣒⣓⢐⢠⣜⠀⡃⢘⣓<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿⣿⣿⣿⣿⣿⢸⣿⣿⣿⣾⣿⣿⠀⣿⢸⣿<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣒⣿⣺⣿<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿⣿⣿⣿⣿⣿⣿⣿⣿⣤⣿⣼⣿<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿⣿⣿⣿⣿⣿⣿⣿⣿<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿⣿⣿⣿⣿⣿⣿<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿⣿⣿⣿⣿<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿⣿⣿<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿</pre> | <pre>⠑⢝⠀⠀⢸⠔⡵⢊⡀⡂⠀⠀⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣕⢕<br>⠀⠀⠑⢄⠀⠳⠡⢡⣒⣃⢣⠯⠆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠌<br>⠀⠀⠀⠀⠑⢄⢶⡘⣂⡎⢎⡭⠯⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠶⠴<br>⠀⠀⠀⠀⠀⠀⠙⢎⣷⣏⢷⣯⡫⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠛⡛<br>⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿⢼⣧⣧⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠤⡤<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢿⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣭⣯<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢄⠀⠀⢄⠀⠤⢠⠀⢄⢀⢀⠀⡀⠀⠀⢟⢟<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⠊⠀⣂⠅⢓⣤⡄⠢⠠⠀⠌⠉⢀⢁<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢄⠉⣋⠀⢁⢰⢔⢔⢠⡖⢥⠁⠃<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢤⠘⣶⡂⠠⡀⣡⠭⣤⢓⢗<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢷⡇⡇⣢⣢⠂⣯⣷⣶<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢕⢟⢝⣒⠭⠭⡭<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠑⢝⣿⣿⡭⡯<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿⣿⣿<br>⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢿</pre> |


[travis-img]: https://travis-ci.org/fredrikekre/Metis.jl.svg?branch=master
[travis-url]: https://travis-ci.org/fredrikekre/Metis.jl

[codecov-img]: http://codecov.io/github/fredrikekre/Metis.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/fredrikekre/Metis.jl?branch=master

[metis-url]: http://glaros.dtc.umn.edu/gkhome/metis/metis/overview