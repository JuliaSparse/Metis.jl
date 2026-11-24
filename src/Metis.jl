# SPDX-License-Identifier: MIT

module Metis

using SparseArrays
using LinearAlgebra: Hermitian, Symmetric
using METIS_jll: libmetis

# Metis C API: Clang.jl auto-generated bindings and some manual methods
include("LibMetis.jl")
using .LibMetis
using .LibMetis: idx_t, @check

# Global options array (deprecated). It is only kept since existing code refers to it: other
# packages pass it to the C API with 1-based arrays, so it keeps 1-based numbering, and change
# options in it before calling `partition` etc., so the Julia interface still uses it when no
# `Options` is passed (see `metis_options`). Remove it in the next breaking release.
const options = fill(Cint(-1), METIS_NOPTIONS)
options[Int(METIS_OPTION_NUMBERING)+1] = 1

# Metis has global state: every call reseeds the C library's random number generator, and
# most calls install and then restore process-wide signal handlers. Calls running at the
# same time can therefore give non-reproducible results or leave the wrong signal handlers
# installed, so they are serialized with this lock.
const metis_lock = ReentrantLock()

"""
    Metis.Options(; kwargs...)

Options for Metis, passed to `Metis.partition`, `Metis.permutation` and `Metis.separator` with
the `options` keyword argument, e.g.

```julia
Metis.partition(G, 4; options = Metis.Options(seed = 42, ufactor = 30))
```

Each keyword argument sets one Metis option, see section 5.4 "Options array" in the Metis
manual for their meaning and values. Values can be integers, `Bool`s or the corresponding
constants, e.g. `objtype = Metis.METIS_OBJTYPE_VOL`. Options that aren't given use the Metis
default. Metis throws a `MetisError` for an invalid value. The options, and the functions that
use them, are:

| Keyword     | Metis option             | Used by                                                         |
|:----------- |:------------------------ |:--------------------------------------------------------------- |
| `ptype`     | `METIS_OPTION_PTYPE`     | Only the mesh partitioning functions of the C API               |
| `objtype`   | `METIS_OPTION_OBJTYPE`   | `partition`; must not be set for `permutation` and `separator`  |
| `ctype`     | `METIS_OPTION_CTYPE`     | All                                                             |
| `iptype`    | `METIS_OPTION_IPTYPE`    | `partition` with `RECURSIVE`, `permutation`, `separator`        |
| `rtype`     | `METIS_OPTION_RTYPE`     | `permutation`, `separator`                                      |
| `dbglvl`    | `METIS_OPTION_DBGLVL`    | All                                                             |
| `niter`     | `METIS_OPTION_NITER`     | All                                                             |
| `ncuts`     | `METIS_OPTION_NCUTS`     | `partition`                                                     |
| `seed`      | `METIS_OPTION_SEED`      | All                                                             |
| `no2hop`    | `METIS_OPTION_NO2HOP`    | All                                                             |
| `minconn`   | `METIS_OPTION_MINCONN`   | `partition` with `KWAY`                                         |
| `contig`    | `METIS_OPTION_CONTIG`    | `partition` with `KWAY`                                         |
| `compress`  | `METIS_OPTION_COMPRESS`  | `permutation`, `separator`                                      |
| `ccorder`   | `METIS_OPTION_CCORDER`   | `permutation`                                                   |
| `pfactor`   | `METIS_OPTION_PFACTOR`   | `permutation`                                                   |
| `nseps`     | `METIS_OPTION_NSEPS`     | `permutation`, `separator`                                      |
| `ufactor`   | `METIS_OPTION_UFACTOR`   | All                                                             |
| `numbering` | `METIS_OPTION_NUMBERING` | Only the C API, see below                                       |

An `Options` stores the options array in the layout Metis uses, and can be passed directly to
the C API functions, e.g. `Metis.METIS_NodeND`, with `numbering = 1` for 1-based arrays.
`Metis.partition`, `Metis.permutation` and `Metis.separator` pass 0-based arrays to Metis, and
throw an `ArgumentError` for `numbering = 1`.

When `options` is passed, the global `Metis.options` array isn't used. Without it, options set
in the global array are still used, but setting them there is deprecated in favour of
`Metis.Options`.
"""
struct Options
    values::Vector{Cint}
    function Options(; kwargs...)
        values = fill(Cint(-1), METIS_NOPTIONS)
        for (name, value) in kwargs
            values[option_index(name)] = value
        end
        return new(values)
    end
end

# The options that can be set with `Options`: all options that the Metis library reads. Not
# included are the options after `numbering` (`help`, `tpwgts`, `ncommon`, `nooutput`,
# `balance`, `gtype` and `ubvec`), which are only used by the Metis command-line programs.
const OPTION_NAMES = (
    :ptype, :objtype, :ctype, :iptype, :rtype, :dbglvl, :niter, :ncuts, :seed,
    :no2hop, :minconn, :contig, :compress, :ccorder, :pfactor, :nseps, :ufactor, :numbering,
)

# The index in the options array of each option in `OPTION_NAMES`
const OPTION_INDICES = map(OPTION_NAMES) do name
    return Int(getfield(LibMetis, Symbol("METIS_OPTION_", uppercase(String(name))))) + 1
end

function option_index(name::Symbol)
    i = findfirst(==(name), OPTION_NAMES)
    if i === nothing
        throw(ArgumentError("unknown Metis option `$(name)`"))
    end
    return OPTION_INDICES[i]
end

# Pass the options array when an `Options` is passed to the C API
Base.cconvert(::Type{Ptr{idx_t}}, opts::Options) = opts.values

# Setting options in the global `options` array, instead of with `Options`, is deprecated.
# The array is still used, but changing an option from its default gives a deprecation
# warning. Numbering is 1-based by default since the array is meant for the C API.
function check_global_options()
    numbering = option_index(:numbering)
    for (i, value) in pairs(options)
        if value != (i == numbering ? Cint(1) : Cint(-1))
            Base.depwarn(
                "setting options in the global `Metis.options` array is deprecated, use the " *
                "`options` keyword argument with `Metis.Options` instead, e.g. " *
                "`Metis.partition(G, n; options = Metis.Options(seed = 42))`",
                :check_global_options,
            )
            return
        end
    end
    return
end

# Return the options for a call to Metis with the arrays from `zero_based_adjacency`, i.e. with
# 0-based numbering. Without the `options` keyword argument (`nothing`) this is a copy of the
# global options, so that code that sets options there keeps working (deprecated, see
# `check_global_options`). This is the only use of the global array in the Julia interface, and
# it should be removed in the next breaking release. Call it while holding `metis_lock`.
function metis_options(::Nothing)
    o = copy(options)
    o[option_index(:numbering)] = 0
    return o
end

# Metis only reads the options array in the functions used by the Julia interface, so `opts`
# is passed as is. The Julia interface passes 0-based arrays, so 1-based numbering (which is
# for using `Options` with the C API) is an error.
function metis_options(opts::Options)
    if !(opts.values[option_index(:numbering)] in (-1, 0))
        throw(ArgumentError(
            "`numbering = 1` is only for the C API, `Metis.partition`, `Metis.permutation` " *
            "and `Metis.separator` pass 0-based arrays to Metis"
        ))
    end
    return opts
end

# Julia interface
"""
    Metis.Graph

1-based CSR representation of a graph as defined in
section 5.5 "Graph data structure" in the Metis manual.

Metis requires a graph with symmetric structure and no self edges or duplicate edges: every
edge `i → j` must be stored once, with a matching edge `j → i` of the same weight. Metis doesn't
validate its input and can crash, hang or give invalid results otherwise. `Metis.graph` checks
this when converting other types, but a `Graph` constructed directly is passed to Metis as is.
"""
struct Graph
    nvtxs::idx_t
    xadj::Vector{idx_t}
    adjncy::Union{Vector{idx_t}, Ptr{Nothing}}
    vwgt::Union{Vector{idx_t}, Ptr{Nothing}}
    adjwgt::Union{Vector{idx_t}, Ptr{Nothing}}
    function Graph(nvtxs, xadj, adjncy, vwgt=C_NULL, adjwgt=C_NULL)
        return new(nvtxs, xadj, adjncy, vwgt, adjwgt)
    end
end

# Convert the weight `w` of the edge `(i, j)` to the integer type used by Metis.
# Checks `isreal` since a complex Hermitian matrix can store real integers as `Complex`.
function edge_weight(w, i, j)
    if !(isreal(w) && isinteger(real(w)) && real(w) > 0)
        throw(ArgumentError("weights must be positive integers, got weight $(w) for edge ($(i), $(j))"))
    end
    return idx_t(real(w))
end

# Find a stored entry (i, j) off the diagonal of the square matrix G without a stored entry
# (j, i), or, if `values` is true, with a different value than the entry (j, i). Return
# `(i, j, :pattern)` or `(i, j, :value)` for the first one found, and `nothing` if there is
# none. Unlike `issymmetric` this considers stored zeros, and unlike comparing with the
# transpose it doesn't allocate.
function find_asymmetry(G::SparseMatrixCSC, values::Bool)
    rows = rowvals(G)
    vals = nonzeros(G)
    for j in axes(G, 2), k in nzrange(G, j)
        i = rows[k]
        i == j && continue
        r = nzrange(G, i)
        l = searchsortedfirst(view(rows, r), j)
        (l <= length(r) && rows[r[l]] == j) || return (i, j, :pattern)
        values && vals[r[l]] != vals[k] && return (i, j, :value)
    end
    return nothing
end

"""
    Metis.graph(G::SparseMatrixCSC; weights=false, check=true)

Construct the 1-based CSR representation of the square sparse matrix `G`.
Every stored entry off the diagonal becomes an edge, including explicitly stored zeros.

Metis requires the stored entries to have a symmetric pattern, and with `weights=true` also
symmetric values, i.e. the same weight in both directions of each edge. If `check` is `true`
this is checked and an `ArgumentError` is thrown otherwise. A stored zero at `(i, j)` without a
stored entry at `(j, i)` breaks the pattern, and can be removed with `SparseArrays.dropzeros`.
If `check` is `false` nothing is checked, and Metis can crash or give invalid results for a
matrix that doesn't satisfy this.

If `weights=true` the entries of the matrix are used as edge weights. These must be positive
integers, otherwise an `ArgumentError` is thrown. This includes explicitly stored zeros, which
can be removed with `SparseArrays.dropzeros`. Diagonal entries are ignored.
"""
function graph(G::SparseMatrixCSC; weights::Bool=false, check::Bool=true,
               check_hermitian::Union{Bool, Nothing}=nothing)
    if check_hermitian !== nothing
        Base.depwarn("the `check_hermitian` keyword argument is deprecated, use `check` instead", :graph)
        check = check_hermitian
    end
    size(G, 1) == size(G, 2) || throw(ArgumentError("matrix must be square, got size $(size(G))"))
    if check
        asymmetry = find_asymmetry(G, weights)
        if asymmetry !== nothing
            i, j, kind = asymmetry
            if kind === :pattern
                throw(ArgumentError(
                    "the stored entries of the matrix must have a symmetric pattern, since each " *
                    "of them becomes an edge, but entry ($(i), $(j)) is stored and entry " *
                    "($(j), $(i)) isn't. One-sided stored zeros can be removed with " *
                    "`SparseArrays.dropzeros`."
                ))
            else
                throw(ArgumentError(
                    "with `weights=true` the matrix must be symmetric, since its values are the " *
                    "edge weights, but entry ($(i), $(j)) is $(G[i, j]) and entry ($(j), $(i)) " *
                    "is $(G[j, i])"
                ))
            end
        end
    end
    N = size(G, 1)
    xadj = Vector{idx_t}(undef, N+1)
    xadj[1] = 1
    adjncy = Vector{idx_t}(undef, nnz(G))
    vwgt = C_NULL # TODO: Vertex weights could be passed as input argument
    adjwgt = weights ? Vector{idx_t}(undef, nnz(G)) : C_NULL
    adjncy_i = 0
    @inbounds for j in 1:N
        n_rows = 0
        for k in G.colptr[j] : (G.colptr[j+1] - 1)
            i = G.rowval[k]
            i == j && continue # skip self edges
            n_rows += 1
            adjncy_i += 1
            adjncy[adjncy_i] = i
            if weights
                adjwgt[adjncy_i] = edge_weight(G.nzval[k], i, j)
            end
        end
        xadj[j+1] = xadj[j] + n_rows
    end
    resize!(adjncy, adjncy_i)
    weights && resize!(adjwgt, adjncy_i)
    return Graph(idx_t(N), xadj, adjncy, vwgt, adjwgt)
end

const HermOrSymCSC{Tv, Ti} = Union{
    Hermitian{Tv, SparseMatrixCSC{Tv, Ti}}, Symmetric{Tv, SparseMatrixCSC{Tv, Ti}},
}

using SparseArrays: nzrangeup, nzrangelo

"""
    Metis.graph(G::Union{Hermitian, Symmetric}; weights::Bool=false)

Construct the 1-based CSR representation of the `Hermitian` or `Symmetric` wrapped sparse
matrix `G`.
Weights are not currently supported for this method so passing `weights=true` will throw an
error.
"""
function graph(H::HermOrSymCSC; weights::Bool=false)
    # This method is derived from the method `SparseMatrixCSC(::HermOrSymCSC)` from
    # SparseArrays.jl
    # (https://github.com/JuliaSparse/SparseArrays.jl/blob/313a04f4a78bbc534f89b6b4d9c598453e2af17c/src/sparseconvert.jl#L124-L173)
    # with MIT license
    # (https://github.com/JuliaSparse/SparseArrays.jl/blob/main/LICENSE.md).
    weights && throw(ArgumentError("weights not supported yet"))
    # Extract data
    A = H.data
    upper = H.uplo == 'U'
    rowval = rowvals(A)
    m, n = size(A)
    @assert m == n
    # New colptr for the full matrix
    newcolptr = Vector{idx_t}(undef, n + 1)
    newcolptr[1] = 1
    # SparseArrays.nzrange for the upper/lower part excluding the diagonal
    nzrng = if upper
        (A, col) -> nzrangeup(A, col, #=exclude diagonal=# true)
    else
        (A, col) -> nzrangelo(A, col, #=exclude diagonal=# true)
    end
    # If the upper part is stored we loop forward, otherwise backwards
    colrange = upper ? (1:1:n) : (n:-1:1)
    @inbounds for col in colrange
        r = nzrng(A, col)
        # Number of entries in the stored part of this column, excluding the diagonal entry
        newcolptr[col + 1] = length(r)
        # Increment columnptr corresponding to the stored rows
        for k in r
            row = rowval[k]
            @assert upper ? row < col : row > col
            @assert row != col # Diagonal entries should not be here
            newcolptr[row + 1] += 1
        end
    end
    # Accumulate the colptr and allocate new rowval
    cumsum!(newcolptr, newcolptr)
    nz = newcolptr[n + 1] - 1
    newrowval = Vector{idx_t}(undef, nz)
    # Populate the rowvals
    @inbounds for col = 1:n
        newk = newcolptr[col]
        for k in nzrng(A, col)
            row = rowval[k]
            @assert col != row
            newrowval[newk] = row
            newk += 1
            ni = newcolptr[row]
            newrowval[ni] = col
            newcolptr[row] = ni + 1
        end
        newcolptr[col] = newk
    end
    # Shuffle back the colptrs
    @inbounds for j = n:-1:1
        newcolptr[j+1] = newcolptr[j]
    end
    newcolptr[1] = 1
    # Return Graph
    N = n
    xadj = newcolptr
    adjncy = newrowval
    vwgt = C_NULL
    adjwgt = C_NULL
    return Graph(idx_t(N), xadj, adjncy, vwgt, adjwgt)
end

# Return 0-based copies of the adjacency arrays of `G` to pass to Metis. With 1-based
# numbering Metis converts the caller's arrays to 0-based and back in place, so they would be
# invalid for anyone else reading them during the call. METIS_ComputeVertexSeparator also
# ignores the numbering option and always expects 0-based arrays. `adjncy` can be `C_NULL`
# for a graph without edges.
function zero_based_adjacency(G::Graph)
    xadj = G.xadj .- idx_t(1)
    adjncy = G.adjncy isa Vector ? G.adjncy .- idx_t(1) : G.adjncy
    return xadj, adjncy
end

"""
    perm, iperm = Metis.permutation(G; options = nothing)

Compute the fill reducing permutation `perm`
and its inverse `iperm` of `G`.
See `Metis.Options` for the `options` keyword argument.
"""
permutation(G; options::Union{Options, Nothing} = nothing) = permutation(graph(G); options)

function permutation(G::Graph; options::Union{Options, Nothing} = nothing)
    perm = Vector{idx_t}(undef, G.nvtxs)
    iperm = Vector{idx_t}(undef, G.nvtxs)
    options === nothing && check_global_options()
    xadj, adjncy = zero_based_adjacency(G)
    @lock metis_lock begin
        @check METIS_NodeND(Ref{idx_t}(G.nvtxs), xadj, adjncy, G.vwgt, metis_options(options), perm, iperm)
    end
    perm .+= idx_t(1)
    iperm .+= idx_t(1)
    return perm, iperm
end

"""
    Metis.Algorithm

The partitioning algorithm for `Metis.partition`, passed with the `alg` keyword argument:
 - `Metis.Algorithm.KWAY`: multilevel k-way partitioning (`METIS_PartGraphKway`)
 - `Metis.Algorithm.RECURSIVE`: multilevel recursive bisection (`METIS_PartGraphRecursive`)
"""
baremodule Algorithm
    using Base: @enum
    @enum T KWAY RECURSIVE
end

# Convert the `alg` keyword argument of `partition`, which can also be a `Symbol` (deprecated)
to_algorithm(alg::Algorithm.T) = alg
function to_algorithm(alg)
    if alg === :KWAY
        enum = Algorithm.KWAY
    elseif alg === :RECURSIVE
        enum = Algorithm.RECURSIVE
    else
        throw(ArgumentError("unknown algorithm $(repr(alg))"))
    end
    Base.depwarn(
        "passing `alg` as a `Symbol` is deprecated, use `alg = Metis.Algorithm.$(alg)` instead",
        :to_algorithm,
    )
    return enum
end

"""
    Metis.partition(G, n; alg = Metis.Algorithm.KWAY, options = nothing)

Partition the graph `G` in `n` parts, with the algorithm given by `alg`, see
`Metis.Algorithm`.

An `ArgumentError` is thrown for any other `alg`, if `n` is less than 1, and if `n` is larger
than the number of vertices.
See `Metis.Options` for the `options` keyword argument.
"""
partition(G, nparts; alg = Algorithm.KWAY, options::Union{Options, Nothing} = nothing) =
    partition(graph(G), nparts; alg, options)

function partition(G::Graph, nparts::Integer; alg = Algorithm.KWAY, options::Union{Options, Nothing} = nothing)
    alg = to_algorithm(alg)
    nparts < 1 && throw(ArgumentError("the number of parts must be at least 1, got $(nparts)"))
    if nparts > 1 && nparts > G.nvtxs
        throw(ArgumentError("the number of parts ($(nparts)) must not be larger than the number of vertices ($(G.nvtxs))"))
    end
    part = Vector{idx_t}(undef, G.nvtxs)
    nparts == 1 && return fill!(part, 1) # https://github.com/JuliaSparse/Metis.jl/issues/49
    edgecut = fill(idx_t(0), 1)
    options === nothing && check_global_options()
    xadj, adjncy = zero_based_adjacency(G)
    @lock metis_lock begin
        if alg == Algorithm.RECURSIVE
            @check METIS_PartGraphRecursive(Ref{idx_t}(G.nvtxs), Ref{idx_t}(1), xadj, adjncy, G.vwgt, C_NULL, G.adjwgt,
                                            Ref{idx_t}(nparts), C_NULL, C_NULL, metis_options(options), edgecut, part)
        else
            @check METIS_PartGraphKway(Ref{idx_t}(G.nvtxs), Ref{idx_t}(1), xadj, adjncy, G.vwgt, C_NULL, G.adjwgt,
                                       Ref{idx_t}(nparts), C_NULL, C_NULL, metis_options(options), edgecut, part)
        end
    end
    part .+= idx_t(1)
    return part
end

"""
    Metis.separator(G; options = nothing)

Compute a vertex separator of the graph `G`.
See `Metis.Options` for the `options` keyword argument.

Return a vector with one label per vertex: `1` and `2` for the two parts and `3` for the
separator. No edge joins part `1` and part `2`. A part or the separator can be empty, in
which case its label is missing from the result.
"""
separator(G; options::Union{Options, Nothing} = nothing) = separator(graph(G); options)

function separator(G::Graph; options::Union{Options, Nothing} = nothing)
    part = Vector{idx_t}(undef, G.nvtxs)
    sepsize = fill(idx_t(0), 1)
    options === nothing && check_global_options()
    xadj, adjncy = zero_based_adjacency(G)
    @lock metis_lock begin
        @check METIS_ComputeVertexSeparator(Ref{idx_t}(G.nvtxs), xadj, adjncy, G.vwgt, metis_options(options), sepsize, part)
    end
    part .+= idx_t(1)
    return part
end


end # module
