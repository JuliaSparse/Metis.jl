# SPDX-License-Identifier: MIT

module Metis

using SparseArrays
using LinearAlgebra: Hermitian, Symmetric
using METIS_jll: libmetis

# Metis C API: Clang.jl auto-generated bindings and some manual methods
include("LibMetis.jl")
using .LibMetis
using .LibMetis: idx_t, @check

# Global options array -- should probably do something better...
# Other packages pass this to the C API with 1-based arrays, and change other options in it
# before calling `partition` etc., so it keeps 1-based numbering. `metis_options` makes the
# copy that the Julia interface passes to Metis.
const options = fill(Cint(-1), METIS_NOPTIONS)
options[Int(METIS_OPTION_NUMBERING)+1] = 1

# Metis has global state: every call reseeds the C library's random number generator, and
# most calls install and then restore process-wide signal handlers. Calls running at the
# same time can therefore give non-reproducible results or leave the wrong signal handlers
# installed, so they are serialized with this lock.
const metis_lock = ReentrantLock()

# Return a copy of the global options with 0-based numbering, to use with the arrays from
# `zero_based_adjacency`. Call it while holding `metis_lock`.
function metis_options()
    opts = copy(options)
    opts[Int(METIS_OPTION_NUMBERING)+1] = 0
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

if VERSION < v"1.10"
    # From https://github.com/JuliaSparse/SparseArrays.jl/blob/313a04f4a78bbc534f89b6b4d9c598453e2af17c/src/linalg.jl#L1106-L1117
    # MIT license: https://github.com/JuliaSparse/SparseArrays.jl/blob/main/LICENSE.md
    function nzrangeup(A, i, excl=false)
        r = nzrange(A, i); r1 = r.start; r2 = r.stop
        rv = rowvals(A)
        @inbounds r2 < r1 || rv[r2] <= i - excl ? r : r1:(searchsortedlast(view(rv, r1:r2), i - excl) + r1-1)
    end
    function nzrangelo(A, i, excl=false)
        r = nzrange(A, i); r1 = r.start; r2 = r.stop
        rv = rowvals(A)
        @inbounds r2 < r1 || rv[r1] >= i + excl ? r : (searchsortedfirst(view(rv, r1:r2), i + excl) + r1-1):r2
    end
else
    using SparseArrays: nzrangeup, nzrangelo
end

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
    perm, iperm = Metis.permutation(G)

Compute the fill reducing permutation `perm`
and its inverse `iperm` of `G`.
"""
permutation(G) = permutation(graph(G))

function permutation(G::Graph)
    perm = Vector{idx_t}(undef, G.nvtxs)
    iperm = Vector{idx_t}(undef, G.nvtxs)
    xadj, adjncy = zero_based_adjacency(G)
    lock(metis_lock) do
        @check METIS_NodeND(Ref{idx_t}(G.nvtxs), xadj, adjncy, G.vwgt, metis_options(), perm, iperm)
    end
    perm .+= idx_t(1)
    iperm .+= idx_t(1)
    return perm, iperm
end

"""
    Metis.partition(G, n; alg = :KWAY)

Partition the graph `G` in `n` parts.
The partition algorithm is defined by the `alg` keyword:
 - :KWAY: multilevel k-way partitioning
 - :RECURSIVE: multilevel recursive bisection

An `ArgumentError` is thrown for any other `alg`, if `n` is less than 1, and if `n` is larger
than the number of vertices.
"""
partition(G, nparts; alg = :KWAY) = partition(graph(G), nparts, alg = alg)

function partition(G::Graph, nparts::Integer; alg = :KWAY)
    alg === :RECURSIVE || alg === :KWAY || throw(ArgumentError("unknown algorithm $(repr(alg))"))
    nparts < 1 && throw(ArgumentError("the number of parts must be at least 1, got $(nparts)"))
    if nparts > 1 && nparts > G.nvtxs
        throw(ArgumentError("the number of parts ($(nparts)) must not be larger than the number of vertices ($(G.nvtxs))"))
    end
    part = Vector{idx_t}(undef, G.nvtxs)
    nparts == 1 && return fill!(part, 1) # https://github.com/JuliaSparse/Metis.jl/issues/49
    edgecut = fill(idx_t(0), 1)
    xadj, adjncy = zero_based_adjacency(G)
    lock(metis_lock) do
        if alg === :RECURSIVE
            @check METIS_PartGraphRecursive(Ref{idx_t}(G.nvtxs), Ref{idx_t}(1), xadj, adjncy, G.vwgt, C_NULL, G.adjwgt,
                                            Ref{idx_t}(nparts), C_NULL, C_NULL, metis_options(), edgecut, part)
        else
            @check METIS_PartGraphKway(Ref{idx_t}(G.nvtxs), Ref{idx_t}(1), xadj, adjncy, G.vwgt, C_NULL, G.adjwgt,
                                       Ref{idx_t}(nparts), C_NULL, C_NULL, metis_options(), edgecut, part)
        end
    end
    part .+= idx_t(1)
    return part
end

"""
    Metis.separator(G)

Compute a vertex separator of the graph `G`.
"""
separator(G) = separator(graph(G))

function separator(G::Graph)
    part = Vector{idx_t}(undef, G.nvtxs)
    sepsize = fill(idx_t(0), 1)
    xadj, adjncy = zero_based_adjacency(G)
    lock(metis_lock) do
        @check METIS_ComputeVertexSeparator(Ref{idx_t}(G.nvtxs), xadj, adjncy, G.vwgt, metis_options(), sepsize, part)
    end
    part .+= idx_t(1)
    return part
end

# Compatibility for Julias that doesn't support package extensions
if !(isdefined(Base, :get_extension))
    include("../ext/MetisGraphs.jl")
    include("../ext/MetisLightGraphs.jl")
end

end # module
