# SPDX-License-Identifier: MIT

module Metis

using SparseArrays
using LinearAlgebra: ishermitian
using METIS_jll: libmetis

# Metis C API: Clang.jl auto-generated bindings and some manual methods
include("LibMetis.jl")
using .LibMetis
using .LibMetis: idx_t, @check

# Global options array -- should probably do something better...
const options = fill(Cint(-1), METIS_NOPTIONS)
options[Int(METIS_OPTION_NUMBERING)+1] = 1

# Julia interface
"""
    Metis.Graph

1-based CSR representation of a graph as defined in
section 5.5 "Graph data structure" in the Metis manual.
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

"""
    Metis.graph(G::SparseMatrixCSC; weights=false, check_hermitian=true)

Construct the 1-based CSR representation of the sparse matrix `G`.
If `check_hermitian` is `false` the matrix is not checked for being hermitian
before constructing the graph.
If `weights=true` the entries of the matrix are used as edge weights.
"""
function graph(G::SparseMatrixCSC; weights::Bool=false, check_hermitian::Bool=true)
    if check_hermitian
        ishermitian(G) || throw(ArgumentError("matrix must be Hermitian"))
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
                adjwgt[adjncy_i] = G.nzval[k]
            end
        end
        xadj[j+1] = xadj[j] + n_rows
    end
    resize!(adjncy, adjncy_i)
    weights && resize!(adjwgt, adjncy_i)
    return Graph(idx_t(N), xadj, adjncy, vwgt, adjwgt)
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
    @check METIS_NodeND(Ref{idx_t}(G.nvtxs), G.xadj, G.adjncy, G.vwgt, options, perm, iperm)
    return perm, iperm
end

"""
    Metis.partition(G, n; alg = :KWAY)

Partition the graph `G` in `n` parts.
The partition algorithm is defined by the `alg` keyword:
 - :KWAY: multilevel k-way partitioning
 - :RECURSIVE: multilevel recursive bisection
"""
partition(G, nparts; alg = :KWAY) = partition(graph(G), nparts, alg = alg)

function partition(G::Graph, nparts::Integer; alg = :KWAY)
    part = Vector{idx_t}(undef, G.nvtxs)
    nparts == 1 && return fill!(part, 1) # https://github.com/JuliaSparse/Metis.jl/issues/49
    edgecut = fill(idx_t(0), 1)
    if alg === :RECURSIVE
        @check METIS_PartGraphRecursive(Ref{idx_t}(G.nvtxs), Ref{idx_t}(1), G.xadj, G.adjncy, G.vwgt, C_NULL, C_NULL,
                                        Ref{idx_t}(nparts), C_NULL, C_NULL, options, edgecut, part)
    elseif alg === :KWAY
        @check METIS_PartGraphKway(Ref{idx_t}(G.nvtxs), Ref{idx_t}(1), G.xadj, G.adjncy, G.vwgt, C_NULL, C_NULL,
                                   Ref{idx_t}(nparts), C_NULL, C_NULL, options, edgecut, part)
    else
        throw(ArgumentError("unknown algorithm $(repr(alg))"))
    end
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
    # METIS_ComputeVertexSeparator segfaults with 1-based indexing
    xadj = G.xadj .- idx_t(1)
    adjncy = G.adjncy .- idx_t(1)
    @check METIS_ComputeVertexSeparator(Ref{idx_t}(G.nvtxs), xadj, adjncy, G.vwgt, options, sepsize, part)
    part .+= 1
    return part
end

# Compatibility for Julias that doesn't support package extensions
if !(isdefined(Base, :get_extension))
    include("../ext/MetisGraphs.jl")
    include("../ext/MetisLightGraphs.jl")
end

end # module
