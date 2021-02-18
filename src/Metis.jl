module Metis

using SparseArrays
using LinearAlgebra
import LightGraphs
import SimpleWeightedGraphs
using METIS_jll: libmetis

# Metis C API
include("metis_h.jl")
const options = fill(Cint(-1), METIS_NOPTIONS)
options[METIS_OPTION_NUMBERING] = 1

# Julia interface
"""
    Metis.Graph

1-based CSR representation of a graph as defined in
section 5.5 "Graph data structure" in the Metis manual.
"""
struct Graph
    nvtxs::idx_t
    xadj::Vector{idx_t}
    adjncy::Vector{idx_t}
    vwgt::Vector{idx_t}
    adjwgt::Vector{idx_t}
    Graph(nvtxs, xadj, adjncy) = new(nvtxs, xadj, adjncy)
    Graph(nvtxs, xadj, adjncy, vwgt) = new(nvtxs, xadj, adjncy, vwgt)
    Graph(nvtxs, xadj, adjncy, vwgt, adjwgt) = new(nvtxs, xadj, adjncy, vwgt, adjwgt)
end

"""
    Metis.graph(G::SparseMatrixCSC; check_hermitian=true)

Construct the 1-based CSR representation of the sparse matrix `G`.
If `check_hermitian` is `false` the matrix is not checked for being hermitian
before constructing the graph.
"""
function graph(G::SparseMatrixCSC; check_hermitian=true)
    if check_hermitian
        ishermitian(G) || throw(ArgumentError("matrix must be Hermitian"))
    end
    N = size(G, 1)
    xadj = Vector{idx_t}(undef, N+1)
    xadj[1] = 1
    adjncy = Vector{idx_t}(undef, nnz(G))
    adjncy_i = 0
    @inbounds for j in 1:N
        n_rows = 0
        for k in G.colptr[j] : (G.colptr[j+1] - 1)
            i = G.rowval[k]
            if i != j # don't include diagonal elements
                n_rows += 1
                adjncy_i += 1
                adjncy[adjncy_i] = i
            end
        end
        xadj[j+1] = xadj[j] + n_rows
    end
    resize!(adjncy, adjncy_i)
    return Graph(idx_t(N), xadj, adjncy)
end

"""
    graph(G::LightGraphs.AbstractSimpleGraph)

Construct the 1-based CSR representation of the `LightGraphs` graph `G`.
"""
function graph(G::LightGraphs.AbstractSimpleGraph)
    N = LightGraphs.nv(G)
    xadj = Vector{idx_t}(undef, N+1)
    xadj[1] = 1
    adjncy = Vector{idx_t}(undef, 2*LightGraphs.ne(G))
    adjncy_i = 0
    for j in 1:N
        ne = 0
        for i in LightGraphs.outneighbors(G, j)
            ne += 1
            adjncy_i += 1
            adjncy[adjncy_i] = i
        end
        xadj[j+1] = xadj[j] + ne
    end
    resize!(adjncy, adjncy_i)
    return Graph(idx_t(N), xadj, adjncy)
end

"""
    graph(G::SimpleWeightedGraphs.AbstractSimpleWeightedGraph)

Construct the 1-based CSR representation of the `SimpleWeightedGraphs` graph `G`.
The edge weights should be positive integers.
"""
function graph(G::SimpleWeightedGraphs.AbstractSimpleWeightedGraph)
    N = SimpleWeightedGraphs.nv(G)
    xadj = Vector{idx_t}(undef, N+1)
    xadj[1] = 1
    adjncy = Vector{idx_t}(undef, 2*SimpleWeightedGraphs.ne(G))
    vwgt = ones(idx_t, N)
    adjwgt = Vector{idx_t}(undef, 2*SimpleWeightedGraphs.ne(G))
    adjncy_i = 0
    for j in 1:N
        ne = 0
        for i in SimpleWeightedGraphs.outneighbors(G, j)
            ne += 1
            adjncy_i += 1
            adjncy[adjncy_i] = i
            edgeweight = SimpleWeightedGraphs.get_weight(G, j, i)
            if isinteger(edgeweight) && edgeweight > 0
                adjwgt[adjncy_i] = idx_t(edgeweight)
            else
                throw(ErrorException("edge weight should be a positive integer (edge ($j, $i) with weight $edgeweight)"))
            end
        end
        xadj[j+1] = xadj[j] + ne
    end
    resize!(adjncy, adjncy_i)
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
    vwgt = isdefined(G, :vwgt) ? G.vwgt : C_NULL
    METIS_NodeND(G.nvtxs, G.xadj, G.adjncy, vwgt, options, perm, iperm)
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
    vwgt = isdefined(G, :vwgt) ? G.vwgt : C_NULL
    adjwgt = isdefined(G, :adjwgt) ? G.adjwgt : C_NULL
    edgecut = fill(idx_t(0), 1)
    if alg === :RECURSIVE
        METIS_PartGraphRecursive(G.nvtxs, idx_t(1), G.xadj, G.adjncy, vwgt, C_NULL, adjwgt,
                                 idx_t(nparts), C_NULL, C_NULL, options, edgecut, part)
    elseif alg === :KWAY
        METIS_PartGraphKway(G.nvtxs, idx_t(1), G.xadj, G.adjncy, vwgt, C_NULL, adjwgt,
                            idx_t(nparts), C_NULL, C_NULL, options, edgecut, part)
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
    vwgt = isdefined(G, :vwgt) ? G.vwgt : C_NULL
    # METIS_ComputeVertexSeparator segfaults with 1-based indexing
    xadj = G.xadj .- idx_t(1)
    adjncy = G.adjncy .- idx_t(1)
    METIS_ComputeVertexSeparator(G.nvtxs, xadj, adjncy, vwgt, options, sepsize, part)
    part .+= 1
    return part
end

end # module
