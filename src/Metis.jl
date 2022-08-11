module Metis

using SparseArrays
using LinearAlgebra
import LightGraphs, Graphs
using METIS_jll: libmetis

# Metis C API
include("metis_h.jl")
const options = fill(idx_t(-1), METIS_NOPTIONS)
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
            if i != j # don't include diagonal elements
                n_rows += 1
                adjncy_i += 1
                adjncy[adjncy_i] = i
                if weights
                    adjwgt[adjncy_i] = G.nzval[k]
                end
            end
        end
        xadj[j+1] = xadj[j] + n_rows
    end
    resize!(adjncy, adjncy_i)
    weights && resize!(adjwgt, adjncy_i)
    return Graph(idx_t(N), xadj, adjncy, vwgt, adjwgt)
end

"""
    graph(G::Graphs.AbstractSimpleGraph)
    graph(G::LightGraphs.AbstractSimpleGraph)

Construct the 1-based CSR representation of the `(Light)Graphs` graph `G`.
"""
function graph(G::Union{Graphs.AbstractSimpleGraph, LightGraphs.AbstractSimpleGraph})
    N = nv(G)
    xadj = Vector{idx_t}(undef, N+1)
    xadj[1] = 1
    adjncy = Vector{idx_t}(undef, 2*ne(G))
    adjncy_i = 0
    for j in 1:N
        ne = 0
        for i in outneighbors(G, j)
            ne += 1
            adjncy_i += 1
            adjncy[adjncy_i] = i
        end
        xadj[j+1] = xadj[j] + ne
    end
    resize!(adjncy, adjncy_i)
    return Graph(idx_t(N), xadj, adjncy)
end

# Legacy support for LightGraphs
for mod in (:LightGraphs, :Graphs); @eval begin
    nv(G::$(mod).AbstractSimpleGraph) = $(mod).nv(G)
    ne(G::$(mod).AbstractSimpleGraph) = $(mod).ne(G)
    outneighbors(G::$(mod).AbstractSimpleGraph, j) = $(mod).outneighbors(G, j)
end end

"""
    graph(elements::AbstractMatrix{<:Integer}, ncommon=1)

Construct the 1-based CSR representation of the mesh described by the
connectivity matrix `elements`.

`elements[i,j]` contains the `i`-th node of the `j`-th element. Thus the size of
the first dimension `elements` is equal to the number of nodes per element, and
the size of the second dimension is the number of elements.

The nodes are assumed to be numbered consecutively, starting from 1.
"""
function graph(elements::AbstractMatrix{<:Integer}, ncommon=1)
    # Assume the number of nodes in the mesh is equal to the highest node
    # number.
    nnodes = maximum(elements)

    # Number of nodes per element.
    n_npe = size(elements, 1)

    # Number of elements in the mesh.
    ne = size(elements, 2)

    # Get the element connectivity data in the form METIS needs.
    eptr = convert(Array{Metis.idx_t}, collect(0:ne) .* n_npe)
    eind = convert(Array{Metis.idx_t}, elements[:])

    # With 1-based indexing, need to add 1 to eptr.
    eptr .+= 1

    # Create some pointers and stuff that will be needed for the METIS library call.
    ne = Metis.idx_t(ne)
    nn = Metis.idx_t(nnodes)
    ncommon = Metis.idx_t(ncommon)
    num_flag = Metis.options[Metis.METIS_OPTION_NUMBERING]
    xadj = [Ptr{Metis.idx_t}()]
    adjncy = [Ptr{Metis.idx_t}()]

    # Do it!
    Metis.METIS_MeshToDual(ne, nn, eptr, eind, ncommon, num_flag, xadj, adjncy)

    # Create a Julian copy of the xadj array. The length of xadj is one more
    # than the number of graph verticies (here, the number of mesh elements).
    xadj_out = copy(unsafe_wrap(Vector{Metis.idx_t}, xadj[1], ne+1))

    # Free the METIS memory.
    Metis.METIS_Free(xadj[1])

    # Create a Julian copy of the adjncy array.  The length of adjncy is twice
    # the number of graph edges. Turns out that's the last entry in xadj minus
    # 1 if we're using one-based indices.
    adjncy_out = copy(unsafe_wrap(Vector{Metis.idx_t}, adjncy[1], xadj_out[end]-1))

    # Free the METIS memory.
    Metis.METIS_Free(adjncy[1])

    return Graph(ne, xadj_out, adjncy_out)
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
    METIS_NodeND(G.nvtxs, G.xadj, G.adjncy, G.vwgt, options, perm, iperm)
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
    edgecut = fill(idx_t(0), 1)
    if alg === :RECURSIVE
        METIS_PartGraphRecursive(G.nvtxs, idx_t(1), G.xadj, G.adjncy, G.vwgt, C_NULL, C_NULL,
                                 idx_t(nparts), C_NULL, C_NULL, options, edgecut, part)
    elseif alg === :KWAY
        METIS_PartGraphKway(G.nvtxs, idx_t(1), G.xadj, G.adjncy, G.vwgt, C_NULL, C_NULL,
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
    # METIS_ComputeVertexSeparator segfaults with 1-based indexing
    xadj = G.xadj .- idx_t(1)
    adjncy = G.adjncy .- idx_t(1)
    METIS_ComputeVertexSeparator(G.nvtxs, xadj, adjncy, G.vwgt, options, sepsize, part)
    part .+= 1
    return part
end

end # module
