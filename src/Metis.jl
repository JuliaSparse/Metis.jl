module Metis

using SparseArrays
using LinearAlgebra
import LightGraphs
import TriangleMesh
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
    Graph(nvtxs, xadj, adjncy) = new(nvtxs, xadj, adjncy)
    Graph(nvtxs, xadj, adjncy, vwgt) = new(nvtxs, xadj, adjncy, vwgt)
end

"""
    Metis.Mesh

Representation of a mesh as defined in section 5.6 
"Mesh data structure" in the Metis manual.
"""
struct Mesh
    ne::idx_t
    nn::idx_t
    eptr::Vector{idx_t}
    eind::Vector{idx_t}
    Mesh(ne, nn, eptr, eind) = new(ne, nn, eptr, eind)
end

"""
    Mesh(mesh::TriangleMesh.TriMesh)

Constructor of Metis.Mesh using a TriangleMesh.TriMesh.
"""
function Mesh(mesh::TriangleMesh.TriMesh)
    ne = idx_t(mesh.n_cell)
    nn = idx_t(mesh.n_point)
    eptr = Vector{idx_t}(undef, ne+1)
    eind = Vector{idx_t}(undef, length(mesh.cell))
    eptr[:] = [idx_t(3*i) for i in 0:ne]
    eind[:] = reshape(mesh.cell, length(mesh.cell))
    return Mesh(ne, nn, eptr, eind)
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
    edgecut = fill(idx_t(0), 1)
    if alg === :RECURSIVE
        METIS_PartGraphRecursive(G.nvtxs, idx_t(1), G.xadj, G.adjncy, vwgt, C_NULL, C_NULL,
                                 idx_t(nparts), C_NULL, C_NULL, options, edgecut, part)
    elseif alg === :KWAY
        METIS_PartGraphKway(G.nvtxs, idx_t(1), G.xadj, G.adjncy, vwgt, C_NULL, C_NULL,
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

"""
    Metis.mesh_partition(M, n; alg = :DUAL)

Partition the mesh `M` in `n` parts.
The partition algorithm is defined by the `alg` keyword:
 - :DUAL: Partition of the mesh's dual graph
 - :NODAL: Partition of the mesh's nodal graph
"""
function mesh_partition(mesh::TriangleMesh.TriMesh, nparts::Integer; alg = :DUAL)
    M = Mesh(mesh)
    epart = Vector{idx_t}(undef, M.ne)
    npart = Vector{idx_t}(undef, M.nn)
    objval = fill(idx_t(0), 1)
    if alg === :NODAL
        METIS_PartMeshNodal(M.ne, M.nn, M.eptr, M.eind, C_NULL, C_NULL, idx_t(nparts),
                            C_NULL, options, objval, epart, npart)
    elseif alg === :DUAL
        METIS_PartMeshDual(M.ne, M.nn, M.eptr, M.eind, C_NULL, C_NULL, idx_t(2), idx_t(nparts),
                           C_NULL, options, objval, epart, npart)
    else
        throw(ArgumentError("unknown algorithm $(repr(alg))"))
    end
    return epart, npart
end

end # module

