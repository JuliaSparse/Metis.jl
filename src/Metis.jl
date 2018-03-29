__precompile__(true)

module Metis

import Compat: undef

# Load libmetis with BinaryProvider
__init__() = check_deps()
let depsfile = joinpath(@__DIR__, "..", "deps", "deps.jl")
    if isfile(depsfile)
        include(depsfile)
    else
        error("$(depsfile) does not exist, Please re-run Pkg.build(\"Metis\"), and restart Julia.")
    end
end

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
    Metis.graph(G::SparseMatrixCSC)

Construct the 1-based CSR representation of the sparse matrix `G`.
"""
function graph(G::SparseMatrixCSC)
    ishermitian(G) || throw(ArgumentError("matrix must be Hermitian"))
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

end # module
