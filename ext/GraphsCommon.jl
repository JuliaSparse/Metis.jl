# SPDX-License-Identifier: MIT

# The extension module loading this file is required to define `nv`, `ne`, and
# `outneighbors` following the Graphs.jl meanings.

using Metis.LibMetis: idx_t
using Metis: Graph

function graph(G; weights::Bool=false, kwargs...)
    N = nv(G)
    xadj = Vector{idx_t}(undef, N+1)
    xadj[1] = 1
    adjncy = Vector{idx_t}(undef, 2*ne(G))
    vwgt = C_NULL # TODO: Vertex weights could be passed as an input argument
    adjwgt = weights ? Vector{idx_t}(undef, 2*ne(G)) : C_NULL
    adjncy_i = 0
    for j in 1:N
        ne = 0
        for i in outneighbors(G, j)
            i == j && continue # skip self edges
            ne += 1
            adjncy_i += 1
            adjncy[adjncy_i] = i
            if weights
                weight = get_weight(G, i, j)
                if !(isinteger(weight) && weight > 0)
                    error("weights must be positive integers, got weight $(weight) for edge ($(i), $(j))")
                end
                adjwgt[adjncy_i] = weight
            end
        end
        xadj[j+1] = xadj[j] + ne
    end
    resize!(adjncy, adjncy_i)
    weights && resize!(adjwgt, adjncy_i)
    return Graph(idx_t(N), xadj, adjncy, vwgt, adjwgt)
end
