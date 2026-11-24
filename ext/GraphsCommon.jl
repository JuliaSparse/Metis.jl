# SPDX-License-Identifier: MIT

# The extension module loading this file is required to define `nv`, `ne`,
# `outneighbors`, `is_directed`, `has_edge` and `get_weight` following the Graphs.jl and
# SimpleWeightedGraphs.jl meanings.

using Metis.LibMetis: idx_t
using Metis: Graph, edge_weight

function graph(G; weights::Bool=false, check::Bool=true)
    N = nv(G)
    xadj = Vector{idx_t}(undef, N+1)
    xadj[1] = 1
    adjncy = Vector{idx_t}(undef, 2*ne(G))
    vwgt = C_NULL # TODO: Vertex weights could be passed as an input argument
    adjwgt = weights ? Vector{idx_t}(undef, 2*ne(G)) : C_NULL
    adjncy_i = 0
    # Metis requires every edge to have a matching edge of the same weight in the other
    # direction, which is always the case for an undirected graph
    check_reverse = check && is_directed(G)
    for j in 1:N
        ne = 0
        for i in outneighbors(G, j)
            i == j && continue # skip self edges
            if check_reverse && !has_edge(G, i, j)
                throw(ArgumentError(
                    "Metis requires an edge in both directions, but the directed graph has an " *
                    "edge $(j) → $(i) without an edge $(i) → $(j)"
                ))
            end
            ne += 1
            adjncy_i += 1
            adjncy[adjncy_i] = i
            if weights
                w = get_weight(G, j, i)
                if check_reverse && w != get_weight(G, i, j)
                    throw(ArgumentError(
                        "Metis requires the same weight in both directions, but the directed " *
                        "graph has weight $(w) for edge $(j) → $(i) and weight " *
                        "$(get_weight(G, i, j)) for edge $(i) → $(j)"
                    ))
                end
                adjwgt[adjncy_i] = edge_weight(w, j, i)
            end
        end
        xadj[j+1] = xadj[j] + ne
    end
    resize!(adjncy, adjncy_i)
    weights && resize!(adjwgt, adjncy_i)
    return Graph(idx_t(N), xadj, adjncy, vwgt, adjwgt)
end
