# SPDX-License-Identifier: MIT

module MetisGraphs

using Graphs: AbstractSimpleGraph, ne, nv, outneighbors
using Metis: Metis

"""
    graph(G::Graphs.AbstractSimpleGraph; weights=false) :: Metis.Graph

Construct the 1-based CSR representation of the (unweighted) graph `G`.
If `weights=true` all edges get the weight 1.
"""
Metis.graph(G::AbstractSimpleGraph; kwargs...) = graph(G; kwargs...)

get_weight(G, i, j) = 1

include("GraphsCommon.jl")

end # module MetisGraphs
