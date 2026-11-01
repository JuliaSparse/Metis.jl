# SPDX-License-Identifier: MIT

module MetisSimpleWeightedGraphs

using Graphs: ne, nv, outneighbors
using Metis: Metis
using SimpleWeightedGraphs: AbstractSimpleWeightedGraph, get_weight

"""
    Metis.graph(G::SimpleWeightedGraphs.AbstractSimpleWeightedGraph; weights=true) :: Metis.Graph

Construct the 1-based CSR representation of the weighted graph `G`.
"""
function Metis.graph(G::AbstractSimpleWeightedGraph; weights::Bool=true, kwargs...)
    return graph(G; weights=weights, kwargs...)
end

include("GraphsCommon.jl")

end # module MetisSimpleWeightedGraphs
