# SPDX-License-Identifier: MIT

module MetisSimpleWeightedGraphs

using Graphs: ne, nv, outneighbors
using Metis: Metis
using SimpleWeightedGraphs: AbstractSimpleWeightedGraph, get_weight

"""
    graph(G::SimpleWeightedGraphs.AbstractSimpleGraph; weights=true) :: Metis.Graph

Construct the 1-based CSR representation of the weighted graph `G`.
"""
function Metis.graph(G::AbstractSimpleWeightedGraph; weights::Bool=true, kwargs...)
    return graph(G; weights=weights, kwargs...)
end

# METIS_NodeND and METIS_ComputeVertexSeparator don't take edge weights, so don't
# require the weights of `G` to be valid for Metis.
Metis.permutation(G::AbstractSimpleWeightedGraph) = Metis.permutation(graph(G; weights=false))
Metis.separator(G::AbstractSimpleWeightedGraph) = Metis.separator(graph(G; weights=false))

include("GraphsCommon.jl")

end # module MetisSimpleWeightedGraphs
