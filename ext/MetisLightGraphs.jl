# SPDX-License-Identifier: MIT

module MetisLightGraphs

using LightGraphs: AbstractSimpleGraph, ne, nv, outneighbors
using Metis: Metis

"""
    graph(G::LightGraphs.AbstractSimpleGraph; weights=false) :: Metis.Graph

Construct the 1-based CSR representation of the (unweighted) graph `G`.
If `weights=true` all edges get the weight 1.
"""
Metis.graph(G::AbstractSimpleGraph; kwargs...) = graph(G; kwargs...)

get_weight(G, i, j) = 1

include("GraphsCommon.jl")

end # module MetisLightGraphs
