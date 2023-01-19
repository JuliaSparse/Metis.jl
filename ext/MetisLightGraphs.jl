# SPDX-License-Identifier: MIT

module MetisLightGraphs

using LightGraphs: AbstractSimpleGraph, ne, nv, outneighbors
using Metis: Metis

"""
    graph(G::LightGraphs.AbstractSimpleGraph) :: Metis.Graph

Construct the 1-based CSR representation of the (unweighted) graph `G`.
"""
Metis.graph(G::AbstractSimpleGraph; kwargs...) = graph(G; kwargs...)

include("GraphsCommon.jl")

end # module MetisLightGraphs
