# SPDX-License-Identifier: MIT

module MetisGraphs

using Graphs: AbstractSimpleGraph, ne, nv, outneighbors
using Metis: Metis

"""
    graph(G::Graphs.AbstractSimpleGraph) :: Metis.Graph

Construct the 1-based CSR representation of the (unweighted) graph `G`.
"""
Metis.graph(G::AbstractSimpleGraph; kwargs...) = graph(G; kwargs...)

include("GraphsCommon.jl")

end # module MetisGraphs
