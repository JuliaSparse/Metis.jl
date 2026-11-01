# SPDX-License-Identifier: MIT

module MetisLightGraphs

using LightGraphs: AbstractSimpleGraph, has_edge, is_directed, ne, nv, outneighbors
using Metis: Metis

"""
    Metis.graph(G::LightGraphs.AbstractSimpleGraph; weights=false, check=true) :: Metis.Graph

Construct the 1-based CSR representation of the (unweighted) graph `G`.
If `weights=true` all edges get the weight 1.
If `check=true` an `ArgumentError` is thrown if `G` is directed and has an edge without a
matching edge in the other direction, since Metis requires edges in both directions. With
`check=false` this isn't checked.
"""
Metis.graph(G::AbstractSimpleGraph; weights::Bool=false, check::Bool=true) = graph(G; weights, check)

get_weight(G, i, j) = 1

include("GraphsCommon.jl")

end # module MetisLightGraphs
