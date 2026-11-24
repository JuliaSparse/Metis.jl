# SPDX-License-Identifier: MIT

module MetisSimpleWeightedGraphs

using Graphs: has_edge, is_directed, ne, nv, outneighbors
using Metis: Metis
using SimpleWeightedGraphs: AbstractSimpleWeightedGraph, get_weight

"""
    Metis.graph(G::SimpleWeightedGraphs.AbstractSimpleWeightedGraph; weights=true, check=true) :: Metis.Graph

Construct the 1-based CSR representation of the weighted graph `G`.
If `check=true` an `ArgumentError` is thrown if `G` is directed and has an edge without a
matching edge in the other direction, or, if `weights=true`, with a different weight than the
matching edge, since Metis requires edges in both directions with the same weight. With
`check=false` this isn't checked.
"""
function Metis.graph(G::AbstractSimpleWeightedGraph; weights::Bool=true, check::Bool=true)
    return graph(G; weights, check)
end

# METIS_NodeND and METIS_ComputeVertexSeparator don't take edge weights, so don't
# require the weights of `G` to be valid for Metis.
Metis.permutation(G::AbstractSimpleWeightedGraph; options::Metis.Options = Metis.Options()) =
    Metis.permutation(graph(G; weights=false); options)
Metis.separator(G::AbstractSimpleWeightedGraph; options::Metis.Options = Metis.Options()) =
    Metis.separator(graph(G; weights=false); options)

include("GraphsCommon.jl")

end # module MetisSimpleWeightedGraphs
