# SPDX-License-Identifier: MIT

# The extension module loading this file is required to define `nv`, `ne`, and
# `outneighbors` following the Graphs.jl meanings.

using Metis.LibMetis: idx_t
using Metis: Graph

function graph(G)
    N = nv(G)
    xadj = Vector{idx_t}(undef, N+1)
    xadj[1] = 1
    adjncy = Vector{idx_t}(undef, 2*ne(G))
    adjncy_i = 0
    for j in 1:N
        ne = 0
        for i in outneighbors(G, j)
            ne += 1
            adjncy_i += 1
            adjncy[adjncy_i] = i
        end
        xadj[j+1] = xadj[j] + ne
    end
    resize!(adjncy, adjncy_i)
    return Graph(idx_t(N), xadj, adjncy)
end
