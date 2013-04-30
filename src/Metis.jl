## Interface to the Metis graph-partitioning algorithms,
##   http://glaros.dtc.umn.edu/gkhome/views/metis
## Version 5.1.0 or later of the metis library is required.
## The library should be compiled with IDXTYPEWIDTH = 32 (the default)

module Metis
    using Graphs

    export nodeND

    include("metis_h.jl")

    const metis_options = -ones(Int32, METIS_NOPTIONS) # defaults indicated by -1
    metis_options[METIS_OPTION_NUMBERING] = int32(1)   # 1-based numbering of vertices

    function nodeND{T<:Integer}(al::GenericAdjacencyList{T,Range1{T},Vector{Vector{T}}})
        if is_directed(al) error("nodeND applies only to undirected graphs") end
        if first(al.vertices) != 1 error("Vertices must be numbered from 1") end
        n = num_vertices(al)
        nz = num_edges(al)
        perm = Array(Int32, n)
        iperm = Array(Int32, n)
        cp = int32(cumsum(vcat(1, map(length, al.adjlist))))
        rv = Array(Int32, nz)
        pos = 1
        for v in al.adjlist, el in v
            rv[pos] = el
            pos += 1
        end
        try 
            err = ccall((:METIS_NodeND,:libmetis), Int32,
                        (Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
                         Ptr{Int32}, Ptr{Int32}),
                        &int32(n), cp, rv, C_NULL, metis_options, perm, iperm)
            if (err != METIS_OK) error("METIS_NodeND returned error code $err") end
        catch
            error("See http://glaros.dtc.umn.edu/gkhome/views/metis to install libmetis")
        end
        perm, iperm
    end
end
