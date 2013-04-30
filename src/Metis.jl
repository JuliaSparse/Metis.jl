## Interface to the Metis graph-partitioning algorithms,
##   http://glaros.dtc.umn.edu/gkhome/views/metis
## Version 5.1.0 or later of the metis library is required.
## The library should be compiled with IDXTYPEWIDTH = 32 (the default)

module Metis
    using Graphs

    export
        nodeND,
        partGraphKway,
        partGraphRecursive

    include("metis_h.jl")

    const metis_options = -ones(Int32, METIS_NOPTIONS) # defaults indicated by -1
    metis_options[METIS_OPTION_NUMBERING] = int32(1)   # 1-based numbering of vertices

    ## Create the 1-based CSR representation of the adjacency list
    function mkadj{T<:Integer}(al::GenericAdjacencyList{T,Range1{T},Vector{Vector{T}}})
        if is_directed(al) error("Metis functions require undirected graphs") end
        if first(al.vertices) != 1 error("Vertices must be numbered from 1") end
        length(al.adjlist), int32(cumsum(vcat(1, map(length, al.adjlist)))), int32(vcat(al.adjlist...))
    end
        
    function nodeND{T<:Integer}(al::GenericAdjacencyList{T,Range1{T},Vector{Vector{T}}})
        n, xadj, adjncy = mkadj(al)
        perm = Array(Int32, n)
        iperm = Array(Int32, n)
        err = ccall((:METIS_NodeND,:libmetis), Int32,
                    (Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
                     Ptr{Int32}, Ptr{Int32}),
                    &int32(n), xadj, adjncy, C_NULL, metis_options, perm, iperm)
        if (err != METIS_OK) error("METIS_NodeND returned error code $err") end
        perm, iperm
    end

    function partGraphKway{T<:Integer}(al::GenericAdjacencyList{T,Range1{T},Vector{Vector{T}}},
                                       nparts::Integer)
        n, xadj, adjncy = mkadj(al)
        part = Array(Int32, n)
        objval = zeros(Int32, 1)
        err = ccall((:METIS_PartGraphKway,:libmetis), Int32,
                    (Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
                     Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
                     Ptr{Int32}, Ptr{Int32}, Ptr{Int32}),
                    &int32(n), &one(Int32), xadj, adjncy, C_NULL, C_NULL, C_NULL, &int32(nparts),
                    C_NULL, C_NULL, metis_options, objval, part)
        if (err != METIS_OK) error("METIS_PartGraphKWay returned error code $err") end
        objval, part
    end

    function partGraphRecursive{T<:Integer}(al::GenericAdjacencyList{T,Range1{T},Vector{Vector{T}}},
                                            nparts::Integer)
        n, xadj, adjncy = mkadj(al)
        part = Array(Int32, n)
        objval = zeros(Int32, 1)
        err = ccall((:METIS_PartGraphKway,:libmetis), Int32,
                    (Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
                     Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
                     Ptr{Int32}, Ptr{Int32}, Ptr{Int32}),
                    &int32(n), &one(Int32), xadj, adjncy, C_NULL, C_NULL, C_NULL, &int32(nparts),
                    C_NULL, C_NULL, metis_options, objval, part)
        if (err != METIS_OK) error("METIS_PartGraphKWay returned error code $err") end
        objval, part
    end

    function testgraph(nm::ASCIIString)
        pathnm = joinpath(Pkg.dir(), "Metis", "deps", "metis-5.1.0", "graphs", string(nm, ".graph"))
        ff = open(pathnm, "r")
        nvert, nedge = map(int, split(readline(ff)))
        adjlist = Array(Vector{Int32}, nvert)
        for i in 1:nvert adjlist[i] = map(int32, split(readline(ff))) end
        GenericAdjacencyList{Int32,Range1{Int32},Vector{Vector{Int32}}}(false,
                                                                        one(Int32):int32(nvert),
                                                                        nedge,
                                                                        adjlist)
    end
end
