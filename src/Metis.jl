## Interface to the Metis graph-partitioning algorithms,
##   http://glaros.dtc.umn.edu/gkhome/views/metis
## Version 5.1.0 or later of the metis library is required.
## The library should be compiled with IDXTYPEWIDTH = 32 (the default)

module Metis
    if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
        include("../deps/deps.jl")      # define libmetis and check it can be dlopen'd
    else
        error("Metis package not properly installed. Please run Pkg.build(\"Metis\")")
    end

    using Graphs,Compat                 # for AdjacencyList types
    using LightGraphs                   # metisform
    export
        nodeND,                         # determine fill-reducing permutation
        vertexSep,                      # single separator
        partGraphKway,
        partGraphRecursive

    include("metis_h.jl")               # define constants
    include("util.jl")                  # metisform and testgraph functions

    const metis_options = -ones(Int32, METIS_NOPTIONS) # defaults indicated by -1

    function nodeND(x,verbose::Integer=0)
        n,xadj,adjncy = metisform(x)
        metis_options[METIS_OPTION_DBGLVL] = verbose
        perm = Array(Cint, n)
        iperm = Array(Cint, n)
        err = ccall((:METIS_NodeND,libmetis), Cint,
                    (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
                     Ptr{Cint}, Ptr{Cint}, Ptr{Cint}),
                    &n, xadj, adjncy, C_NULL, metis_options, perm, iperm)
        err == METIS_OK || error("METIS_NodeND returned error code $err")
        perm .+ one(Cint), iperm .+ one(Cint)
    end

    function vertexSep(x,verbose::Integer=0)
        n,xadj,adjncy = metisform(x)

        metis_options[METIS_OPTION_DBGLVL] = verbose
        sepSize = zeros(Cint, 1)
        part = Array(Cint, n)

        err = ccall((:METIS_ComputeVertexSeparator,libmetis), Cint,
                    (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
                     Ptr{Cint}, Ptr{Cint}, Ptr{Cint}),
                    &n, xadj, adjncy, C_NULL, metis_options,
                    sepSize, part)
        err == METIS_OK || error("METIS_ComputeVertexSeparator returned error code $err")

        sizes = zeros(Cint,3)
        for i=1:n
          sizes[part[i]+1] += 1
        end
        sizes, part
    end

    function partGraphKway(x, nparts::Integer)
        n, xadj, adjncy = metisform(x)
        part = Array(Cint, n)
        objval = zeros(Cint, 1)
        err = ccall((:METIS_PartGraphKway,libmetis), Int32,
                    (Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
                     Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
                     Ptr{Int32}, Ptr{Int32}, Ptr{Int32}),
                    &n, &one(Int32), xadj, adjncy, C_NULL, C_NULL, C_NULL, &convert(Int32,nparts),
                    C_NULL, C_NULL, metis_options, objval, part)
        err == METIS_OK || error("METIS_PartGraphKWay returned error code $err")
        objval[1], part .+ one(Cint)
    end

    function partGraphRecursive(x, nparts::Integer)
        n, xadj, adjncy = metisform(x)
        part = Array(Int32, n)
        objval = zeros(Int32, 1)
        err = ccall((:METIS_PartGraphKway,libmetis), Int32,
                    (Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
                     Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
                     Ptr{Int32}, Ptr{Int32}, Ptr{Int32}),
                    &n, &one(Int32), xadj, adjncy, C_NULL, C_NULL, C_NULL, &convert(Int32,nparts),
                    C_NULL, C_NULL, metis_options, objval, part)
        err == METIS_OK || error("METIS_PartGraphKWay returned error code $err")
        objval[1], part .+ one(Cint)
    end
end
