## Interface to the Metis graph-partitioning algorithms,
##   http://glaros.dtc.umn.edu/gkhome/views/metis
## Version 5.1.0 or later of the metis library is required.
## The library should be compiled with IDXTYPEWIDTH = 32 and 
## REALTYPEWIDTH = 32 (the default values)

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

"""
```
partGraphKway(G, nparts,
              [, adjwgt, vwgt, vsize,
               tpwgts, ubvec, options]) -> objval, part
```

Partition a graph with multilevel k-way partitioning. See Metis 
documentation for more detail.

Inputs:

* `G` : Undirected graph. Can be an adjacency matrix (sparse or dense), a Graphs package adjacency list, or a LightGraphs package graph.
* `nparts` : The number of parts to partition the graph.
* `adjwgt` : Boolean indicating if edges are weighted. If true, `G` must be an adjacency matrix.
* `vwgt` : The weights of the vertices.
* `vsize` : The size of the vertices for computing the total communication volume.
* `tpwgts` : Array of size `nparts` x `ncon` that specifies the desired weight for each partition and constraint.
* `ubvec` : Array of size `ncon` that specifies the allowed load imbalance tolerance for each constraint.
* `options` : Array of options.

Outputs:

* `objval` : The edge-cut or the total communication volume of the partitioning solution.
* `part` : The partition vector of the graph.

"""
    function partGraphKway(G, nparts::Integer;
                           adjwgt::Bool=false,
                           vwgt::Array=[], vsize::Array=[],
                           tpwgts::Array=[], ubvec::Array=[],
                           options::Array=[])

        # Get adjacency structure of graph
        if adjwgt
            typeof(G)<:SparseMatrixCSC || error("weighted graphs must be represented in CSC format")
            n, xadj, adjncy, _adjwgt = metisform_weighted(G)
        else
            n, xadj, adjncy = metisform(G)
            _adjwgt = C_NULL
        end

        # Check parameters are valid
        nparts > 1 || error("nparts must be greater than one")
        length(vwgt)==0 || size(vwgt)==(n,) || error("vwgt must have n entries")
        length(vsize)==0 || size(vsize)==(n,) || error("vsize must have n entries")
        if length(tpwgts)>0
            ncon = convert(Cint, size(tpwgts)[1])
            size(tpwgts)==(ncon,nparts) || error("tpwgts must be an ncon x nparts array")
        else
            ncon = one(Cint)
        end
        length(ubvec)==0 || size(ubvec)==(ncon,) || error("ubvec must have ncon entries")
        length(options)==0 || size(options)==(METIS_NOPTIONS,) || error("options must have METIS_NOPTIONS entries")

        # Initialize partition parameters
        _vwgt = (length(vwgt)>0) ? round(Cint, vwgt) : C_NULL
        _vsize = (length(vsize)>0) ? round(Cint, vsize) : C_NULL
        _tpwgts = (length(tpwgts)>0) ? convert(Array{Cfloat}, tpwgts) : C_NULL
        _ubvec = (length(ubvec)>0) ? convert(Array{Cfloat}, ubvec) : C_NULL
        _options = (length(options)>0) ? convert(Array{Cint}, options) : C_NULL

        # Allocate memory for outputs
        part = Array(Cint, n)
        objval = zeros(Cint, 1)

        # Call Metis partitioner
        err = ccall((:METIS_PartGraphKway,libmetis), Cint,
                    (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
                     Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cfloat}, Ptr{Cfloat},
                     Ptr{Cint}, Ptr{Cint}, Ptr{Cint}),
                    &n, &ncon, xadj, adjncy,
                    _vwgt, _vsize, _adjwgt, &convert(Cint,nparts),
                    _tpwgts, _ubvec, _options, objval, part)
        err == METIS_OK || error("METIS_PartGraphKWay returned error code $err")

        # Return results
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
