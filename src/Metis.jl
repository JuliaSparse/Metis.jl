## Interface to the Metis graph-partitioning algorithms,
##   http://glaros.dtc.umn.edu/gkhome/views/metis
## Version 5.1.0 or later of the metis library is required.
## The library should be compiled with IDXTYPEWIDTH = 32 (the default)

module Metis
    if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
        include("../deps/deps.jl")      # load the library and define libmetis
    else
        error("Metis package not properly installed. Please run Pkg.build(\"Metis\")")
    end

    using Graphs                        # for AdjacencyList types

    export
        nodeND,                         # determine fill-reducing permutation
        vertexSep,                      # single separator
        partGraphKway,
        partGraphRecursive

    include("metis_h.jl")

    const metis_options = -ones(Int32, METIS_NOPTIONS) # defaults indicated by -1
    metis_options[METIS_OPTION_NUMBERING] = 1 # 1-based numbering of vertices

    ## Create the 1-based CSR representation of the adjacency list
    function mkadj{T<:Integer}(al::GenericAdjacencyList{T,Range1{T},Vector{Vector{T}}})
        !is_directed(al) || error("Metis functions require undirected graphs")
        first(al.vertices) == 1 || error("Vertices must be numbered from 1")
        length(al.adjlist), int32(cumsum(vcat(1, map(length, al.adjlist)))), int32(vcat(al.adjlist...))
    end
        
    function nodeND(m::SparseMatrixCSC,verbose::Integer=0)
        issym(m) || isHermitian(m) || error("m must be symmetric or Hermitian")

        ## copy m.rowval and m.colptr to Int32 vectors dropping diagonal elements
        rowvs = sizehint(Cint[],nnz(m))
        colpt = ones(Cint,m.n+1)
        for col in 1:m.n
            colcount = 0
            for k in m.colptr[col] : (m.colptr[col+1] - 1)
                rv = m.rowval[k]
                if rv != col
                    colcount += 1
                    push!(rowvs,rv)
                end
            end
            colpt[col+1] = colpt[col] + colcount
        end

        metis_options[METIS_OPTION_DBGLVL] = verbose
        perm = Array(Cint, m.n)
        iperm = Array(Cint, m.n)
        err = ccall((:METIS_NodeND,libmetis), Cint,
                    (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
                     Ptr{Cint}, Ptr{Cint}, Ptr{Cint}),
                    &int32(m.n), colpt, rowvs, C_NULL, metis_options, perm, iperm)
        err == METIS_OK || error("METIS_NodeND returned error code $err")
        perm, iperm
    end

    function nodeND{T<:Integer}(al::GenericAdjacencyList{T,Range1{T},Vector{Vector{T}}})
        n, xadj, adjncy = mkadj(al)
        perm = Array(Int32, n)
        iperm = Array(Int32, n)
        err = ccall((:METIS_NodeND,libmetis), Int32,
                    (Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
                     Ptr{Int32}, Ptr{Int32}),
                    &int32(n), xadj, adjncy, C_NULL, metis_options, perm, iperm)
        if (err != METIS_OK) error("METIS_NodeND returned error code $err") end
        perm, iperm
    end

    function vertexSep{Tv}(m::SparseMatrixCSC{Tv,Cint},verbose::Integer)
        # We must make a copy of m so that we can check for structural symmetry
        # as well as drop entries to satisfy METIS's xadj and adjncy constraints
        mMod = copy(m) 

        # Test for structural symmetry
        nz = mMod.nzval 
        @inbounds for i in 1:length(nz)
            nz[i] = one(Tv)
        end
        # check symmetry of structure
        issym(mMod) || ishermitian(mMod) || 
          error("mMod must be symmetric or Hermitian")

        # Drop self-edges
        Base.SparseMatrix.fkeep!(mMod,(i,j,x,other) -> i != j, None)

        metis_options[METIS_OPTION_DBGLVL] = verbose
        n = convert(Cint, mMod.n)
        sepSize = zeros(Cint, 1)
        part = Array(Cint, n)

        # Since the ParMETIS internal routines do not seem to support one-based
        # numbering, we must manually decrease the entries of mMod.xadj and
        # mMod.rowval (which is destructive)
        for i=1:mMod.colptr[n+1]-1
          mMod.rowval[i] -= 1;
        end
        for i=1:n+1
          mMod.colptr[i] -= 1; 
        end

        err = ccall((:METIS_ComputeVertexSeparator,libmetis), Cint,
                    (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
                     Ptr{Cint}, Ptr{Cint}, Ptr{Cint}),
                    &n, mMod.colptr, mMod.rowval, C_NULL, metis_options, 
                    sepSize, part)
        err == METIS_OK || error("METIS_ComputeVertexSeparator returned error code $err")

        sizes = zeros(Cint,3)
        for i=1:n
          sizes[part[i]+1] += 1
        end
        sizes, part
    end

    vertexSep{Tv,Ti}(m::SparseMatrixCSC{Tv,Ti},verbose::Integer=0) = vertexSep(convert(SparseMatrixCSC{Tv,Cint},m),verbose)

    function vertexSep{T<:Integer}(al::GenericAdjacencyList{T,Range1{T},Vector{Vector{T}}})
        n, xadj, adjncy = mkadj(al)
        sepSize = zeros(Int32, 1)
        part = Array(Int32, n)

        # Since the ParMETIS internal routines do not seem to support one-based
        # numbering, we must manually decrease the entries of xadj and adjncy
        for i=1:length(adjncy)
          adjncy[i] -= 1;
        end
        for i=1:n+1
          xadj[i] -= 1; 
        end

        err = ccall((:METIS_ComputeVertexSeparator,libmetis), Int32,
                    (Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
                     Ptr{Int32}, Ptr{Int32}),
                    &int32(n), xadj, adjncy, C_NULL, metis_options, 
                    sepSize, part)
        if (err != METIS_OK) error("METIS_ComputeVertexSeparator returned error code $err") end

        sizes = zeros(Int32,3)
        for i=1:n
          sizes[part[i]+1] += 1
        end
        sizes, part
    end

    function partGraphKway{T<:Integer}(al::GenericAdjacencyList{T,Range1{T},Vector{Vector{T}}},
                                       nparts::Integer)
        n, xadj, adjncy = mkadj(al)
        part = Array(Int32, n)
        objval = zeros(Int32, 1)
        err = ccall((:METIS_PartGraphKway,libmetis), Int32,
                    (Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
                     Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
                     Ptr{Int32}, Ptr{Int32}, Ptr{Int32}),
                    &int32(n), &one(Int32), xadj, adjncy, C_NULL, C_NULL, C_NULL, &int32(nparts),
                    C_NULL, C_NULL, metis_options, objval, part)
        if (err != METIS_OK) error("METIS_PartGraphKWay returned error code $err") end
        objval[1], part
    end

    function partGraphRecursive{T<:Integer}(al::GenericAdjacencyList{T,Range1{T},Vector{Vector{T}}},
                                            nparts::Integer)
        n, xadj, adjncy = mkadj(al)
        part = Array(Int32, n)
        objval = zeros(Int32, 1)
        err = ccall((:METIS_PartGraphKway,libmetis), Int32,
                    (Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
                     Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
                     Ptr{Int32}, Ptr{Int32}, Ptr{Int32}),
                    &int32(n), &one(Int32), xadj, adjncy, C_NULL, C_NULL, C_NULL, &int32(nparts),
                    C_NULL, C_NULL, metis_options, objval, part)
        if (err != METIS_OK) error("METIS_PartGraphKWay returned error code $err") end
        objval[1], part
    end

    function testgraph(nm::ASCIIString)
        pathnm = Pkg.dir("Metis", "graphs", string(nm, ".graph"))
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
