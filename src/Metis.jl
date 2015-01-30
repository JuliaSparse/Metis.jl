## Interface to the Metis graph-partitioning algorithms,
##   http://glaros.dtc.umn.edu/gkhome/views/metis
## Version 5.1.0 or later of the metis library is required.
## The library should be compiled with IDXTYPEWIDTH = 32 (the default)

module Metis
    include("../deps/deps.jl")          # load library

    using Graphs                        # for AdjacencyList types

    export
        nodeND,                         # determine fill-reducing permutation
        nodeND!,                        # mutating version
        vertexSep,                      # single separator
        vertexSep!,
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
        
    function nodeND!{Tv}(m::SparseMatrixCSC{Tv,Cint},verbose::Integer)
                                        # check symmetry of structure
        nz = m.nzval                    # nzval made uniform
        @inbounds for i in 1:length(nz)
            nz[i] = one(Tv)
        end
        issym(m) || ishermitian(m) || error("m must be symmetric or Hermitian")
                                        # drop diagonal elements
        Base.SparseMatrix.fkeep!(m,(i,j,x,other) -> i != j, None)
        metis_options[METIS_OPTION_DBGLVL] = verbose
        n = convert(Cint, m.n)
        perm = Array(Cint, n)
        iperm = Array(Cint, n)
        err = ccall((:METIS_NodeND,libmetis), Cint,
                    (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
                     Ptr{Cint}, Ptr{Cint}, Ptr{Cint}),
                    &n, m.colptr, m.rowval, C_NULL, metis_options, perm, iperm)
        err == METIS_OK || error("METIS_NodeND returned error code $err")
        perm, iperm
    end

    nodeND{Tv}(m::SparseMatrixCSC{Tv,Cint},verbose::Integer) = nodeND!(copy(m),verbose)

    nodeND!{Tv}(m::SparseMatrixCSC{Tv,Cint}) = nodeND!(m,0) # default to no output

    nodeND{Tv}(m::SparseMatrixCSC{Tv,Cint}) = nodeND!(copy(m))

    nodeND{Tv,Ti}(m::SparseMatrixCSC{Tv,Ti}) = nodeND!(convert(SparseMatrixCSC{Tv,Cint},m))

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

    function vertexSep!{Tv}(m::SparseMatrixCSC{Tv,Cint},verbose::Integer)
                                        # check symmetry of structure
        nz = m.nzval                    # nzval made uniform
        @inbounds for i in 1:length(nz)
            nz[i] = one(Tv)
        end
        issym(m) || ishermitian(m) || error("m must be symmetric or Hermitian")
                                        # drop diagonal elements
        Base.SparseMatrix.fkeep!(m,(i,j,x,other) -> i != j, None)
        metis_options[METIS_OPTION_DBGLVL] = verbose
        n = convert(Cint, m.n)
        sepSize = Cint()
        part = Array(Cint, n)
        err = ccall((:METIS_ComputeVertexSeparator,libmetis), Cint,
                    (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
                     Ptr{Cint}, Ptr{Cint}, Ptr{Cint}),
                    &n, m.colptr, m.rowval, C_NULL, metis_options, 
                    &sepSize, part)
        err == METIS_OK || error("METIS_ComputeVertexSeparator returned error code $err")
        sepSize, part
    end

    vertexSep{Tv}(m::SparseMatrixCSC{Tv,Cint},verbose::Integer) = vertexSep!(copy(m),verbose)

    vertexSep!{Tv}(m::SparseMatrixCSC{Tv,Cint}) = vertexSep!(m,0) # default to no output

    vertexSep{Tv}(m::SparseMatrixCSC{Tv,Cint}) = vertexSep!(copy(m))

    vertexSep{Tv,Ti}(m::SparseMatrixCSC{Tv,Ti}) = vertexSep!(convert(SparseMatrixCSC{Tv,Cint},m))

    function vertexSep{T<:Integer}(al::GenericAdjacencyList{T,Range1{T},Vector{Vector{T}}})
        n, xadj, adjncy = mkadj(al)
        sepSize = Cint()
        part = Array(Int32, n)
        err = ccall((:METIS_ComputeVertexSeparator,libmetis), Int32,
                    (Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
                     Ptr{Int32}, Ptr{Int32}),
                    &int32(n), xadj, adjncy, C_NULL, metis_options, 
                    &sepSize, part)
        if (err != METIS_OK) error("METIS_ComputeVertexSeparator returned error code $err") end
        sepSize, part
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
