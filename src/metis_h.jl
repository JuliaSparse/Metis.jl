const idx_t = Cint
const real_t = Cfloat

const METIS_NOPTIONS = 40

## Return codes
const METIS_OK           = Cint(1)  # Returned normally
const METIS_ERROR_INPUT  = Cint(-2) # Returned due to erroneous inputs and/or options
const METIS_ERROR_MEMORY = Cint(-3) # Returned due to insufficient memory
const METIS_ERROR        = Cint(-4) # Some other errors

struct MetisError <: Exception code::Cint end
function Base.showerror(io::IO, me::MetisError)
    print(io, "MetisError: ")
    if me.code == Metis.METIS_ERROR_INPUT
        print(io, "input error")
    elseif me.code == Metis.METIS_ERROR_MEMORY
        print(io, "could not allocate the required memory")
    else
        print(io, "unknown error")
    end
    print(io, " (error code $(me.code)).")
end

## Operation type codes
const METIS_OP_PMETIS = Cint(0)
const METIS_OP_KMETIS = Cint(1)
const METIS_OP_OMETIS = Cint(2)

## (1-based) positions in options vector
const METIS_OPTION_PTYPE     = 1
const METIS_OPTION_OBJTYPE   = 2
const METIS_OPTION_CTYPE     = 3
const METIS_OPTION_IPTYPE    = 4
const METIS_OPTION_RTYPE     = 5
const METIS_OPTION_DBGLVL    = 6
const METIS_OPTION_NITER     = 7
const METIS_OPTION_NCUTS     = 8
const METIS_OPTION_SEED      = 9
const METIS_OPTION_NO2HOP    = 10
const METIS_OPTION_MINCONN   = 11
const METIS_OPTION_CONTIG    = 12
const METIS_OPTION_COMPRESS  = 13
const METIS_OPTION_CCORDER   = 14
const METIS_OPTION_PFACTOR   = 15
const METIS_OPTION_NSEPS     = 16
const METIS_OPTION_UFACTOR   = 17
const METIS_OPTION_NUMBERING = 18
const METIS_OPTION_HELP      = 19
const METIS_OPTION_TPWGTS    = 20
const METIS_OPTION_NCOMMON   = 21
const METIS_OPTION_NOOUTPUT  = 22
const METIS_OPTION_BALANCE   = 23
const METIS_OPTION_GTYPE     = 24
const METIS_OPTION_UBVEC     = 25

## Partitioning Schemes
const METIS_PTYPE_RB   = Cint(0)
const METIS_PTYPE_KWAY = Cint(1)

## Graph types for meshes
const METIS_GTYPE_DUAL  = Cint(0)
const METIS_GTYPE_NODAL = Cint(1)

## Coarsening Schemes
const METIS_CTYPE_RM   = Cint(0)
const METIS_CTYPE_SHEM = Cint(1)

## Initial partitioning schemes
const METIS_IPTYPE_GROW    = Cint(0)
const METIS_IPTYPE_RANDOM  = Cint(1)
const METIS_IPTYPE_EDGE    = Cint(2)
const METIS_IPTYPE_NODE    = Cint(3)
const METIS_IPTYPE_METISRB = Cint(4)

## Refinement schemes
const METIS_RTYPE_FM        = Cint(0)
const METIS_RTYPE_GREEDY    = Cint(1)
const METIS_RTYPE_SEP2SIDED = Cint(2)
const METIS_RTYPE_SEP1SIDED = Cint(3)

## Debug levels (bit positions)
const METIS_DBG_INFO       = Cint(1)    # Shows various diagnostic messages
const METIS_DBG_TIME       = Cint(2)    # Perform timing analysis
const METIS_DBG_COARSEN    = Cint(4)    # Show the coarsening progress
const METIS_DBG_REFINE     = Cint(8)    # Show the refinement progress
const METIS_DBG_IPART      = Cint(16)   # Show info on initial partitioning
const METIS_DBG_MOVEINFO   = Cint(32)   # Show info on vertex moves during refinement
const METIS_DBG_SEPINFO    = Cint(64)   # Show info on vertex moves during sep refinement
const METIS_DBG_CONNINFO   = Cint(128)  # Show info on minimization of subdomain connectivity
const METIS_DBG_CONTIGINFO = Cint(256)  # Show info on elimination of connected components
const METIS_DBG_MEMORY     = Cint(2048) # Show info related to wspace allocation

## Types of objectives
const METIS_OBJTYPE_CUT  = Cint(0)
const METIS_OBJTYPE_VOL  = Cint(1)
const METIS_OBJTYPE_NODE = Cint(2)


## Metis C API
function METIS_PartGraphRecursive(nvtxs, ncon, xadj, adjncy, vwgt, vsize, adjwgt, nparts,
                                  tpwgts, ubvec, options, edgecut, part)
    r = ccall((:METIS_PartGraphRecursive, libmetis), Cint,
              (Ref{idx_t}, Ref{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t},
               Ptr{idx_t}, Ref{idx_t}, Ptr{real_t}, Ptr{real_t}, Ptr{idx_t}, Ptr{idx_t},
               Ptr{idx_t}),
              nvtxs, ncon, xadj, adjncy, vwgt, vsize, adjwgt, nparts, tpwgts, ubvec,
              options, edgecut, part)
    r == METIS_OK || throw(MetisError(r))
    return
end

function METIS_PartGraphKway(nvtxs, ncon, xadj, adjncy, vwgt, vsize, adjwgt, nparts,
                             tpwgts, ubvec, options, edgecut, part)
    r = ccall((:METIS_PartGraphKway, libmetis), Cint,
              (Ref{idx_t}, Ref{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t},
               Ptr{idx_t}, Ref{idx_t}, Ptr{real_t}, Ptr{real_t}, Ptr{idx_t}, Ptr{idx_t},
               Ptr{idx_t}),
              nvtxs, ncon, xadj, adjncy, vwgt, vsize, adjwgt, nparts, tpwgts, ubvec,
              options, edgecut, part)
    r == METIS_OK || throw(MetisError(r))
    return
end

# function METIS_MeshToDual(ne, nn, eptr, eind, ncommon, numflag, r_xadj, r_adjncy)
#     r = ccall((:METIS_MeshToDual, libmetis), Cint,
#               (Ref{idx_t}, Ref{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ref{idx_t}, Ref{idx_t},
#                Ptr{idx_t}, Ptr{idx_t}),
#               ne, nn, eptr, eind, ncommon, numflag, r_xadj, r_adjncy)
#     r == METIS_OK || throw(MetisError(r))
#     return
# end

# function METIS_MeshToNodal(ne, nn, eptr, eind, numflag, r_xadj, r_adjncy)
#     r = ccall((:METIS_MeshToNodal, libmetis), Cint,
#               (Ref{idx_t}, Ref{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ref{idx_t}, Ptr{idx_t},
#                Ptr{idx_t}),
#               ne, nn, eptr, eind, numflag, r_xadj, r_adjncy)
#     r == METIS_OK || throw(MetisError(r))
#     return
# end

# function METIS_PartMeshNodal(ne, nn, eptr, eind, vwgt, vsize, nparts, tpwgts, options,
#                              objval, epart, npart)
#     r = ccall((:METIS_PartMeshNodal, libmetis), Cint,
#               (Ref{idx_t}, Ref{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t},
#                Ref{idx_t}, Ptr{real_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}),
#               ne, nn, eptr, eind, vwgt, vsize, nparts, tpwgts, options, objval, epart,
#               npart)
#     r == METIS_OK || throw(MetisError(r))
#     return
# end

# function METIS_PartMeshDual(ne, nn, eptr, eind, vwgt, vsize, ncommon, nparts, tpwgts,
#                             options, objval, epart, npart)
#     r = ccall((:METIS_PartMeshDual, libmetis), Cint,
#               (Ref{idx_t}, Ref{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t},
#                Ref{idx_t}, Ref{idx_t}, Ptr{real_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t},
#                Ptr{idx_t}),
#               ne, nn, eptr, eind, vwgt, vsize, ncommon, nparts, tpwgts, options, objval,
#               epart, npart)
#     r == METIS_OK || throw(MetisError(r))
#     return
# end

function METIS_NodeND(nvtxs, xadj, adjncy, vwgt, options, perm, iperm)
    r = ccall((:METIS_NodeND, libmetis), Cint,
              (Ref{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t},
               Ptr{idx_t}),
              nvtxs, xadj, adjncy, vwgt, options, perm, iperm)
    r == METIS_OK || throw(MetisError(r))
    return
end

# function METIS_Free(ptr)
#     r = ccall((:METIS_Free, libmetis), Cint, (Ptr{idx_t},), ptr)
#     r == METIS_OK || throw(MetisError(r))
#     return
# end

# function METIS_SetDefaultOptions(options)
#     r = ccall((:METIS_SetDefaultOptions, libmetis), Cint, (Ptr{idx_t},), options)
#     r == METIS_OK || throw(MetisError(r))
#     return
# end

# function METIS_NodeNDP(nvtxs, xadj, adjncy, vwgt, npes, options, perm, iperm, sizes)
#     r = ccall((:METIS_NodeNDP, libmetis), Cint,
#               (Ref{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t},
#                Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}),
#               nvtxs, xadj, adjncy, vwgt, npes, options, perm, iperm, sizes)
#     r == METIS_OK || throw(MetisError(r))
#     return
# end

function METIS_ComputeVertexSeparator(nvtxs, xadj, adjncy, vwgt, options, sepsize, part)
    r = ccall((:METIS_ComputeVertexSeparator, libmetis), Cint,
              (Ref{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t},
               Ptr{idx_t}),
              nvtxs, xadj, adjncy, vwgt, options, sepsize, part)
    r == METIS_OK || throw(MetisError(r))
    return
end

# function METIS_NodeRefine(nvtxs, xadj, vwgt, adjncy, where, hmarker, ubfactor)
#     r = ccall((:METIS_NodeRefine, libmetis), Cint,
#               (Ref{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t},
#                Ptr{real_t}),
#               nvtxs, xadj, vwgt, adjncy, where, hmarker, ubfactor)
#     r == METIS_OK || throw(MetisError(r))
#     return
# end
