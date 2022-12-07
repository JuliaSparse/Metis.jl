using METIS_jll
export METIS_jll

using CEnum



const idx_t = Int32

const real_t = Cfloat

function METIS_PartGraphRecursive(nvtxs, ncon, xadj, adjncy, vwgt, vsize, adjwgt, nparts, tpwgts, ubvec, options, edgecut, part)
    ccall((:METIS_PartGraphRecursive, libmetis), Cint, (Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{real_t}, Ptr{real_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}), nvtxs, ncon, xadj, adjncy, vwgt, vsize, adjwgt, nparts, tpwgts, ubvec, options, edgecut, part)
end

function METIS_PartGraphKway(nvtxs, ncon, xadj, adjncy, vwgt, vsize, adjwgt, nparts, tpwgts, ubvec, options, edgecut, part)
    ccall((:METIS_PartGraphKway, libmetis), Cint, (Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{real_t}, Ptr{real_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}), nvtxs, ncon, xadj, adjncy, vwgt, vsize, adjwgt, nparts, tpwgts, ubvec, options, edgecut, part)
end

function METIS_MeshToDual(ne, nn, eptr, eind, ncommon, numflag, r_xadj, r_adjncy)
    ccall((:METIS_MeshToDual, libmetis), Cint, (Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{Ptr{idx_t}}, Ptr{Ptr{idx_t}}), ne, nn, eptr, eind, ncommon, numflag, r_xadj, r_adjncy)
end

function METIS_MeshToNodal(ne, nn, eptr, eind, numflag, r_xadj, r_adjncy)
    ccall((:METIS_MeshToNodal, libmetis), Cint, (Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{Ptr{idx_t}}, Ptr{Ptr{idx_t}}), ne, nn, eptr, eind, numflag, r_xadj, r_adjncy)
end

function METIS_PartMeshNodal(ne, nn, eptr, eind, vwgt, vsize, nparts, tpwgts, options, objval, epart, npart)
    ccall((:METIS_PartMeshNodal, libmetis), Cint, (Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{real_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}), ne, nn, eptr, eind, vwgt, vsize, nparts, tpwgts, options, objval, epart, npart)
end

function METIS_PartMeshDual(ne, nn, eptr, eind, vwgt, vsize, ncommon, nparts, tpwgts, options, objval, epart, npart)
    ccall((:METIS_PartMeshDual, libmetis), Cint, (Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{real_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}), ne, nn, eptr, eind, vwgt, vsize, ncommon, nparts, tpwgts, options, objval, epart, npart)
end

function METIS_NodeND(nvtxs, xadj, adjncy, vwgt, options, perm, iperm)
    ccall((:METIS_NodeND, libmetis), Cint, (Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}), nvtxs, xadj, adjncy, vwgt, options, perm, iperm)
end

function METIS_Free(ptr)
    ccall((:METIS_Free, libmetis), Cint, (Ptr{Cvoid},), ptr)
end

function METIS_SetDefaultOptions(options)
    ccall((:METIS_SetDefaultOptions, libmetis), Cint, (Ptr{idx_t},), options)
end

function METIS_NodeNDP(nvtxs, xadj, adjncy, vwgt, npes, options, perm, iperm, sizes)
    ccall((:METIS_NodeNDP, libmetis), Cint, (idx_t, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, idx_t, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}), nvtxs, xadj, adjncy, vwgt, npes, options, perm, iperm, sizes)
end

function METIS_ComputeVertexSeparator(nvtxs, xadj, adjncy, vwgt, options, sepsize, part)
    ccall((:METIS_ComputeVertexSeparator, libmetis), Cint, (Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}), nvtxs, xadj, adjncy, vwgt, options, sepsize, part)
end

function METIS_NodeRefine(nvtxs, xadj, vwgt, adjncy, where, hmarker, ubfactor)
    ccall((:METIS_NodeRefine, libmetis), Cint, (idx_t, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, Ptr{idx_t}, real_t), nvtxs, xadj, vwgt, adjncy, where, hmarker, ubfactor)
end

@cenum rstatus_et::Int32 begin
    METIS_OK = 1
    METIS_ERROR_INPUT = -2
    METIS_ERROR_MEMORY = -3
    METIS_ERROR = -4
end

@cenum moptype_et::UInt32 begin
    METIS_OP_PMETIS = 0
    METIS_OP_KMETIS = 1
    METIS_OP_OMETIS = 2
end

@cenum moptions_et::UInt32 begin
    METIS_OPTION_PTYPE = 0
    METIS_OPTION_OBJTYPE = 1
    METIS_OPTION_CTYPE = 2
    METIS_OPTION_IPTYPE = 3
    METIS_OPTION_RTYPE = 4
    METIS_OPTION_DBGLVL = 5
    METIS_OPTION_NITER = 6
    METIS_OPTION_NCUTS = 7
    METIS_OPTION_SEED = 8
    METIS_OPTION_NO2HOP = 9
    METIS_OPTION_MINCONN = 10
    METIS_OPTION_CONTIG = 11
    METIS_OPTION_COMPRESS = 12
    METIS_OPTION_CCORDER = 13
    METIS_OPTION_PFACTOR = 14
    METIS_OPTION_NSEPS = 15
    METIS_OPTION_UFACTOR = 16
    METIS_OPTION_NUMBERING = 17
    METIS_OPTION_HELP = 18
    METIS_OPTION_TPWGTS = 19
    METIS_OPTION_NCOMMON = 20
    METIS_OPTION_NOOUTPUT = 21
    METIS_OPTION_BALANCE = 22
    METIS_OPTION_GTYPE = 23
    METIS_OPTION_UBVEC = 24
end

@cenum mptype_et::UInt32 begin
    METIS_PTYPE_RB = 0
    METIS_PTYPE_KWAY = 1
end

@cenum mgtype_et::UInt32 begin
    METIS_GTYPE_DUAL = 0
    METIS_GTYPE_NODAL = 1
end

@cenum mctype_et::UInt32 begin
    METIS_CTYPE_RM = 0
    METIS_CTYPE_SHEM = 1
end

@cenum miptype_et::UInt32 begin
    METIS_IPTYPE_GROW = 0
    METIS_IPTYPE_RANDOM = 1
    METIS_IPTYPE_EDGE = 2
    METIS_IPTYPE_NODE = 3
    METIS_IPTYPE_METISRB = 4
end

@cenum mrtype_et::UInt32 begin
    METIS_RTYPE_FM = 0
    METIS_RTYPE_GREEDY = 1
    METIS_RTYPE_SEP2SIDED = 2
    METIS_RTYPE_SEP1SIDED = 3
end

@cenum mdbglvl_et::UInt32 begin
    METIS_DBG_INFO = 1
    METIS_DBG_TIME = 2
    METIS_DBG_COARSEN = 4
    METIS_DBG_REFINE = 8
    METIS_DBG_IPART = 16
    METIS_DBG_MOVEINFO = 32
    METIS_DBG_SEPINFO = 64
    METIS_DBG_CONNINFO = 128
    METIS_DBG_CONTIGINFO = 256
    METIS_DBG_MEMORY = 2048
end

@cenum mobjtype_et::UInt32 begin
    METIS_OBJTYPE_CUT = 0
    METIS_OBJTYPE_VOL = 1
    METIS_OBJTYPE_NODE = 2
end

const IDXTYPEWIDTH = 32

const REALTYPEWIDTH = 32

const METIS_VER_MAJOR = 5

const METIS_VER_MINOR = 1

const METIS_VER_SUBMINOR = 0

const METIS_NOPTIONS = 40

