const idx_t = Cint
const real_t = Cfloat

const METIS_NOPTIONS = 40

## Return codes
const METIS_OK           = Cint(1)  # Returned normally
const METIS_ERROR_INPUT  = Cint(-2) # Returned due to erroneous inputs and/or options
const METIS_ERROR_MEMORY = Cint(-3) # Returned due to insufficient memory
const METIS_ERROR        = Cint(-4) # Some other errors

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
