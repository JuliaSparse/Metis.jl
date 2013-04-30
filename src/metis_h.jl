const METIS_NOPTIONS = 40
## Return codes
const METIS_OK           = int32(1)  # normal return
const METIS_ERROR_INPUT  = int32(-2) # erroneous inputs and/or options
const METIS_ERROR_MEMORY = int32(-3) # insufficient memory 
const METIS_ERROR        = int32(-4) # Other errors
## Operation type codes
const METIS_OP_PMETIS = int32(0)
const METIS_OP_KMETIS = int32(1)
const METIS_OP_OMETIS = int32(2)
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
const METIS_PTYPE_RB   = int32(0)
const METIS_PTYPE_KWAY = int32(1)
## Graph types for meshes
const METIS_GTYPE_DUAL  = int32(0)
const METIS_GTYPE_NODAL = int32(1)            
## Coarsening Schemes
const METIS_CTYPE_RM   = int32(0)
const METIS_CTYPE_SHEM = int32(1)
## Initial partitioning schemes
const METIS_IPTYPE_GROW    = int32(0)
const METIS_IPTYPE_RANDOM  = int32(1)
const METIS_IPTYPE_EDGE    = int32(2)
const METIS_IPTYPE_NODE    = int32(3)
const METIS_IPTYPE_METISRB = int32(4)
## Refinement schemes
const METIS_RTYPE_FM        = int32(0)
const METIS_RTYPE_GREEDY    = int32(1)
const METIS_RTYPE_SEP2SIDED = int32(2)
const METIS_RTYPE_SEP1SIDED = int32(3)
## Debug levels (bit positions)
const METIS_DBG_INFO       = int32(1) # Shows various diagnostic messages
const METIS_DBG_TIME       = int32(2) # Perform timing analysis
const METIS_DBG_COARSEN    = int32(4) # Show the coarsening progress
const METIS_DBG_REFINE     = int32(8) # Show the refinement progress
const METIS_DBG_IPART      = int32(16) # Show info on initial partitioning
const METIS_DBG_MOVEINFO   = int32(32) # Show info on vertex moves during refinement
const METIS_DBG_SEPINFO    = int32(64) # Show info on vertex moves during sep refinement
const METIS_DBG_CONNINFO   = int32(128) # Show info on minimization of subdomain connectivity
const METIS_DBG_CONTIGINFO = int32(256) # Show info on elimination of connected components
const METIS_DBG_MEMORY     = int32(2048) # Show info related to wspace allocation
## Types of objectives
const METIS_OBJTYPE_CUT  = int32(0)
const METIS_OBJTYPE_VOL  = int32(1)
const METIS_OBJTYPE_NODE = int32(2)
