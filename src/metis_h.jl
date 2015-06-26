const METIS_NOPTIONS = 40
## Return codes
const METIS_OK           = @compat Int32(1)  # normal return
const METIS_ERROR_INPUT  = @compat Int32(-2) # erroneous inputs and/or options
const METIS_ERROR_MEMORY = @compat Int32(-3) # insufficient memory
const METIS_ERROR        = @compat Int32(-4) # Other errors
## Operation type codes
const METIS_OP_PMETIS = @compat Int32(0)
const METIS_OP_KMETIS = @compat Int32(1)
const METIS_OP_OMETIS = @compat Int32(2)
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
const METIS_PTYPE_RB   = @compat Int32(0)
const METIS_PTYPE_KWAY = @compat Int32(1)
## Graph types for meshes
const METIS_GTYPE_DUAL  = @compat Int32(0)
const METIS_GTYPE_NODAL = @compat Int32(1)
## Coarsening Schemes
const METIS_CTYPE_RM   = @compat Int32(0)
const METIS_CTYPE_SHEM = @compat Int32(1)
## Initial partitioning schemes
const METIS_IPTYPE_GROW    = @compat Int32(0)
const METIS_IPTYPE_RANDOM  = @compat Int32(1)
const METIS_IPTYPE_EDGE    = @compat Int32(2)
const METIS_IPTYPE_NODE    = @compat Int32(3)
const METIS_IPTYPE_METISRB = @compat Int32(4)
## Refinement schemes
const METIS_RTYPE_FM        = @compat Int32(0)
const METIS_RTYPE_GREEDY    = @compat Int32(1)
const METIS_RTYPE_SEP2SIDED = @compat Int32(2)
const METIS_RTYPE_SEP1SIDED = @compat Int32(3)
## Debug levels (bit positions)
const METIS_DBG_INFO       = @compat Int32(1) # Shows various diagnostic messages
const METIS_DBG_TIME       = @compat Int32(2) # Perform timing analysis
const METIS_DBG_COARSEN    = @compat Int32(4) # Show the coarsening progress
const METIS_DBG_REFINE     = @compat Int32(8) # Show the refinement progress
const METIS_DBG_IPART      = @compat Int32(16) # Show info on initial partitioning
const METIS_DBG_MOVEINFO   = @compat Int32(32) # Show info on vertex moves during refinement
const METIS_DBG_SEPINFO    = @compat Int32(64) # Show info on vertex moves during sep refinement
const METIS_DBG_CONNINFO   = @compat Int32(128) # Show info on minimization of subdomain connectivity
const METIS_DBG_CONTIGINFO = @compat Int32(256) # Show info on elimination of connected components
const METIS_DBG_MEMORY     = @compat Int32(2048) # Show info related to wspace allocation
## Types of objectives
const METIS_OBJTYPE_CUT  = @compat Int32(0)
const METIS_OBJTYPE_VOL  = @compat Int32(1)
const METIS_OBJTYPE_NODE = @compat Int32(2)
