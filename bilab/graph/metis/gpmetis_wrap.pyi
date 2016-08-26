
#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

"""
  interface for metis C header
"""

import cython
from cython.view cimport array
# import array
from libc.stdint cimport uintptr_t, int64_t
from libc.stdlib cimport malloc, free
from bilab.graph import Graph, Vertex

from cpython cimport PyObject, Py_INCREF

# Import the Python-level symbols of numpy
import numpy as np

# Import the C-level symbols of numpy
cimport numpy as np

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()

DTYPE = np.int
ctypedef np.int_t DTYPE_t
ctypedef np.float64_t DTYPE_d
cdef int METIS_OPTION_PTYPE
ctypedef int64_t idx_t
ctypedef double real_t

# GKlib/gk_struct.h
#cdef extern from "gk_struct.h" nogil:
#  #*************************************************************************
#  #*! The following data structure stores information about a memory 
#  #  allocation operation that can either be served from gk_mcore_t or by
#  #  a gk_malloc if not sufficient workspace memory is available.
#  #*************************************************************************
#  cdef struct gk_mop_t:
#    int type
#    ssize_t nbytes
#    void *ptr
#
#  cdef struct gk_mcore_t:
#    #Workspace information
#    size_t coresize  #!< The amount of core memory that has been allocated
#    size_t corecpos  #*!< Index of the first free location in core
#    void *core       #*!< Pointer to the core itself
#
#    # These are for implementing a stack-based allocation scheme using both
#    # core and also dynamically allocated memory
#    size_t nmops  #*!< The number of maop_t entries that have been allocated
#    size_t cmop   #*!< Index of the first free location in maops
#    gk_mop_t *mops  #*!< The array recording the maop_t operations
#
#    # These are for keeping various statistics for wspacemalloc
#    size_t num_callocs  #*!< The number of core mallocs
#    size_t num_hallocs  #*!< The number of heap mallocs
#    size_t size_callocs #*!< The total # of bytes in core mallocs
#    size_t size_hallocs #*!< The total # of bytes in heap mallocs
#    size_t cur_callocs  #*!< The current # of bytes in core mallocs
#    size_t cur_hallocs  #*!< The current # of bytes in heap mallocs
#    size_t max_callocs  #*!< The maximum # of bytes in core mallocs at any given time
#    size_t max_hallocs  #*!< The maximum # of bytes in heap mallocs at any given time

# libmetis/struct.h
cdef extern from "libmetis/struct.h":

  cdef struct ckrinfo_t:
      idx_t id     #!< The internal degree of a vertex (sum of weights)
      idx_t ed     #!< The total external degree of a vertex
      idx_t nnbrs  #!< The number of neighboring subdomains
      idx_t inbr   #!< The index in the cnbr_t array where the nnbrs list 
                        #  of neighbors is stored

  cdef struct vkrinfo_t:
      idx_t nid  #!< The internal degree of a vertex (count of edges)
      idx_t ned  #!< The total external degree of a vertex (count of edges)
      idx_t gv   #!< The volume gain of moving that vertex
      idx_t nnbrs  #!< The number of neighboring subdomains
      idx_t inbr  #!< The index in the vnbr_t array where the nnbrs list
                       # of neighbors is stored

  cdef struct nrinfo_t:
      idx_t edegrees[2]

  cdef struct graph_t:
      idx_t nvtxs
      idx_t nedges
      idx_t *xadj
      idx_t *vwgt
      idx_t *vsize   # Vertex sizes for min-volume formulation
      idx_t *adjncy  # Array that stores the adjacency lists of nvtxs
      idx_t *adjwgt  # Array that stores the weights of the adjacency lists

      idx_t *tvwgt   # The sum of the vertex weights in the graph
      real_t *invtvwgt # The inverse of the sum of the vertex weights in the graph

      # These are to keep track control if the corresponding fields correspond to
      # application or library memory
      int free_xadj
      int free_vwgt
      int free_vsize
      int free_adjncy
      int free_adjwgt

      idx_t *label

      idx_t *cmap

      # Partition parameters
      idx_t mincut
      idx_t minvol
      idx_t *where
      idx_t *pwgts
      idx_t nbnd
      idx_t *bndptr
      idx_t *bndind
   
      # Bisection refinement parameters
      idx_t *id
      idx_t *ed
   
      # K-way refinement parameters
      ckrinfo_t *ckrinfo   #!< The per-vertex cut-based refinement info
      vkrinfo_t *vkrinfo   #!< The per-vertex volume-based refinement info
   
      # Node refinement information
      nrinfo_t *nrinfo

      graph_t *coarser
      graph_t *finer

# bilab/graph/metis/libmetis/metis.h
cdef extern from "metis.h":
  # Return codes
  ctypedef enum rstatus_et:
    METIS_OK = 1
    METIS_ERROR_INPUT = -2
    METIS_ERROR_MEMORY = -3
    METIS_ERROR = -4

  # ! Operation type codes
  ctypedef enum moptype_et:
    METIS_OP_PMETIS
    METIS_OP_KMETIS
    METIS_OP_OMETIS

  # define default Metis options 
  ctypedef enum moptions_et:
    METIS_OPTION_PTYPE,
    METIS_OPTION_OBJTYPE,
    METIS_OPTION_CTYPE,
    METIS_OPTION_IPTYPE,
    METIS_OPTION_RTYPE,
    METIS_OPTION_DBGLVL,
    METIS_OPTION_NITER,
    METIS_OPTION_NCUTS,
    METIS_OPTION_SEED,
    METIS_OPTION_NO2HOP,
    METIS_OPTION_MINCONN,
    METIS_OPTION_CONTIG,
    METIS_OPTION_COMPRESS,
    METIS_OPTION_CCORDER,
    METIS_OPTION_PFACTOR,
    METIS_OPTION_NSEPS,
    METIS_OPTION_UFACTOR,
    METIS_OPTION_NUMBERING,

    # Used for command-line parameter purposes
    METIS_OPTION_HELP,
    METIS_OPTION_TPWGTS,
    METIS_OPTION_NCOMMON,
    METIS_OPTION_NOOUTPUT,
    METIS_OPTION_BALANCE,
    METIS_OPTION_GTYPE,
    METIS_OPTION_UBVEC

  # *! Partition Schemes
  ctypedef enum mptype_et:
    METIS_PTYPE_RB,
    METIS_PTYPE_KWAY

  # *! Graph type for meshes
  ctypedef enum mgtype_et:
    METIS_GTYPE_DUAL,
    METIS_GTYPE_NODAL

  # *! Coarsening Schemes
  ctypedef enum mctype_et:
    METIS_CTYPE_RM,
    METIS_CTYPE_SHEM

  # *! Initial partitioning schemes
  ctypedef enum miptype_et:
    METIS_IPTYPE_GROW,
    METIS_IPTYPE_RANDOM,
    METIS_IPTYPE_EDGE,
    METIS_IPTYPE_NODE,
    METIS_IPTYPE_METISRB

  # *! Refinement schemes
  ctypedef enum mrtype_et:
    METIS_RTYPE_FM,
    METIS_RTYPE_GREEDY,
    METIS_RTYPE_SEP2SIDED,
    METIS_RTYPE_SEP1SIDED

  # *! Debug levels
  ctypedef enum mdbglvl_et:
    METIS_DBG_INFO = 1, # *!< Shows various diagnostic messages
    METIS_DBG_TIME = 2, # *!< Perform timing analysis
    METIS_DBG_COARSEN = 4, # *!< Show the coarsening progress
    METIS_DBG_REFINE = 8, # *!< Show the refinement progress
    METIS_DBG_IPART = 16, # *!< Show info on initial partitioning
    METIS_DBG_MOVEINFO = 32, # *!< Show info on vertex moves during refinement
    METIS_DBG_SEPINFO = 64, # *!< Show info on vertex moves during sep refinement
    METIS_DBG_CONNINFO = 128, # *!< Show info on minimization of subdomain connectivity
    METIS_DBG_CONTIGINFO = 256, # *!< Show info on elimination of connected components
    METIS_DBG_MEMORY = 2048 # *!< Show info related to wspace allocation

  # *! Types of objectives
  ctypedef enum mobjtype_et:
    METIS_OBJTYPE_CUT,
    METIS_OBJTYPE_VOL,
    METIS_OBJTYPE_NODE

# bilab/graph/metis/metisbin.h
cdef extern from "metisbin.h":
    int METIS_VER_MAJOR
    int METIS_VER_MINOR
    int METIS_VER_SUBMINOR

    int IDXTYPEWIDTH
    int REALTYPEWIDTH

    int METIS_NOPTIONS

    # # implemented in libmetis/pmetis.c
    # int METIS_PartGraphRecursive(
    #     idx_t *nvtxs,   # number of vertices in the graph
    #     idx_t *ncon,    # number of edges in the graph
    #     idx_t *xadj,    # adjacency structure of the graph. size: nvtxs + 1
    #     idx_t *adjncy,  # adjacency structure of the graph. size: 2*nedges
    #     idx_t *vwgt,    # (NULL) the weights of the vertices.
    #     idx_t *vsize,   # (NULL) the size of the vertices for computing the total communication volume
    #     idx_t *adjwgt,  # (NULL) the weights of the edges.
    #     idx_t *nparts,  # The number of parts to partition the graph
    #     real_t *tpwgts, # An array of size nparts x ncon
    #     real_t *ubvec,  #(NULL) An array of size ncon
    #     idx_t *options, #(NULL) This is the array of options
    #     idx_t *objval, #
    #     idx_t *part     # The is a vector of size nvtxs that upon successful
    #                     # completion stores the partition vector of the graph.
    #                     # The numbering of this vector starts from either 0 or 1,
    #                     # depending on the value of options[METIS_OPTION_NUMBERING]
    #     )
  
    # implemented in libmetis/kmetis.c
    # For k-way clustering, the appropriate options are::
    #    objtype   = 'cut' or 'vol'
    #    ctype     = 'rm' or 'shem'
    #    iptype    = 'grow', 'random', 'edge', 'node'
    #    rtype     = 'fm', 'greedy', 'sep2sided', 'sep1sided'
    #    ncuts     = integer, number of cut attempts (default = 1)
    #    niter     = integer, number of iterations (default = 10)
    #    ufactor   = integer, maximum load imbalance of (1+x)/1000
    #    minconn   = bool, minimize degree of subdomain graph
    #    contig    = bool, force contiguous partitions
    #    seed      = integer, RNG seed
    #    numbering = 0 (C-style) or 1 (Fortran-style) indices
    #    dbglvl    = Debug flag bitfield
    # int METIS_PartGraphKway(
    #     idx_t  *nvtxs,
    #     idx_t  *ncon,
    #     idx_t  *xadj,
    #     idx_t  *adjncy,
    #     idx_t  *vwgt,
    #     idx_t  *vsize,
    #     idx_t  *adjwgt,
    #     idx_t  *nparts,
    #     real_t *tpwgts,
    #     real_t *ubvec,
    #     idx_t  *options,
    #     idx_t  *objval,
    #     idx_t  *part)

cdef extern from "auxapi.c":
  int METIS_SetDefaultOptions(idx_t *options)
  int METIS_Free(void *ptr)

cdef extern from "pmetis.c":
  # implemented in libmetis/pmetis.c
  int METIS_PartGraphRecursive(
      idx_t *nvtxs,   # number of vertices in the graph
      idx_t *ncon,    # number of balancing constraints. should be at least 1.
      idx_t *xadj,    # adjacency structure of the graph. size: nvtxs + 1
      idx_t *adjncy,  # adjacency structure of the graph. size: 2*nedges
      idx_t *vwgt,    # (NULL) the wights of the vertices.
      idx_t *vsize,   # (NULL) the size of the vertices for computing the total communication volume
      idx_t *adjwgt,  # (NULL) the weights of the edges.
      idx_t *nparts,  # The number of parts to partition the graph
      real_t *tpwgts, # An array of size nparts x ncon
      real_t *ubvec,  #(NULL) An array of size ncon
      idx_t *options, #(NULL) This is the array of options
      idx_t *objval, #
      idx_t *part     # The is a vector of size nvtxs that upon successful
                      # completion stores the partition vector of the graph.
                      # The numbering of this vector starts from either 0 or 1,
                      # depending on the value of options[METIS_OPTION_NUMBERING]
      )

cdef extern from "kmetis.c":
  # implemented in libmetis/kmetis.c
  int METIS_PartGraphKway(
      idx_t  *nvtxs,
      idx_t  *ncon,
      idx_t  *xadj,
      idx_t  *adjncy,
      idx_t  *vwgt,
      idx_t  *vsize,
      idx_t  *adjwgt,
      idx_t  *nparts,
      real_t *tpwgts,
      real_t *ubvec,
      idx_t  *options,
      idx_t  *objval,
      idx_t  *part)

cdef extern from "memory.c":
  int gk_malloc_init()
  void gk_malloc_cleanup(int showstats)
  void *gk_malloc(size_t nbytes, char *msg)


cdef extern from "gklib_rename.h":
    idx_t* ismalloc(int n, int i, char* msg)
    real_t* rsmalloc(int n, double i, char* msg)
    idx_t* imalloc(int i, char* msg)

# export to python
_METIS_VER_MAJOR = METIS_VER_MAJOR
_METIS_VER_MINOR = METIS_VER_MINOR
_METIS_VER_SUBMINOR = METIS_VER_SUBMINOR
_IDXTYPEWIDTH = IDXTYPEWIDTH
_REALTYPEWIDTH = REALTYPEWIDTH
_METIS_NOPTIONS = METIS_NOPTIONS

