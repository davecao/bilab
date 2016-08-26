
#cython: cdivision=True
#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False

"""
  Python interface for metis
  convert bilab.graph's adjacency list to numeric format

METIS format:
  adjacency list =
          [                                xadj[0] = 0
           [1, 2, 3, 4],        vertex 0 - xadj[1] = 4
           [0],                 vertex 1 - xadj[2] = 5
           [0],                 vertex 2 - xadj[3] = 6
           [0],                 vertex 3 - xadj[4] = 7
           [0, 5],              vertex 4 - xadj[5] = 9
           [4, 6],              vertex 5 - xadj[6] = 11
           [13, 5, 7],          vertex 6 - xadj[7] = 14
           [8, 6],              vertex 7 - xadj[8] = 16
           [9, 10, 11, 12, 7],  vertex 8 - xadj[9] = 21
           [8],                 vertex 9 - xadj[10] = 22
           [8],                 vertex 10 - xadj[11] = 23
           [8],                 vertex 11 - xadj[12] = 24
           [8],                 vertex 12 - xadj[13] = 25
           [14, 6],             vertex 13 - xadj[14] = 27
           [13, 15],            vertex 14 - xadj[15] = 29
           [16, 17, 18, 14],    vertex 15 - xadj[16] = 33
           [15],                vertex 16 - xadj[17] = 34
           [15],                vertex 17 - xadj[18] = 35
           [15]                 vertex 18 - xadj[19] = 36
          ]

For metis:
  xadj = [0,4,5,6,7,9,11,14,16,21,22,23,24,25,27,29,33,34,35,36]
  adjncy = [
    1,2,3,4,
    0,
    0,
    0,
    0,5,
    4,6,
    13,5,7,
    8,6,
    9,10,11,12,7,
    8,
    8,
    8,
    8,
    14,6,
    13,15,
    16,17,18,14,
    15,
    15,
    15]
"""

include "gpmetis_wrap.pyi"


cpdef void PrintInfo(np.ndarray[DTYPE_t, ndim=1] opt):
    print("METIS Options:")
    print("METIS_OPTION_PTYPE({}): {}".format(
          METIS_OPTION_PTYPE, opt[METIS_OPTION_PTYPE]))
    print("METIS_OPTION_OBJTYPE({}): {}".format(
          METIS_OPTION_OBJTYPE, opt[METIS_OPTION_OBJTYPE]))
    print("METIS_OPTION_CTYPE({}): {}".format(
          METIS_OPTION_CTYPE, opt[METIS_OPTION_CTYPE]))
    print("METIS_OPTION_IPTYPE({}): {}".format(
          METIS_OPTION_IPTYPE, opt[METIS_OPTION_IPTYPE]))
    print("METIS_OPTION_RTYPE({}): {}".format(
          METIS_OPTION_RTYPE, opt[METIS_OPTION_RTYPE]))
    print("METIS_OPTION_DBGLVL({}): {}".format(
          METIS_OPTION_DBGLVL, opt[METIS_OPTION_DBGLVL]))
    print("METIS_OPTION_NITER({}): {}".format(
          METIS_OPTION_NITER, opt[METIS_OPTION_NITER]))
    print("METIS_OPTION_NCUTS({}): {}".format(
          METIS_OPTION_NCUTS, opt[METIS_OPTION_NCUTS]))
    print("METIS_OPTION_SEED({}): {}".format(
          METIS_OPTION_SEED, opt[METIS_OPTION_SEED]))
    print("METIS_OPTION_NO2HOP({}): {}".format(
          METIS_OPTION_NO2HOP, opt[METIS_OPTION_NO2HOP]))
    print("METIS_OPTION_MINCONN({}): {}".format(
          METIS_OPTION_MINCONN, opt[METIS_OPTION_MINCONN]))
    print("METIS_OPTION_CONTIG({}): {}".format(
          METIS_OPTION_CONTIG, opt[METIS_OPTION_CONTIG]))
    print("METIS_OPTION_COMPRESS({}): {}".format(
          METIS_OPTION_COMPRESS, opt[METIS_OPTION_COMPRESS]))
    print("METIS_OPTION_CCORDER({}): {}".format(
          METIS_OPTION_CCORDER, opt[METIS_OPTION_CCORDER]))
    print("METIS_OPTION_PFACTOR({}): {}".format(
          METIS_OPTION_PFACTOR, opt[METIS_OPTION_PFACTOR]))
    print("METIS_OPTION_NSEPS({}): {}".format(
          METIS_OPTION_NSEPS, opt[METIS_OPTION_NSEPS]))
    print("METIS_OPTION_UFACTOR({}): {}".format(
          METIS_OPTION_UFACTOR, opt[METIS_OPTION_UFACTOR]))
    print("METIS_OPTION_NUMBERING({}): {}".format(
          METIS_OPTION_NUMBERING, opt[METIS_OPTION_NUMBERING]))
    print("METIS_OPTION_HELP({}): {}".format(
          METIS_OPTION_HELP, opt[METIS_OPTION_HELP]))
    print("METIS_OPTION_TPWGTS({}): {}".format(
          METIS_OPTION_TPWGTS, opt[METIS_OPTION_TPWGTS]))
    print("METIS_OPTION_NCOMMON({}): {}".format(
          METIS_OPTION_NCOMMON, opt[METIS_OPTION_NCOMMON]))
    print("METIS_OPTION_NOOUTPUT({}): {}".format(
          METIS_OPTION_NOOUTPUT, opt[METIS_OPTION_NOOUTPUT]))
    print("METIS_OPTION_BALANCE({}): {}".format(
          METIS_OPTION_BALANCE, opt[METIS_OPTION_BALANCE]))
    print("METIS_OPTION_GTYPE({}): {}".format(
          METIS_OPTION_GTYPE, opt[METIS_OPTION_GTYPE]))
    print("METIS_OPTION_UBVEC({}): {}".format(
          METIS_OPTION_UBVEC, opt[METIS_OPTION_UBVEC]))

    print("METIS_PTYPE_RB: {}".format(METIS_PTYPE_RB))
    print("METIS_PTYPE_KWAY: {}".format(METIS_PTYPE_KWAY))
    print("METIS_IPTYPE_GROW: {}".format(METIS_IPTYPE_GROW))
    print("METIS_IPTYPE_RANDOM: {}".format(METIS_IPTYPE_RANDOM))
    print("METIS_RTYPE_GREEDY: {}".format(METIS_RTYPE_GREEDY))


def gpmetis(G, nparts=2, has_vwgt=False, has_edwgt=False,
            ptype="kway", ctype="shem", rtype="fm",
            no2hop=False, dbg_mode=False):

    if not isinstance(G, Graph):
        raise ValueError("Argument g is not bilab.graph.Graph")

    # ------ declare C parameters -------

    # The number of vertices in the graph
    cdef idx_t nvtxs

    # The number of edges in the graph
    cdef idx_t nedges

    # The number of balancing constraints. It should be at least 1
    cdef idx_t ncon

    # The adjacency structure of the graph, vertices's index
    # cdef np.ndarray[DTYPE_t, ndim = 1] xadj
    cdef idx_t *xadj

    # adjacency structure of the graph, index of vertices
    # cdef np.ndarray[DTYPE_t, ndim = 1] adjncy
    cdef idx_t *adjncy

    # The weights of the vertices (NULL)
    # cdef np.ndarray[DTYPE_t, ndim = 1] vwgt
    cdef idx_t *vwgt

    # edge weight: (NULL)
    # cdef np.ndarray[DTYPE_t, ndim = 1] adjwgt
    cdef idx_t *adjwgt

    # The size of the vertices for computing the total communication volumne
    # (NULL)
    # cdef np.ndarray[DTYPE_t, ndim = 1] vsize
    cdef idx_t *vsize

    # An array of size: nparts * ncon (NULL)
    # This is an array of size nparts×ncon that specifies the desired weight
    # for each partition and constraint. The target partition weight for the
    # ith partition and jth constraint is specified at tpwgts[i*ncon+j] (the
    # numbering for both partitions and constraints starts from 0). For each
    # constraint, the sum of the tpwgts[] entries must be 1.0 (i.e.,
    # \sum_{i} tpwgts[i ∗ ncon + j] = 1.0).
    # A NULL value can be passed to indicate that the graph should be equally
    # divided among the partitions.
    # cdef np.ndarray[DTYPE_d, ndim = 1] tpwgts
    cdef real_t *tpwgts

    # This is an array of size ncon that specifies the allowed load imbalance
    # tolerance for each constraint. For the ith partition and jth constraint
    # the allowed weight is the ubvec[j]*tpwgts[i*ncon+j] fraction of the jth’s
    # constraint total weight. The load imbalances must be greater than 1.0.
    # A NULL value can be passed indicating that the load imbalance tolerance
    # for each constraint should be 1.001 (for ncon=1) or 1.01 (for ncon 1).
    # cdef np.ndarray[DTYPE_d, ndim = 1] ubvec
    cdef real_t *ubvec

    # This is the array of options as described in Section 5.4.
    # METIS_PartGraphRecursive - The following options are valid
    #    METIS_OPTION_CTYPE, METIS_OPTION_IPTYPE, METIS_OPTION_RTYPE,
    #    METIS_OPTION_NO2HOP, METIS_OPTION_NCUTS, METIS_OPTION_NITER,
    #    METIS_OPTION_SEED, METIS_OPTION_UFACTOR, METIS_OPTION_NUMBERING,
    #    METIS_OPTION_DBGLVL
    #
    # METIS PartGraphKway
    #    METIS_OPTION_OBJTYPE, METIS_OPTION_CTYPE, METIS_OPTION_IPTYPE,
    #    METIS_OPTION_RTYPE, METIS_OPTION_NO2HOP, METIS_OPTION_NCUTS,
    #    METIS_OPTION_NITER, METIS_OPTION_UFACTOR, METIS_OPTION_MINCONN,
    #    METIS_OPTION_CONTIG, METIS_OPTION_SEED, METIS_OPTION_NUMBERING,
    #    METIS_OPTION_DBGLVL
    cdef np.ndarray[DTYPE_t, ndim = 1] options

    # Upon successful completion, this variable stores the edge-cut or the
    # total communication volume of the partitioning solution. The value
    # returned depends on the partitioning’s objective function.
    cdef idx_t objval = 0

    # This is a vector of size nvtxs that upon successful completion stores
    # the partition vector of the graph. The numbering of this vector starts
    # from either 0 or 1, depending on the value of
    # options[METIS OPTION NUMBERING].
    # cdef np.ndarray[DTYPE_t, ndim = 1] part
    cdef idx_t *part
    cdef idx_t[:] part_view

    cdef idx_t nparts_c

    # Function returns: local variables
    cdef int status
    cdef int i
    xadj_inx = {}

    # initialization
    nvtxs = G.get_num_vertices()
    ncon = 1
    nedges = G.get_num_edges()

    # xadj = np.zeros(nvtxs + 1, dtype=DTYPE)
    # adjncy = np.zeros(nedges * 2, dtype=DTYPE)
    # vwgt = np.ones(nvtxs, dtype=DTYPE)
    # # vwgt_c = &vwgt[0]  # np.ones(nvtxs + 1, dtype=DTYPE)
    # vsize = np.ones(nvtxs, dtype=DTYPE)
    # # vsize_c = &vsize[0]  # np.ones(nvtxs + 1, dtype=DTYPE)
    # adjwgt = np.ones(nedges * 2, dtype=DTYPE)
    # # adjwgt_c = &adjwgt[0]  # np.zeros(2 * ncon, dtype=DTYPE)
    # #target partition weights.
    # #If no specified the weights are set to 1/nparts
    # tpwgts = 1.0 / nparts * np.ones(nparts * ncon, dtype=np.float64)
    # # tpwgts_c = &tpwgts[0]  # np.zeros(nparts * ncon, dtype=np.float64)
    # if ncon == 1:
    #   ubvec = 1.001 * np.ones(ncon + 1, dtype=np.float64)
    # # ubvec_c = NULL   # np.zeros(ncon, dtype=np.float64)
    options = np.ones(_METIS_NOPTIONS, dtype=DTYPE)

    # part = np.zeros(nvtxs, dtype=DTYPE)
    # Allocate memory
    # xadj = <idx_t *>malloc((nvtxs + 1) * sizeof(idx_t))
    # adjncy = <idx_t *>malloc((nedges * 2) * sizeof(idx_t))
    # vwgt = <idx_t *>malloc((ncon * nvtxs) * sizeof(idx_t))
    # adjwgt = <idx_t *>malloc((nedges * 2) * sizeof(idx_t))
    # vsize = <idx_t *>malloc(nvtxs * sizeof(idx_t))
    # tpwgts = <real_t *>malloc((nparts * ncon) * sizeof(real_t))
    ubvec = NULL

    xadj = ismalloc(nvtxs + 1, 0, "ReadGraph: xadj")
    adjncy = imalloc(nedges * 2, "ReadGraph: adjncy")
    vwgt = ismalloc(nvtxs, 1, "ReadGraph: vwgt")
    adjwgt = ismalloc(nedges * 2, 1, "ReadGraph: adjwgt")
    vsize = ismalloc(nvtxs, 1, "ReadGraph: vsize")
    tpwgts = rsmalloc(nparts * ncon, -1.0, "ReadTpwgts: tpwgts")

    # options = <idx_t *>malloc(_METIS_NOPTIONS * sizeof(idx_t))
    part = <idx_t *>malloc(nvtxs * sizeof(idx_t))
    part_view = <idx_t[:nvtxs]>part
    nparts_c = nparts

    # set METIS default options, -1
    METIS_SetDefaultOptions(&options[0])

    for i in range(nparts):
        for j in range(ncon):
            tpwgts[i * ncon + j] = 1.0 / nparts
    # prepare metis input
    i = 1
    # 1. python code to generate index for each vertex
    for v in G:
        xadj_inx[v.name] = (i, v.get_num_neighbors())
        i += 1

    # n_inx: v.name
    # val: (number, num_of_neighbors) a tuple
    sorted_inx = sorted([(val, k) for (k, val) in xadj_inx.items()])

    neighbor_counter = 0
    for val, n_inx in sorted_inx:
        frm = G[n_inx]
        for n in frm.get_connections():
            if G.get_edge(frm, n):
                # if dbg_mode:
                #  print("{}".format(xadj_inx[n.name][0]))
                # Metis start from zero
                adjncy[neighbor_counter] = xadj_inx[n.name][0] - 1
                neighbor_counter += 1
            else:
                print("Edge: {} - {} not found.".format(frm.name, n.name))
        xadj[val[0]] = neighbor_counter
        # if dbg_mode:
        #    print("{} - {}".format(val[0], neighbor_counter))

    # debug mode
    if dbg_mode:
        print("% #vertices #edges:")
        print(" {} {}".format(nvtxs, nedges))
        print("% Here begins the list of vertex neighbors for each vertex")
        for i in range(1, nvtxs + 1):
            start = xadj[i - 1]
            num = xadj[i] - xadj[i - 1]
            # s = " {}".format(i)
            s = " "
            for j in range(start, start + num):
                s += " {}  ".format(adjncy[j] + 1)
            print("{}".format(s.rstrip()))

    # set default options
    # optype_enum.METIS_OPTION_PTYPE      Specifies the partitioning method.
    #      METIS_PTYPE_RB(0)       Multilevel recursive bisectioning.
    #      METIS_PTYPE_KWAY(1)     Multilevel k-way partitioning.
    if ptype == "kway":
        options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY
    elif ptype == "rb":
        options[METIS_OPTION_PTYPE] = METIS_PTYPE_RB

    # optype_enum.METIS_OPTION_OBJTYPE  Specifies the type of objective.
    #            METIS_OBJTYPE_CUT(0) Edge-cut minimization.
    #            METIS_OBJTYPE_VOL(1) Total communication volume minimization.
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT

    # optype_enum.METIS_OPTION_CTYPE   Specifies the matching scheme to be used
    #           METIS_CTYPE_RM(0)      Random matching.
    #           METIS_CTYPE_SHEM(1)    Sorted heavy-edge matching.
    if ctype == "shem":
        options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM
    else:
        options[METIS_OPTION_CTYPE] = METIS_CTYPE_RM

    # optype_enum.METIS_OPTION_IPTYPE
    #   Determines the algorithm used during initial partitioning.
    # iptype_enum.METIS_IPTYPE_GROW(0)
    #     Grows a bisection using a greedy strategy.
    # iptype_enum.METIS_IPTYPE_RANDOM(1)
    #     Computes a bisection at random followed by a refinement.
    # iptype_enum.METIS_IPTYPE_EDGE(2)
    #     Derives a separator from an edge cut.
    # iptype_enum.METIS_IPTYPE_NODE(3)
    #     Grows a bisection using a greedy node-based strategy.
    if options[METIS_OPTION_IPTYPE] == -1:
        if options[METIS_OPTION_PTYPE] == METIS_PTYPE_RB:
            if ncon == 1:
                options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_GROW
            else:
                options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_RANDOM
        else:
            options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_METISRB

    # optype_enum.METIS_OPTION_RTYPE
    #       Determines the algorithm used for refinement.
    # rtype_enum.METIS_RTYPE_FM(0)    FM-based cut refinement.
    # rtype_enum.METIS_RTYPE_GREEDY(1)  Greedy-based cut and volume refinement.
    # rtype_enum.METIS_RTYPE_SEP2SIDED(2) Two-sided node FM refinement.
    # rtype_enum.METIS_RTYPE_SEP1SIDED(3) One-sided node FM refinement.
    options[METIS_OPTION_RTYPE] = METIS_RTYPE_FM

    # optype_enum.METIS_OPTION_NO2HOP Specifies that the coarsening will not
    #                             perform any 2-hop matchings when the standard
    #                             matching approach fails to sufficiently
    #                             coarsen the graph. The 2-hop matching is very
    #                             effective for graphs with power-law degree
    #                             distributions.
    #         0                   Performs a 2-hop matching.
    #         1                   Does not perform a 2-hop matching.
    if no2hop:
        options[METIS_OPTION_NO2HOP] = 1
    else:
        options[METIS_OPTION_NO2HOP] = 0

    # optype_enum.METIS_OPTION_NCUTS  Specifies the number of different
    #                             partitionings that it will compute. The final
    #                             partitioning is the one that achieves the
    #                             best edgecut or communication volume. Default
    #                             is 1.
    options[METIS_OPTION_NCUTS] = 1

    # optype_enum.METIS_OPTION_NITER
    #     Specifies the number of iterations for the
    #     refinement algorithms at each stage of the
    #     uncoarsening process. Default is 10.
    options[METIS_OPTION_NITER] = 10

    # optype_enum.METIS_OPTION_UFACTOR
    #  Specifies the maximum allowed load imbalance
    #  among the partitions. A value of x indicates
    #  that the allowed load imbalance is
    #  (1+x)/1000. The load imbalance for the jth
    #  constraint is defined to be
    #  max_i(w(j,i)/t(j,i)), where w(j,i) is the
    #  fraction of the overall weight of the jth
    #  constraint that is assigned to the ith
    #  partition and t(j,i) is the desired target
    #  weight of the jth constraint for the ith
    #  partition (i.e., that specifies via -tpwgts).
    #  For recursive bisectioning, the default value
    #  is 1 (i.e., load imbalance of 1.001) and for
    #  k-way partitioning, the default value is 30
    #  (i.e., load imbalance of 1.03).
    if options[METIS_OPTION_PTYPE] == METIS_PTYPE_RB:
        options[METIS_OPTION_UFACTOR] = 1
    else:
        options[METIS_OPTION_UFACTOR] = 30

    # optype_enum.METIS_OPTION_MINCONN
    #       Specifies that the partitioning routines
    #       should try to minimize the maximum degree of
    #       the subdomain graph, i.e., the graph in
    #       which each partition is a node, and edges
    #       connect subdomains with a shared interface.
    if options[METIS_OPTION_PTYPE] == METIS_PTYPE_KWAY:
        options[METIS_OPTION_MINCONN] = 1
    else:
        options[METIS_OPTION_MINCONN] = 0

    # optype_enum.METIS_OPTION_CONTIG  Specifies that the partitioning routines
    #              should try to produce partitions that are
    #              contiguous. Note that if the input graph is
    #              not connected this option is ignored.
    #          0                       Does not force contiguous partitions.
    #          1                       Forces contiguous partitions.
    options[METIS_OPTION_CONTIG] = 0

    # optype_enum.METIS_OPTION_SEED
    #  Specifies the seed for the random number generator.
    options[METIS_OPTION_SEED] = -1

    # optype_enum.METIS_OPTION_NUMBERING
    # options_c[METIS_OPTION_NUMBERING] = -1

    # optype_enum.METIS_OPTION_DBGLVL
    #                           Specifies the amount of progress/debugging
    #                           information will be printed during the
    #                           execution of the algorithms. The default
    #                           value is 0 (no debugging/progress
    #                           information). A non-zero value can be
    #                           supplied that is obtained by a bit-wise OR of
    #                           the following values.
    #   METIS_DBG_INFO(1)       Prints various diagnostic messages.
    #   METIS_DBG_TIME(2)       Performs timing analysis.
    #   METIS_DBG_COARSEN(4)    Displays various statistics during
    #                           coarsening.
    #   METIS_DBG_REFINE(8)     Displays various statistics during
    #                           refinement.
    #   METIS_DBG_IPART(16)     Displays various statistics during initial
    #                           partitioning.
    #   METIS_DBG_MOVEINFO(32)  Displays detailed information about vertex
    #                           moves during refinement.
    #   METIS_DBG_SEPINFO(64)   Displays information about vertex separators.
    #   METIS_DBG_CONNINFO(128) Displays information related to the
    #                           minimization of subdomain connectivity.
    #   METIS_DBG_CONTIGINFO(256)
    #                           Displays information related to the
    #                                  elimination of connected components.
    if dbg_mode:
        options[METIS_OPTION_DBGLVL] = METIS_DBG_INFO
    else:
        options[METIS_OPTION_DBGLVL] = 0

    # if dbg_mode:
    #     PrintInfo(options)

    gk_malloc_init()
    if options[METIS_OPTION_PTYPE] == METIS_PTYPE_KWAY:
        #status = METIS_PartGraphKway(
        #              &nvtxs, &ncon,
        #              &xadj[0], &adjncy[0],
        #              &vwgt[0], &vsize[0], &adjwgt[0],
        #              &nparts,
        #              &tpwgts[0], &ubvec[0],
        #              &options[0],
        #              &objval,
        #              &part[0])
        status = METIS_PartGraphKway(
                      &nvtxs, &ncon,
                      xadj, adjncy,
                      vwgt, vsize, adjwgt,
                      &nparts_c,
                      tpwgts, ubvec,
                      &options[0],
                      &objval,
                      part)
    elif options[METIS_OPTION_PTYPE] == METIS_PTYPE_RB:
        # status = METIS_PartGraphRecursive(
        #           &nvtxs, &ncon,
        #           &xadj[0], &adjncy[0],
        #           &vwgt[0], &vsize[0], &adjwgt[0],
        #           &nparts,
        #           &tpwgts[0], &ubvec[0],
        #           &options[0],
        #           &objval,
        #           &part[0])
        status = METIS_PartGraphRecursive(
                  &nvtxs, &ncon,
                  xadj, adjncy,
                  vwgt, vsize, adjwgt,
                  &nparts_c,
                  tpwgts, ubvec,
                  &options[0],
                  &objval,
                  part)

    if not status == METIS_OK:
        print("*Metis returned with an error")
        if status == METIS_ERROR_INPUT:
            print("The input error occurred")
        elif status == METIS_ERROR_MEMORY:
            print("Could not allocate the required memory")
        elif status == METIS_ERROR:
            print("Some error")

    # prepare return for groups

    # release memory
    free(xadj)
    free(adjncy)
    free(vwgt)
    free(adjwgt)
    free(vsize)
    free(tpwgts)
    free(part)
    # for val, n_inx in sorted_inx:
    #     print("{}:{}".format(n_inx, part_view[val[0] - 1]))
    gk_malloc_cleanup(0)

    return sorted_inx, objval, list(part_view)
