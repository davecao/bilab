# -*- coding: utf-8 -*-
__all__ = []

try:
    from bilab.graph.metis import _gpmetis_wrap
except ImportError as exc:
    print("Error: failed to import module {}".format(exc))
    raise ImportError('Could not load _gpmetis_wrap.so')
else:
    pass

METIS_VER_MAJOR = _gpmetis_wrap._METIS_VER_MAJOR
METIS_VER_MINOR = _gpmetis_wrap._METIS_VER_MINOR
METIS_VER_SUBMINOR = _gpmetis_wrap._METIS_VER_SUBMINOR

IDXTYPEWIDTH = _gpmetis_wrap._IDXTYPEWIDTH
REALTYPEWIDTH = _gpmetis_wrap._REALTYPEWIDTH

__version__ = "METIS v{}.{}.{}".format(METIS_VER_MAJOR, METIS_VER_MINOR,
                                       METIS_VER_SUBMINOR)


def gp_partition(g, nparts=2, has_vwgt=False, has_edwgt=False,
                 ptype="kway", ctype="shem", rtype="fm", optNiter=10,
                 no2hop=False, dbg_mode=False):
    """ Interface to METIS_PartGraphRecursive and METIS_PartGraphKway
        in METIS package.

    Args:
        g (object):  object of bilab.graph.Graph

    Kwargs:
        nparts (int):
        ptype (str): 'kway', 'rb'.
        ctype (str): 'shem', 'rm'.
        rtype (str): 'fm', 'greedy', 'sep2sided', 'sep1sided'.
        optNiter (int): iteration steps for refinement, default is 10.
        no2hop (boolean): True, False.
        dbg_mode (boolean): True, False.
    """
    name_inx, objval, groups = _gpmetis_wrap.gpmetis(
        g, nparts=nparts,
        has_vwgt=has_vwgt, has_edwgt=has_edwgt,
        ptype=ptype, ctype=ctype, rtype=rtype, optNiter=optNiter,
        no2hop=no2hop, dbg_mode=dbg_mode)

    part_result = {}
    for inx, n in name_inx:
        part_result[n] = groups[inx[0] - 1]
    #    print("{}: {}".format(n, part_result[n]))
    return part_result
