#!/usr/bin python
# -*- coding: utf-8 -*-
"""
    This app is based on tree comparison of ete3 package.
    Please refer to ete3_compare in ete3's tools 
"""
__author__ = "Wei Cao"
__contact__ = "davecao@bi.a.u-tokyo.ac.jp"
__date__ = "2015/08/30"
__version__ = "0.1"
__copyright__ = """
    Feel free to do whatever you like with this code.
    """
import os
import sys
import uuid
import time
import numpy as np
from signal import signal, SIGPIPE, SIG_DFL
from optparse import OptionParser, make_option

try:
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import *
    matplotlib_version = matplotlib.__version__
except ImportError:
    raise ImportError('Matplotlib is a required package')

try:
    import ete3
    from ete3 import Tree, TreeStyle
except ImportError:
    raise ImportError('ete3 is a required package')

#import bilab
from bilab.utilities.logger import Console

# reset
signal(SIGPIPE, SIG_DFL)
# global
LOGGER = None

def parse_cmd(argv):
    """
        Parse command line arguments
    """
    option_list = [

        make_option("--source", dest="source",
                  help="The location of a source tree in newwick [REQUIRED]."),

        make_option("--ref", dest="reference",
                  help="The location of a reference tree in newwick [REQUIRED]."),

        make_option("--method", dest="method", default="RF",
                    help="The method for tree comparison: "\
                       "  RF: Robinson-Foulds. Only supported now."),

        make_option("--unrooted", dest="unrooted", 
                    action = "store_true", default=True,
                    help="The input trees are unrooted"),

        make_option("--treefmt", dest="treefmt", default=1,type="int",
                    help="""The newick format of input trees.
                            1 - flexible with internal node names.
                            7 - leaf branches + all names         
                            8 - all names                         
                            9 - leaf names                        
                            100 - topology only"""),
#        make_option("--treeko", dest="treeko", 
#                    action = "store_true", default=False,
#                    help="The input trees are unrooted"),

        make_option("--show_mismatches", dest="show_mismatches", 
                    action = "store_true", default=False,
                    help="show mismatches between the input trees."),

        make_option("--show_matches", dest="show_matches",
                    action = "store_true", default=False,
                    help="show matches between the input trees."),

        make_option("--show_edges", dest="show_edges",
                    action = "store_true", default=False,
                    help="show matches between the input trees."),

        make_option("--outfmt", dest="outfmt", default='txt',
                    help="output format: txt or tabulate. default is txt."),

        make_option("--out", dest="out", default=None,
                    help="the name of output file."),

        make_option("--log", dest='logfilename', default=None,
                    help="The name of a log file. "+ \
                    "If not specified, the name will" + \
                    " be composed of pdbid.log"),

        make_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=False,
                  help="print verbose info")
    ]

    usage = 'usage: %prog [options] --source tree1.nw' \
            ' --ref tree.nw --unrooted --show_matches' \
            ' --show_mismatches --outfmt txt --log tree1_vs_tree2.log' \
            ' -v'

    parser = OptionParser(option_list=option_list,
                           usage=usage,
                           version=__version__)
    options, arguments = parser.parse_args(argv)
    outfmt_list = ['txt', 'tabulate']
    if arguments == 0:
        print ("Error: no arguments found")
        parser.print_help()
        sys.exit(1)

    if not options.source:
        print ("Error: do not specify tree1")
        parser.print_help()
        sys.exit(1)

    if not options.reference:
        print ("Error: do not specify tree2")
        parser.print_help()
        sys.exit(1)
    
    if options.outfmt not in outfmt_list:
        print ("Error: unknown specified output format.")
        parser.print_help()
        sys.exit(1)

    if options.out is None: 
        print ("Error: No output file specified.")
        parser.print_help()
        sys.exit(1)

    return options

def loadTree(tree, format=1):
    """ 
        Load a tree in newick format

    Args:
        tree (string): the file name of a tree stored in newick format
    
    Return:
        A tree object generated by ete3 pacakge
    
    About format
    .. table::
          ======  ==============================================
          FORMAT  DESCRIPTION
          ======  ==============================================
          0        flexible with support values
          1        flexible with internal node names
          2        all branches + leaf names + internal supports
          3        all branches + all names
          4        leaf branches + leaf names
          5        internal and leaf branches + leaf names
          6        internal branches + leaf names
          7        leaf branches + all names
          8        all names
          9        leaf names
          100      topology only
          ======  ==============================================
    """

    with open (tree, "r") as treefile:
        tStr=treefile.read().replace('\n', '')
    t = Tree(tStr, format=format)
    return t

def print_(src, ref, src_name, ref_name, unrooted=True):
    if unrooted:
        for tag, part in [("src: %s"%src_name, src), ("ref: %s"%ref_name, ref)]:
            print("%s\t%s" %(tag, 
                '\t'.join(
                map(lambda x: '%s|%s' %(','.join(x[0]), ','.join(x[1])), part))))
    else:
        for tag, part in [("src: %s"%src_name, src), ("ref: %s"%ref_name, ref)]:
            print("%s\t%s" %(tag, '\t'.join([','.join(p) for p in part])))

def main(argv):
    global LOGGER
    # parse command line arguments
    opt = parse_cmd(argv)
    # input
    source = opt.source
    reference = opt.reference
    method = opt.method
    
    #treeko = opt.treeko
    # input bool
    unrooted = opt.unrooted
    treefmt = opt.treefmt
    #output
    outfmt = opt.outfmt
    outputfile = opt.out
    # output bool
    show_matches = opt.show_matches
    show_mismatches = opt.show_mismatches
    show_edges = opt.show_edges

    #log
    if opt.logfilename is None:
        logfile = 'treeCompare.log'
    else:
        logfile = opt.logfilename

    LOGGER = Console("tree", prefix='treeCompare>')
    LOGGER.start(logfile)

    header = ['source', 'ref', 'eff.size', 'nRF',
              'RF', 'maxRF', "%src_branches",
              "%ref_branches", "subtrees", "treekoD" ]


    # Assign the commmand line arguments
    t_src = loadTree(source, format=treefmt)
    t_ref = loadTree(reference, format=treefmt)

    #compare(): 
    #  ref_tree, 
    #  use_collateral=False, 
    #  min_support_source=0.0, min_support_ref=0.0,
    #  has_duplications=False, 
    #  expand_polytomies=False, 
    #  unrooted=False,
    #  max_treeko_splits_to_be_artifact=1000, 
    #  ref_tree_attr='name', 
    #  source_tree_attr='name'
    r = t_src.compare(t_ref, unrooted=unrooted, has_duplications=False)
    if show_mismatches or show_matches or show_edges:
        if show_mismatches:
            src = r['source_edges'] - r['ref_edges']
            ref = r['ref_edges'] - r['source_edges']
            print_(src, ref, source, reference, unrooted=unrooted)
        if show_matches:
            src = r['source_edges'] & r['ref_edges']
            ref = r['ref_edges'] & r['source_edges']
            print_(src, ref, source, reference, unrooted=unrooted)
        if show_edges:
            src = r['source_edges']
            ref = r['ref_edges']
            print_(src, ref, source, reference, unrooted=unrooted)
        
    else:
        data = [source,
                reference,
                r['effective_tree_size'],
                r['norm_rf'],
                r['rf'], r['max_rf'],
                r["source_edges_in_ref"],
                r["ref_edges_in_source"],
                r['source_subtrees'],
                r['treeko_dist']]
        #if r['effective_tree_size'] == 0:
        #    for i in xrange(len(data)):
        #        data[i] = -1
        if outfmt == 'tabulate':
            print('# ' + '\t'.join(header))
            print('\t'.join(map(str, data)))
        elif outfmt == 'txt':
            print('# ' + ' '.join(header))
            print(' '.join(map(str,data)))

    #with open(outputfile, 'w') as fh:
    #    fh.write('# ' + '\t'.join(header))


if __name__ == '__main__':
    main(sys.argv)