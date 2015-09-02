#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Wei Cao"
__contact__ = "davecao@bi.a.u-tokyo.ac.jp"
__date__ = "2015/08/27"
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

import bilab
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
        make_option("--pdbfile", dest="pdbfile",
                  help="The location of the pdb files [REQUIRED]."),

        make_option("--pdbid", dest="pdbid",
                  help="The pdbid in four letters [REQUIRED]."),

        make_option("--output-dir", dest="out_dir", default=".",
                  help="The output directory for storing the "+\
                       "intermidiate results. Default is '.'."),

        make_option("--mapping-sampling", dest="map_sampling", default=100,
                  type="int",
                  help="The sampling times when using NDR , default is 100."),

        make_option("--sel", dest="selection", default='protein and name CA',
                  help="Atom selection, default is 'protein and name CA'."),

        make_option("--draw-density", action="store_true", 
                    dest='draw_density', default=False,
                    help="draw a density map if True. default is False."),

        make_option("--save-density", dest='save_density', default=None,
                    help="save a density map to a txt file." + \
                         " Use pdbid as file name if not specified." + \
                         " e.g. 9mht_density.txt"),

        make_option("--save-mapping", dest='save_mapping', default=None,
                    help="save mapping data in text format." + \
                         " Use pdbid as file name if not specified." + \
                         " e.g. 9mht_mapping.txt"),

        make_option("--build-tree", dest='build_tree', default="random",
                    help="Specify the strategy how to select a node \
                     when building up a tree if True. 'random' or 'binary'. \
                     default is 'random'."),

        make_option("--draw-tree", action="store_true",
                    dest='draw_tree', default=False,
                    help="draw a tree if True. default is False."),

        make_option("--save-tree", dest='save_tree', default=None,
                    help="save a tree in newick format." +\
                         " Use pdbid as file name if not specified." + \
                         " e.g. 9mht_tree.nw"),

        make_option("--save-tree-fmt", dest='save_tree_fmt', default=1,
                    type="int",
                    help="save a tree in newick format.\n" +\
                         "   1 : flexible with internal node names.\n" + \
                         " 100 : topology only"),
#        make_option("--outfmt", dest='outfmt', default='txt',
#                    help="output format: xml or txt. default is txt." + \
#                         "Currently text output is only support."),

        make_option("--log", dest='logfilename', default='unamed',
                    help="The name of a log file. "+ \
                    "If not specified, the name will" + \
                    " be composed of pdbid.log"),

        make_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=False,
                  help="print verbose info"),
    ]

    usage = 'usage: %prog [options] --pdbfile 9mht.pdb' \
            '--pdbid 9mht --draw-density ' \
            '--draw-tree --log 9mht.log'

    parser = OptionParser(option_list=option_list,
                           usage=usage,
                           version=__version__)
    options, arguments = parser.parse_args(argv)
    
    nodeSelection = ['random', 'median']

    if arguments == 0:
        print ("Error: no arguments found")
        parser.print_help()
        sys.exit(1)

    if not options.pdbid:
        print ("Error: do not specify pdbid")
        parser.print_help()
        sys.exit(1)

    if not options.pdbfile:
        print ("Error: do not specify a pdb file")
        parser.print_help()
        sys.exit(1)

    if not options.build_tree in nodeSelection:
        print ("Error: Unknown specified strategy for \
               building a vp tree - {}".format(options.build_tree))
        parser.print_help()
        sys.exit(1)
#    if options.save_density is None:
#        options.save_density = options.pdbid + "_density.txt"
        #print("Error: do not specify the output, --save_density.")
        #parser.print_help()
        #sys.exit(1)

#    if not options.save_tree:
#        options.save_tree = options.pdbid + "_tree.nw"
        #print("Error: do not specify the output, --save_tree.")
        #parser.print_help()
        #sys.exit(1)

#    if not options.save_mapping:
#        options.save_density = options.pdbid + "_density.txt"
        #print("Error: do not specify the output, --save_mapping.")
        #parser.print_help()
        #sys.exit(1)

    #print "Left options %s" % arguments
    return options

def sanitize_id(id):
    return id.strip().replace(" ", "")

class Point(object):

    def __init__(self, identifier=None, name=None, value=0.0):
        super(Point, self).__init__()
        self.__identifier = (str(uuid.uuid1()) if identifier is None else sanitize_id(str(identifier)))
        self.name = name
        self.value = value

    def __str__(self):
        return "{}, name={}, value={}".format(self.__identifier, self.name, self.value)

    def __lt__(self, other):
        return self.value < other.value

    def __sub__(self, other):
        return self.value - other.value


def NDRmapping(X, pdbid, chain, save_mapping_file, no_dims=1, iters=100, 
                  perplexity=30.0, verbose=False):
    """
        Mapping the density matrix to 1d by t-SNE
    Args:
        X (ndarray): the density matrix (N x N)
    Kwargs:
        no_dims (integer) : the final dimension of mapping X
        sampling (integer) : iteration times of performing t-SNE
        perplexity (float) : the perplexity
        verbose (bool) : show the verbosity
    """
    if not isinstance(X, (np.ndarray, np.generic)):
        raise TypeError("X: Inappropriate argument type for {}"
                    .format(X.__class__.__name__))

    N = X.shape[0]
    #dseq = list(chain.getSequence())

    if X.ndim == no_dims:
        print("The dimension of the input and output is equal. No need to do further.")
        return
    
    if (N-1) < (3.0 * perplexity):
        perplexity = (N-1) * 0.3
    
    if verbose:
        start_time = time.time()
    
    Y = bilab.ml.NDR.tSNE.bhtsne(X, no_dims=1, sampling=iters, 
                                perplexity=perplexity, 
                                verbose=verbose)
    if verbose:
        end_time = time.time()
        elapsed_time = "{:.3f}".format(end_time-start_time)
        print("NDRmapping (sampling,{}) : elapsed time:{}".format(iters, \
                                                              elapsed_time))
    results = []
    if save_mapping_file:
        data = []
        dY = [ "{:.6f}".format(i) for i in Y]
        for i, residue in enumerate(chain.iterResidues()):
            res="{}_{}_{}".format(residue.getResname(),residue.getChid(),residue.getResnum())
            results.append(Point(name=res, value=Y[i]))
            resInOne = bilab.utilities.utilities.three2one(residue.getResname())
            data.append([pdbid, residue.getChid(), residue.getResnum(), 
                        residue.getResname(), resInOne, dY[i]])
        header_comment="Mapping ({},{}) - No_dims:{}, Sampling:{}, Perplexity:{}\n".\
                format(N, N, no_dims, iters, perplexity)
        np.savetxt(save_mapping_file, data, fmt="%4s %1s %6s %3s %1s %15s", header=header_comment)
    else:
        for i, residue in enumerate(chain.iterResidues()):
            res="{}_{}_{}".format(residue.getResname(),residue.getChid(),residue.getResnum())
            results.append(Point(name=res, value=Y[i]))
    return results

def saveMat2png(*args, **kwargs):
    """ 
        plot data and save to a png file.
    Args:
        *args : the number of data(list)
    kwargs:
        nrows, ncols : the subplots parameters
        colmap       : the colorbar
        filename     : the name of output graph file
    """
    # get the file name
    fName = kwargs.pop('fName', 'unnamed')
    # the layout of plot
    nrows = kwargs.pop('nrows', 1)
    num_dpi = kwargs.pop('dpi', 300)
    ncols = len(args)
    # colorbar
    colmap = kwargs.pop('colmap', cm.Oranges)
    # create a graph
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols)
    fig.subplots_adjust(wspace=0.22,hspace=0.2)
    if ncols == 1:
        cax = axs.matshow(args[0], cmap=colmap, interpolation='none')
        fig.colorbar(cax, ax=axs,fraction=0.046, pad=0.04)
    else:
        for ax_num in xrange(ncols):
            cax = axs[ax_num].matshow(args[ax_num], cmap=colmap, interpolation='none')
            fig.colorbar(cax, ax=axs[ax_num],fraction=0.046, pad=0.04)
    fig.savefig(fName, dpi=num_dpi) # save to file
    plt.close(fig) # clos the figure

def getDensity(pdbid, ch, save_density_file, drawGraph=False, verbose=False):
    """ 
        Build a distance matrix and normalize elements to [0, 1] 
    
    Args:
        pdbid (string): the pdb id 
        ch (string): chain id of the pdb
        save_density_file (string) : the output file name

    Kwargs:
        outputDir (string) : the output directory if save_density_file is None,
                             default is ".". 
        drawGraph (bool) : plot and save the distance matrix in png if True
        verbose   (bool) : show the verbosity if True
    """
    ch_id = ch.getChid()
    sequence = ch.getSequence()
    distMat = bilab.structure.measure.buildDistMatrix(ch)
    distMat = np.asarray(distMat)

    #normalize distMat
    n_distMat = 1.0*distMat/distMat.max()

    if save_density_file:
        filename = save_density_file
        # save results
        comm_header = "{} {}\n".format(pdbid, ch_id)
        comm_header +="{}\n".format(sequence)
        np.savetxt(filename, n_distMat , delimiter=',', fmt="%8.3f", 
                    header=comm_header)
    if verbose:
        print("Save the density matrix to {}".format(filename))

    # draw a graph
    if drawGraph:
        ext = filename.split('.')[-1]
        graphFileName = filename.replace("."+ext, '.png')
        saveMat2png(distMat, n_distMat, fName=graphFileName)
        if verbose:
            print("draw the density matrix to {}".format(graphFileName))

    return n_distMat

def main(argv):

    global LOGGER
    # parse command line arguments
    opt = parse_cmd(argv)

    # Assign the commmand line arguments
    pdbid = opt.pdbid
    pdb_f = opt.pdbfile
    outputDir = opt.out_dir
    
    sampling_times = opt.map_sampling
    
    selection = opt.selection
    drawDensity = opt.draw_density
    builtTreeStrategy = opt.build_tree
    drawTree = opt.draw_tree
    save_tree_fmt = opt.save_tree_fmt
    verbose = opt.verbose

    if opt.logfilename == 'unamed':
        logfile = pdbid+'.log'
    else:
        logfile = opt.logfilename

    LOGGER = Console(pdbid, prefix=pdbid+'>')
    LOGGER.start(logfile)
    if verbose:
        start_time = time.time()
        LOGGER.timeit(label="Start")
    #Load pdbs
    mol, pdb_header = bilab.structure.parsePDB(pdb_f,title=pdbid, header = True)
    if pdb_header["identifier"] != pdbid:
        print("Input pdbid({}) is inconsistent with the identifier ({}) in the pdb file."
            .format(pdbid, pdb_header["identifier"]))
        print("Here, the function is case sensitive.")
        sys.exit(1)

    # select CA atoms of proteins
    atoms = mol.select(selection)
    # get the number of chains
    hv = atoms.getHierView()
    chains = list(hv)
    for ch in chains:
        ch_id = ch.getChid()
        if opt.save_density is None:
            save_density_file = outputDir + os.sep + pdbid + '_' + \
                            ch_id + '_density.txt'
        else:
            save_density_file = opt.save_density

        if opt.save_mapping is None:
            save_mapping_file = outputDir + os.sep + pdbid + '_' +\
                                ch_id + '_mapping.txt'
        else:
            save_mapping_file = opt.save_mapping

        if opt.save_tree is None:
            save_tree_file = outputDir + os.sep + pdbid + '_' +\
                                ch_id + '_mapping.nw'
        else:
            save_tree_file = opt.save_tree

        # get densityMatrix
        densityMatrix = getDensity(pdbid, ch, save_density_file, 
                                drawGraph=drawDensity, 
                                verbose=verbose)
        map1d = NDRmapping(densityMatrix, pdbid, ch, save_mapping_file, 
                            iters=sampling_times)
        #convert map1d to newick
        t1 = bilab.structure.VPTree(map1d,strategy = builtTreeStrategy)
        t1str = t1.save(format=save_tree_fmt)
        #t.write(format=1, outfile=save_tree_file)
        with open(save_tree_file, "w") as text_file:
            #text_file.write("#{}_{}\n".format(pdbid, ch.getChid()))
            text_file.write("{}".format(t1str))
        if drawTree:
            t = Tree(t1str, format=1)
            ts = TreeStyle()
            ts.show_leaf_name = True
            ts.show_branch_length = True
            ts.show_branch_support = True
            # create a graph
            ext = save_tree_file.split('.')[-1]
            treeGraph = save_tree_file.replace( "."+ext, '.png')
            t.render(treeGraph, w=183, units="mm", dpi=300, tree_style=ts)
    if verbose:
        LOGGER.report(label="Start")
        end_time = time.time()
        elapsed_time = end_time - start_time
        LOGGER.info("Elapsed time: {}".format(elapsed_time))
if __name__ == '__main__':
    main(sys.argv)
