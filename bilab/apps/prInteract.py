#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
prInteract is an application used to find the interacted atoms pairs in a pdb file
"""

from __future__ import print_function
import sys
import os
# import errno
# import re

from datetime import datetime as dt
from xml.etree.ElementTree import Element, SubElement
# from xml.etree.ElementTree import Comment, tostring
from optparse import OptionParser, make_option

__execute_name__ = "prInteract"
__author__ = "Wei Cao"
__contact__ = "davecao@bi.a.u-tokyo.ac.jp"
__date__ = "2014/08/04"
__version__ = "0.1"
__copyright__ = """
    Feel free to do whatever you like with this code.
    Keep this script relative to bilab package
    """

"""
versionInstallation of bilab package:

  git clone https://github.com/davecao/bilab.git

Usage:

./prInteract.py --pdb 9mht.pdb --source protein --target nucleic -v
--target :  not water and hetero
            hetero and (resname 753 or resname DMS or resname GOL)

Note: if --pdbdir is not specified, it will use the current directory
      as resource.
      target option: for example
      protein:
          'ALA', 'ARG', 'ASN', 'ASP', 'ASX', 'CSO', 'CYS', 'GLN', 'GLU',
          'GLX', 'GLY', 'HIP', 'HIS', 'HSD', 'HSE', 'HSP', 'ILE', 'LEU',
          'LYS', 'MET', 'MSE', 'PHE', 'PRO', 'PTR', 'SEC', 'SEP', 'SER',
          'THR', 'TPO', 'TRP', 'TYR', 'VAL', 'XAA', 'XLE'
      nucleic:
          nucleobase: 'GUN', 'ADE', 'CYT', 'THY', 'URA',
          nucleotide: 'DA', 'DC', 'DG', 'DT', 'DU', 'A', 'C', 'G', 'T', 'U'
          nucleoside: 'AMP', 'ADP', 'ATP', 'CDP', 'CTP', 'GMP', 'GDP', 'GTP'
                      'TMP', 'TTP', 'UMP', 'UDP', 'UTP'
      ion: 'AL', 'BA', 'CA', 'CD', 'CL', 'CO', 'CS', 'CU', 'CU1', 'CUA',
           'HG', 'IN', 'IOD', 'K', 'MG', 'MN3', 'NA', 'PB', 'PT', 'RB',
           'TB', 'TL', 'WO4', 'YB', 'ZN'
      lipid: GPE, LPP, OLA, SDS, and STE
            (POPC, LPPC, POPE, DLPE, PCGL, STEA, PALM, OLEO, DMPC
             from CHARMM force field).
      sugar:'GLC', 'GLO', 'BGC'
"""


class PDBInfo(object):
    """
    Store the header info of the input pdb file
    """

    def __init__(self, title, pdbid,
                 experiment_method,
                 dep_date=None,
                 classification=None,
                 resolution=None,
                 version=None):
        """Initialization"""
        self.title = title
        self.pdbid = pdbid
        self.exp_method = experiment_method
        self.deposition_date = dep_date
        self.classification = classification
        self.resolution = resolution
        self.version = version

    def __str__(self):
        return "#{}{}{}#{}{}{}#{}{}{}#{}{}{}#{}{}{}#{}{}".format(
            "PDBID: ", self.pdbid, os.linesep,
            "Experiment method: ", self.exp_method, os.linesep,
            "Deposition date: ", self.deposition_date, os.linesep,
            "Classification: ", self.classification, os.linesep,
            "Resolution: ", self.resolution, os.linesep,
            "Version: ", self.version)


class ATOMInfo(object):
    """ Atom info """
    def __init__(self, name, chId,
                 resName,
                 resNum,
                 insertCode="_",
                 atType="residue"):
        self.at_name = name
        self.chain_id = chId
        self.res_name = resName
        self.res_num = resNum
        # residue or ligand
        self.at_type = atType
        self.insert_code = insertCode

    def __str__(self):
        str_des = None
        if self.at_type == "residue":
            str_des = "{}.{}.{}.{}.{}".format(
                self.res_name,
                self.chain_id,
                self.insert_code,
                self.res_num,
                self.at_name
            )
        elif self.at_type == "ligand":
            str_des = "{}.{}.{}.{}".format(
                self.res_name,
                self.chain_id,
                self.res_num,
                self.at_name
            )
        return str_des

    def to_xml(self, node):
        """
        write to xml format
        """
        node.set('type', self.at_type)
        node.set('resName', self.res_name)
        node.set('chId', self.chain_id)
        node.set('resNum', "{}".format(self.res_num))
        node.set('atName', self.at_name)
        if self.at_type == "residue":
            node.set('InsertCode', self.insert_code)


class AtomIteraction(object):
    """ store a pair of interacted atoms """
    def __init__(self, at1, at2, dist, cov_dist):
        self.atom1 = at1
        self.atom2 = at2
        self.dist = dist
        self.cov_dist = cov_dist
        self.cov_bond = True if cov_dist >= dist else False

    def __str__(self):
        return "{}:{}:{}:{}:{}".format(
            self.atom1,
            self.atom2,
            self.dist,
            self.cov_dist,
            self.cov_bond)

    def to_xml(self, parent):
        """ write to xml """
        # parent is an SubElement object of xml
        parent.set('distance', "{}".format(self.dist))
        parent.set('covalent_dist', "{}".format(self.cov_dist))
        parent.set('cov_bond', "{}".format(self.cov_bond))
        atom1_elem = SubElement(parent, 'ATOM')
        self.atom1.to_xml(atom1_elem)
        atom2_elem = SubElement(parent, 'ATOM')
        self.atom2.to_xml(atom2_elem)


class Interactions(PDBInfo):
    """ Store the interaction pairs at the atomic level """

    def __init__(self, *args, **kwargs):
        """Initialization"""
        self.dist_threshold = kwargs.pop('dist_threshold', None)
        self.create_date = dt.now()
        # call parent class constructor
        super(Interactions, self).__init__(*args, **kwargs)
        self.interact_list = []

    def __str__(self):
        return super(Interactions, self).__str__() + \
               "#Create date: {}".format(self.create_date)

    def set_interaction_pair(self, neighbours):
        """ convert neighbours to list """
        append_self = self.interact_list.append
        for at1, at2, distance in neighbours:
            res1_icode = "_"
            if at1.getIcode():
                res1_icode = at1.getIcode()
            atom1 = ATOMInfo(
                at1.getName(),
                at1.getChid(),
                at1.getResname(),
                at1.getResnum(),
                insertCode=res1_icode,
                atType="residue"
            )
            atom2 = ATOMInfo(
                at2.getName(),
                at2.getChid(),
                at2.getResname(),
                at2.getResnum(),
                atType="ligand"
            )
            element1 = at1.getElement()
            element2 = at2.getElement()
            cov_bond_sum = bilab.chemicals.get_covalent_radius(element1) + \
                bilab.chemicals.get_covalent_radius(element2)
            atom_pair = AtomIteraction(atom1, atom2, distance, cov_bond_sum)
            append_self(atom_pair)

    def write(self,
              outfmt='txt',
              ofile='out'):
        """ neighbours """
        ofs = None
        if isinstance(ofile, bilab.string_types):
            ofs = open(ofile, 'w')
        if outfmt == 'txt':
            self.to_txt(ofile=ofs)
        elif outfmt == 'xml':
            self.to_xml(ofile=ofs)
        else:
            print("Unkown output format.")
            sys.exit(1)
        ofs.close()

    def to_txt(self, ofile='out'):
        """ write the interaction pairs to a text file"""
        h_info = "#---- Generated by {} v{} ----------------".format(
            __execute_name__, __version__)
        print(h_info, end='\n', file=ofile)
        print("{}".format(self), end='\n', file=ofile)
        print("# Note: Column separated by :", end='\n', file=ofile)
        print("# Col1. label started with ATOM", end='\n', file=ofile)
        print("# Col2. residue info. the format is below.", end='\n', file=ofile)
        print("#       resName.chId.InsertCode.at_name - ", end='\n', file=ofile)
        print("#       e.g., ASP.A._.29.CA ", end='\n', file=ofile)
        print("# Col3. ligand info. the format is below.", end='\n', file=ofile)
        print("#       resName.chId.InsertCode.at_name - ", end='\n', file=ofile)
        print("#       e.g., 017.B.201.N1", end='\n', file=ofile)
        print("# Col4. the distance between the two atoms in Col2 and Col3.", end='\n', file=ofile)
        print("# Col5. summation of covalent radii of the two atoms.", end='\n', file=ofile)
        print("# Col6. Covalent bond or not. True/False", end='\n', file=ofile)
        print("#"+"-"*(len(h_info) - 1), end='\n', file=ofile)
        for at_pair in self.interact_list:
            print("ATOM:{}".format(at_pair), end='\n', file=ofile)

    def to_xml(self, ofile='out'):
        """ Write to xml """
        # Generate root element
        root = Element(
            'PDB',
            {'id': self.pdbid,
             'method': self.exp_method,
             'deposition_date': self.deposition_date,
             'classification': self.classification,
             'resolution': "{}".format(self.resolution),
             'version': "{}".format(self.version),
             'dist_threshold': "{}".format(self.dist_threshold),
             'create_date': self.create_date.strftime('%Y-%m-%d %H:%M:%S')
            })
        child = SubElement(root, 'Interactions')
        for at_pair in self.interact_list:
            interaction = SubElement(child, 'Interaction')
            at_pair.to_xml(interaction)
        ofile.write(bilab.utilities.prettify_xml(root))


def parse_cmd(argv):
    """
    Parse command line arguments
    """
    option_list = [
        make_option("--pdb", dest="pdbfile",
                    help="The input pdb file[REQUIRED]."),
        make_option("--source", dest="source",
                    help="The name of a source molecule in the pdb. "
                    "option: protein, nucleic [REQUIRED]."),
        make_option("--target", dest="target",
                    help="The name of target molecule in the pdb."
                    "protein, nucleic or ligand [REQUIRED]."),
        make_option("--distance", dest="distance", action='store',
                    type="float",
                    default=5.0,
                    help="The distance threshold for selecting "
                    "interacting pair [OPTION]. Default is 5.0 angstroms."),
        make_option("--outfmt", dest='outfmt', default='txt',
                    help="output format: xml or txt. default is txt."),
        make_option("--out", dest='outfilename', default='out',
                    help="output file name. If not specified, the name will"
                    " be composed of out.fmt"),
        make_option("-v", "--verbose",
                    action="store_true", dest="verbose", default=False,
                    help="print verbose info"),
    ]

    usage = 'usage: %prog [options] --pdb PDBFILE' \
            '--source protein --target protein/dna/rna'
    parser = OptionParser(formatter=IndentedHelpFormatterWithNL(),
                          option_list=option_list,
                          usage=usage,
                          version=__version__)
    options, arguments = parser.parse_args(argv)


    if not options.pdbfile:
        print("Error: do not specifiy the input pdb file")
        parser.print_help()
        sys.exit(1)

    if not options.source:
        print("Error: do not specifiy source name")
        parser.print_help()
        sys.exit(1)

    if not options.target:
        print("Error: do not specifiy target name")
        parser.print_help()
        sys.exit(1)
    return options


def select_atoms(mol, pdbid, select_cmd):
    """ Select atoms by select_cmd from mol """
    selected_atoms = mol.select(select_cmd)
    if selected_atoms is None:
        print("Error: No atoms found for {}".format(pdbid))
        sys.exit(1)
    return selected_atoms


def main(argv):
    """
    Main entry point
    """
    # parse command line arguments
    opt = parse_cmd(argv)
    # load the pdb file
    mol, header = bilab.structure.parsePDB(opt.pdbfile, header=True)
    pdbid = header['identifier']
    # store pdb info
    inter_atoms_set = Interactions (
        header['title'],
        pdbid,
        header['experiment'],
        dep_date=header['deposition_date'],
        classification=header['classification'],
        resolution=header['resolution'] if 'resolution' in header else None,
        version=header['version'],
        dist_threshold=opt.distance
    )
    # set title of the molecule
    mol.setTitle(pdbid)

    # select source
    source_atoms = select_atoms(mol, pdbid, opt.source)
    target_atoms = select_atoms(mol, pdbid, opt.target)

    if opt.verbose:
        print("Atoms in pdb file: {0}".format(mol.numAtoms()))
        print("Atoms of the source: {0}".format(source_atoms.numAtoms()))
        print("Atoms of the target: {0}".format(target_atoms.numAtoms()))

    # find the interacting pairs within the specified distance
    neighbors = bilab.structure.findNeighbors(
        source_atoms, opt.distance, target_atoms)

    inter_atoms_set.set_interaction_pair(neighbors)
    if opt.outfilename == "out":
        opt.outfilename += "." + opt.outfmt
    inter_atoms_set.write(outfmt=opt.outfmt, ofile=opt.outfilename)


if __name__ == '__main__':
    CURR_VERSION = sys.version_info
    if not CURR_VERSION[:2] == (2, 7):
        print("Python v.{}.{} is not supported".format(
            CURR_VERSION[0], CURR_VERSION[1]))
        sys.exit(0)
    # import bilab package
    try:
        import bilab
        from bilab.utilities import IndentedHelpFormatterWithNL
    except ImportError:
        raise ImportError("This script needs bilab package")
    main(sys.argv)
