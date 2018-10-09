#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
`prInteract` is an application used to find the interacted atoms pairs in
a pdb file

Installation of bilab package:

  git clone https://github.com/davecao/bilab.git

Usage:

$ python -m bilab.apps.prInteract --pdb 9mht.pdb --source "not water" \
    --target "not water" -v

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

from __future__ import print_function

import sys
import os
import argparse

from datetime import datetime as dt
from bilab import jinja2_ENV

__execute_name__ = "struct_interact"
__author__ = "Wei Cao"
__contact__ = "davecao@bi.a.u-tokyo.ac.jp"
__date__ = "2014/08/04"
__version__ = "0.1"
__copyright__ = """Feel free to do whatever you like with this code.
    Keep this script relative to bilab package
    """


class PDBInfo:
    """
    Store the header info of the input pdb file.

    Args:
        title (str): a title for the experiment or analysis that is
                     represented in the entry.
        pdbid (str): 4-letter id

    Kwargs:
        experiment_method (str):
            X-RAY  DIFFRACTION
            FIBER  DIFFRACTION
            NEUTRON  DIFFRACTION
            ELECTRON  CRYSTALLOGRAPHY
            ELECTRON  MICROSCOPY
            SOLID-STATE  NMR
            SOLUTION  NMR
            SOLUTION  SCATTERING
        classification (str): Classify the molecule(s)
        resolution (float): resolution if available
        version (str):
    """
    def __init__(self, title, pdb_id, experiment_method,
                 dep_date=None,
                 classification=None,
                 resolution=None,
                 version=None, **kwargs):
        """Initialization"""
        self.title = title
        self.pdb_id = pdb_id
        self.exp_method = experiment_method
        self.deposition_date = dep_date
        self.classification = classification
        self.resolution = resolution
        self.version = version

    def __str__(self):
        return "{}".format(self.pdb_id)


class ATOMInfo:
    """ Atom info

    Args:
        at_obj (bilab.structure.atomic.Atom) :
            an object of Atom class in bilab.structure.atomic

    Kwargs:
        insert_code (str): Insert code. Default is "_".
        at_type (str): type of the atom. ["residue", "ligand"]
    """
    def __init__(self, at_obj):
        if not isinstance(at_obj, bilab.structure.atomic.Atom):
            print("The should be an object of {}".format())
            sys.exit(0)

        self.at_name = at_obj.getName()
        self.chain_id = at_obj.getChid()
        self.res_name = at_obj.getResname()
        self.res_num = at_obj.getResnum()

        # residue or heteros
        self.at_type = "protein" if at_obj.isprotein else "hetero"
        if at_obj.getIcode():
            self.insert_code = at_obj.getIcode()

    def __hash__(self):
        return hash(self.__str__())

    def __eq__(self, other):
        if isinstance(other, ATOMInfo):
            return self.__hash__() == other.__hash__()
        return NotImplemented

    def __ne__(self, other):
        return not self.__eq__()

    def __str__(self):
        # str_des = None
        str_des = "{}.{}.{}.{}.{}".format(
            self.res_name,
            self.chain_id,
            self.insert_code,
            self.res_num,
            self.at_name)
        return str_des


class AtomInteraction(object):
    """ store a pair of interacted atoms

    Args:
        atom1 (ATOMInfo): the first atom of an interaction pair
        atom2 (ATOMInfo): the second atom of an interaction pair
        dist (float): distance that is less than a (cut-off) threshold.
        cov_dist (float): length of covalent bond

    """
    def __init__(self, at1, at2, dist, cov_dist):
        self.atom1 = at1
        self.atom2 = at2
        self.dist = dist
        self.cov_dist = cov_dist
        self.interact_type = ""
        if at1.at_type == 'protein' and at2.at_type == 'protein':
            self.interact_type = "P2P"
        elif at1.at_type == 'hetero' and at2.at_type == 'hetero':
            self.interact_type = "H2H"
        else:
            self.interact_type = "H2P"

        if at1.chain_id == at2.chain_id and \
                at1.res_name == at2.res_name and \
                at1.res_num == at2.res_num:
            self.interact_subtype = "intra-residue"
        else:
            self.interact_subtype = "inter-residue"

    def __hash__(self):
        return hash(self.atom1.__hash__() + self.atom2.__hash__())

    def __eq__(self, other):
        if isinstance(other, AtomInteraction):
            return self.__hash__() == other.__hash__()
        return NotImplemented

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        if hasattr(other, 'interact_type'):
            return self.interact_type < other.interact_type
        return NotImplemented

    def __le__(self, other):
        if hasattr(other, 'interact_type'):
            return self.interact_type <= other.interact_type
        return NotImplemented

    def __gt__(self, other):
        if hasattr(other, 'interact_type'):
            return self.interact_type > other.interact_type
        return NotImplemented

    def __ge__(self, other):
        if hasattr(other, 'interact_type'):
            return self.interact_type >= other.interact_type
        return NotImplemented

    def __str__(self):
        return "{}:{}:{:.3f}:{}".format(
            self.atom1,
            self.atom2,
            self.dist,
            self.cov_dist)


class Interactions(PDBInfo):
    """
        Store the interaction pairs at the atomic level

    Attributes:
        dist_threshold (float): the threshold of distance between any two atoms
        create_date (time): date of creation
    """
    def __init__(self, *args, **kwargs):
        # call parent class constructor
        super(Interactions, self).__init__(*args, **kwargs)

        self.writers = {'txt': self.to_txt, 'xml': self.to_xml}
        self.dist_threshold = kwargs.pop('dist_threshold', None)
        self.out_fmt = kwargs.pop('out_fmt', "txt")
        self.create_date = dt.now()
        self.interact_list = []

        self.writer = self.writers[self.out_fmt]

    def __str__(self):
        return super(Interactions, self).__str__() + \
               "#Create date: {}".format(self.create_date)

    def set_interaction_pair(self, neighbours):
        """ convert neighbours to list

        Args:
            neighbours (tuple): two interacted atoms(ATOM) and distance

        Returns:
            None
        """
        append_self = self.interact_list.append
        for at1, at2, distance in neighbours:
            atom1 = ATOMInfo(at1)
            atom2 = ATOMInfo(at2)

            element1 = at1.getElement()
            element2 = at2.getElement()

            # length of covalent bond
            cov_bond_sum = bilab.chemicals.get_covalent_radius(element1) + \
                bilab.chemicals.get_covalent_radius(element2)

            # exclude covalent bonds
            if distance > cov_bond_sum:
                atom_pair = AtomInteraction(atom1,
                                            atom2,
                                            distance,
                                            cov_bond_sum)
                append_self(atom_pair)
        # self.interact_list = self.__unique()
        seen = set()
        seen_add = seen.add
        results = [x for x in self.interact_list
                   if x not in seen and not seen_add(x)]
        self.interact_list = sorted(results)

    def write(self, o_file='out'):
        """ neighbours

        Kwargs:
            out_fmt (str): specify the file format. ['txt', 'xml'],
                           Default is "txt"
            o_file (str): specify the output file name. Default is "out"

        Returns:
            None
        """
        ofs = None
        if isinstance(o_file, bilab.string_types):
            ofs = open(o_file, 'w')
        self.writer(o_file=ofs)
        ofs.close()

    def to_txt(self, o_file='out'):
        """ write the interaction pairs to a text file

        Kwargs:
            o_file (str): output file name
        """
        tpl = jinja2_ENV.get_template("report.txt")
        tags = {
            'prog_name': __execute_name__,
            'prog_ver': __version__,
            'title': self.title,
            'deposition_date': self.deposition_date,
            'method': self.exp_method,
            'pdbid': self.pdb_id,
            'classification': self.classification,
            'resolution': self.resolution,
            'version': self.version,
            'dist_threshold': self.dist_threshold,
            'created_date': self.create_date,
            'interactions': self.interact_list
        }
        o_file.write(tpl.render(tags))

    def to_xml(self, o_file='out'):
        """ Write to xml

        Kwargs:
            o_file (str): output file name.
        """
        tpl = jinja2_ENV.get_template("report.xml")
        tags = {
            'prog_name': __execute_name__,
            'prog_ver': __version__,
            'title': self.title,
            'deposition_date': self.deposition_date,
            'method': self.exp_method,
            'pdbid': self.pdb_id,
            'classification': self.classification,
            'resolution': self.resolution,
            'version': self.version,
            'dist_threshold': self.dist_threshold,
            'created_date': self.create_date,
            'interactions': self.interact_list
        }
        o_file.write(tpl.render(tags))


def get_file_ext(file_path):
    """ Extract the file extension

    Args:
        file_path (str): the file path

    Returns:
        str: the file extension
    """
    filename_w_ext = os.path.basename(file_path)
    filename, file_extension = os.path.splitext(filename_w_ext)
    return file_extension


def select_atoms(mol, pdb_id, select_cmd):
    """ Select atoms by select_cmd from mol

    Args:
        mol ():
        pdb_id (str): identifier used in Protein DataBank (PDB)
        select_cmd (str): selection expression

    Returns:
        selections: return selected atoms by expression
    """
    selected_atoms = mol.select(select_cmd)
    if selected_atoms is None:
        print("Error: No atoms found for {}".format(pdb_id))
        sys.exit(1)
    return selected_atoms


def main(cli_opts):
    """
    Main entry point
    """
    # load the pdb file
    mol, header = cli_opts.fileloader(cli_opts.pdbfile, header=True)
    pdb_id = header['identifier']

    # store pdb info
    inter_atoms_set = Interactions(
        header['title'],
        pdb_id,
        header['experiment'],
        dep_date=header['deposition_date'],
        classification=header['classification'],
        resolution=header['resolution'] if 'resolution' in header else None,
        version=header['version'],
        dist_threshold=cli_opts.distance,
        out_fmt=cli_opts.outfmt
    )

    # set title of the molecule
    mol.setTitle(pdb_id)

    # select source
    source_atoms = select_atoms(mol, pdb_id, cli_opts.source)
    target_atoms = select_atoms(mol, pdb_id, cli_opts.target)

    if cli_opts.verbose:
        print("Atoms in pdb file: {0}".format(mol.numAtoms()))
        print("Atoms of the source: {0}".format(source_atoms.numAtoms()))
        print("Atoms of the target: {0}".format(target_atoms.numAtoms()))

    # find the interacting pairs within the specified distance
    neighbors = bilab.structure.findNeighbors(
        source_atoms, cli_opts.distance, target_atoms)

    inter_atoms_set.set_interaction_pair(neighbors)
    if cli_opts.outfilename == "out":
        cli_opts.outfilename += "." + cli_opts.outfmt
    inter_atoms_set.write(o_file=cli_opts.outfilename)


if __name__ == '__main__':
    CURR_VERSION = sys.version_info
    if not CURR_VERSION[:2] > (2, 7):
        print("Python v.{}.{} is not supported".format(
            CURR_VERSION[0], CURR_VERSION[1]))
        sys.exit(0)
    # import bilab package
    try:
        import bilab
    except ImportError:
        raise ImportError("This script needs bilab package")

    # Create cli parser
    parser = argparse.ArgumentParser(
        description='parse PDB or mmcif format.',
        prefix_chars='-+/',
    )

    parser.add_argument("--pdb",
                        action="store",
                        dest='pdbfile',
                        help="The input pdb file[REQUIRED].")

    parser.add_argument("--source",
                        action="store",
                        dest="source",
                        help="The name of a source molecule in the pdb.")

    parser.add_argument("--target",
                        action="store",
                        dest="target",
                        help="The name of target molecule in the pdb.")

    parser.add_argument("--distance",
                        type=float,
                        dest="distance",
                        default=5.0,
                        help="The distance threshold for selecting "
                             "interacting pair [OPTION]. "
                             "Default is 5.0 angstroms.")

    parser.add_argument("--outfmt",
                        dest='outfmt',
                        default='txt',
                        help="output format: xml or txt. default is txt.")

    parser.add_argument("--out",
                        dest='outfilename',
                        default='out',
                        help="output file name. If not specified, "
                             "the name will be composed of out.fmt")

    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        dest="verbose",
                        help="print verbose info")

    parser.add_argument('--version', action='version',
                        version='%(prog)s 1.0')

    options = parser.parse_args()

    if not options.pdbfile:
        print("Error: do not specify an input pdb file")
        parser.print_help()
        sys.exit(1)

    if not options.source:
        print("Error: do not specify source name")
        parser.print_help()
        sys.exit(1)

    if not options.target:
        print("Error: do not specify target name")
        parser.print_help()
        sys.exit(1)

    f_ext = get_file_ext(options.pdbfile)
    if f_ext not in ['.pdb', '.cif']:
        print("Error: do not understand file format from the file extension")
        parser.print_help()
        sys.exit(1)
    if f_ext == '.pdb':
        options.fileloader = bilab.structure.parsePDB
    else:
        options.fileloader = bilab.io.cif_reader

    main(options)
