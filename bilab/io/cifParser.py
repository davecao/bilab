# -*- coding: utf-8 -*-
# cifParser.py

# Copyright (c) 2018, Wei Cao
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# * Neither the name of the copyright holders nor the names of any
#   contributors may be used to endorse or promote products derived
#   from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

"""This module defines functions for parsing and writing `PDB files`_ in mmcif format.

C++ extensions _mmcifio.so provide a dictionary of mmcif data

For example, '1a0s.cif' contains the following keys:

1. atom_site                      
    -- "group_PDB"         : ATOM                              
    -- "id"                : 1                              
    -- "type_symbol"       : N                               
    -- "label_atom_id"     : N                              
    -- "label_alt_id"      : .                              
    -- "label_comp_id"     : SER                              
    -- "label_asym_id"     : A                              
    -- "label_entity_id"   : 1                              
    -- "label_seq_id"      : 1                              
    -- "pdbx_PDB_ins_code" : ?                               
    -- "Cartn_x"           : -47.333                               
    -- "Cartn_y"           : 0.941                              
    -- "Cartn_z"           : 8.834                                
    -- "occupancy"         : 1.00                              
    -- "B_iso_or_equiv"    : 52.56                              
    -- "pdbx_formal_charge": ?                              
    -- "auth_seq_id"       : 71                              
    -- "auth_comp_id"      : SER                               
    -- "auth_asym_id"      : P                               
    -- "auth_atom_id"      : N                              
    -- "pdbx_PDB_model_num": 1  
2. atom_sites
3. atom_type
4. audit_author
5. audit_conform
6. cell                           
    -- "entry_id"        : 1A0S                              
    -- "length_a"        : 112.100                               
    -- "length_b"        : 112.100                               
    -- "length_c"        : 147.000                              
    -- "angle_alpha"     : 90.00                              
    -- "angle_beta"      : 90.00                                
    -- "angle_gamma"     : 120.00                              
    -- "Z_PDB"           : 9                              
    -- "pdbx_unique_axis": ?
7. chem_comp
8. citation
9. citation_author
10. database_2
11. database_PDB_matrix
12. diffrn
13. diffrn_detector
14. diffrn_radiation
15. diffrn_radiation_wavelength
16. diffrn_source
17. entity 
    -- "id"               : 1A0S
    -- "pdbx_description" : 
    -- "pdbx_fragment"    :
    -- "pdbx_ec"          :
    -- "pdbx_mutation"    :
    -- "details"          :
18. entity_poly
    -- "pdbx_strand_id" : P,Q,R
    -- "pdbx_target_identifier" : ?
    -- "pdbx_seq_one_letter_code" :
        [';SGFEFHGYARSGVIMNDSGASTKSGAYITPAGETGGAIGRLGNQADTYVEMNLEHKQTLDNGATTRFKVMVADGQTSYND\nWTASTSDLNVRQAFVELGNLPTFAGPFKGSTLWAGKRFDRDNFDIHWIDSDVVFLAGTGGGIYDVKWNDGLRSNFSLYGR\nNFGDIDDSSNSVQNYILTMNHFAGPLQMMVSGLRAKDNDERKDSNGNLAKGDAANTGVHALLGLHNDSFYGLRDGSSKTA\nLLYGHGLGAEVKGIGSDGALRPGADTWRIASYGTTPLSENWSVAPAMLAQRSKDRYADGDSYQWATFNLRLIQAINQNFA\nLAYEGSYQYMDLKPEGYNDRQAVNGSFYKLTFAPTFKVGSIGDFFSRPEIRFYTSWMDWSKKLNNYASDDALGSDGFNSG\nGEWSFGVQMETWF\n;']   
19. entity_poly_seq

20. entity_src_gen
21. entry
22. exptl
    -- "entry_id"       : 1A0S
    -- "method"         : X-RAY DIFFRACTION
    -- "crystals_number": 1

23. exptl_crystal
24. exptl_crystal_grow
25. pdbx_audit_revision_category
26. pdbx_audit_revision_details
27. pdbx_audit_revision_group
28. pdbx_audit_revision_history
29. pdbx_audit_revision_item
30. pdbx_database_status          
    -- "status_code" : REL                              
    -- "entry_id"    : 1A0S                               
    -- "recvd_initial_deposition_date" : 1997-12-07
31. pdbx_entity_nonpoly
32. pdbx_nonpoly_scheme
33. pdbx_poly_seq_scheme
34. pdbx_struct_assembly
35. pdbx_struct_assembly_gen
36. pdbx_struct_assembly_prop
37. pdbx_struct_conn_angle
38. pdbx_struct_oper_list
39. pdbx_struct_sheet_hbond
40. pdbx_validate_planes
41. pdbx_validate_rmsd_angle
42. pdbx_validate_symm_contact
43. pdbx_validate_torsion
44. pdbx_xplor_file
45. refine
    -- "entry_id"           : 1A0S
    -- "ls_d_res_low"       : 100.
    -- "ls_d_res_high"      : 2.4 <--------REMARK  2 RESOLUTION.  2.4 ANGSTROMS
    -- "ls_R_factor_R_work" : 0.2140000
    -- "ls_R_factor_R_free" : 0.2280000
    -- "B_iso_mean"         : 27.5
46. refine_analyze
47. refine_hist
48. refine_ls_restr
49. refine_ls_restr_ncs
50. refine_ls_shell
51. reflns
52. reflns_shell
53. software
    -- "name" : 
54. source
55. struct
    -- "entry_id"                : 1A0S
    -- "title"                   : SUCROSE-SPECIFIC PORIN
    -- "pdbx_descriptor"         : SUCROSE-SPECIFIC PORIN
    -- "pdbx_model_details"      : ?
    -- "pdbx_CASP_flag"          : ?
    -- "pdbx_model_type_details" : ?
56. struct_asym
57. struct_biol
58. struct_conf
59. struct_conf_type
60. struct_conn
61. struct_conn_type
62. struct_keywords
    -- "entry_id"      : 1A0S
    -- "pdbx_keywords" : OUTER MEMBRANE PROTEIN
    -- "text"          : OUTER MEMBRANE PROTEIN, PORIN
63. struct_ncs_dom
64. struct_ncs_ens
65. struct_ncs_oper
    -- "id"           : 1
    -- "code"         : given
    -- "details"      : ?
    -- "matrix[1][1]" : -0.498390
    -- "matrix[1][2]" :  0.864212
    -- "matrix[1][3]" :  0.068880
    -- "matrix[2][1]" : -0.866499
    -- "matrix[2][2]" : -0.499125
    -- "matrix[2][3]" : -0.007322
    -- "matrix[3][1]" :  0.028052
    -- "matrix[3][2]" : -0.063334
    -- "matrix[3][3]" :  0.997598
    -- "vector[1]"    : -54.44100
    -- "vector[2]"    : -31.93100
    -- "vector[3]"    :   1.00740
66. struct_ref
67. struct_ref_seq
68. struct_sheet
69. struct_sheet_order
70. struct_sheet_range
71. struct_site
72. struct_site_gen
73. symmetry



>> fn = '1a0s.cif'
>> cif_data = bilab.io._mmcifio.reader_mmcif(fn)
>> print("group_PDB: {}".format(cif_data['atom_site']['group_PDB'][0]))
>> print("id:{}".format(cif_data['atom_site']['id'][0]))
>> print("type_symbol: {}".format(cif_data['atom_site']['type_symbol'][0]))
>> print("label_atom_id: {}".format(cif_data['atom_site']['label_atom_id'][0]))
>> print("label_alt_id: {}".format(cif_data['atom_site']['label_alt_id'][0]))
>> print("label_comp_id: {}".format(cif_data['atom_site']['label_comp_id'][0]))
>> print("label_asym_id: {}".format(cif_data['atom_site']['label_asym_id'][0]))
>> print("label_entity_id: {}".format(cif_data['atom_site']['label_entity_id'][0]))
>> print("label_seq_id: {}".format(cif_data['atom_site']['label_seq_id'][0]))
>> print("pdbx_PDB_ins_code: {}".format(cif_data['atom_site']['pdbx_PDB_ins_code'][0]))
>> print("Cartn_x: {}".format(cif_data['atom_site']['Cartn_x'][0]))
>> print("Cartn_y: {}".format(cif_data['atom_site']['Cartn_y'][0]))
>> print("Cartn_z: {}".format(cif_data['atom_site']['Cartn_z'][0]))
>> print("occupancy: {}".format(cif_data['atom_site']['occupancy'][0]))
>> print("B_iso_or_equiv: {}".format(cif_data['atom_site']['B_iso_or_equiv'][0]))
>> print("pdbx_formal_charge: {}".format(cif_data['atom_site']['pdbx_formal_charge'][0]))
>> print("auth_seq_id: {}".format(cif_data['atom_site']['auth_seq_id'][0]))
>> print("auth_comp_id: {}".format(cif_data['atom_site']['auth_comp_id'][0]))
>> print("auth_asym_id: {}".format(cif_data['atom_site']['auth_asym_id'][0]))
>> print("auth_atom_id: {}".format(cif_data['atom_site']['auth_atom_id'][0]))
>> print("pdbx_PDB_model_num: {}".format(cif_data['atom_site']['pdbx_PDB_model_num'][0]))
"""
#from __future__ import division, print_function

from collections import defaultdict
import os.path

import numpy as np

from bilab.aaprop.aminoacids import AminoAcids_interity
from bilab.structure.atomic import AtomGroup
from bilab.structure.atomic import flags
from bilab.structure.atomic import ATOMIC_FIELDS
from bilab.utilities import openFile
from bilab import LOGGER, SETTINGS

from bilab.structure.header import Chemical
from bilab.structure.header import Polymer
from bilab.structure.header import DBRef
from bilab.structure.header import buildBiomolecules
from bilab.structure.header import assignSecstr

try:
    from bilab.io._mmcifio import reader_mmcif as read_cif
except ImportError:
    raise ImportError('read_cif() function could not be imported. '
                      'Reinstall bilab may solve the problem.')

__all__ = ['cif_reader']

__version__ = '2018.07.03'

_PDBSubsets = {'ca': 'ca', 'calpha': 'ca', 'bb': 'bb', 'backbone': 'bb'}


def get_pdb_info(data, *keys):
    pass
                


def cif_reader(cif, **kwargs):
    """Return an :class:`.AtomGroup` and/or dictionary containing header data
    parsed from a PDB file in MMCIF format.

    :arg cif: a filename
    """
    if not os.path.isfile(cif):
        raise IOError('{0} is not a valid filename or a valid PDB (cif)'
                      'identifier.'.format(pdb))
    data = read_cif(cif)
    model = kwargs.get('model')
    header = kwargs.get('header', False)
    assert isinstance(header, bool), 'header must be a boolean'
    chain = kwargs.get('chain')
    subset = kwargs.get('subset')
    altloc = kwargs.get('altloc', '.')
    cif_category = kwargs.get('category','atom_site')
    
    if model is not None:
        if isinstance(model, int):
            if model < 0:
                raise ValueError('model must be greater than 0')
        else:
            raise TypeError('model must be an integer, {0} is invalid'
                            .format(str(model)))
    title_suffix = ''
    if subset:
        try:
            subset = _PDBSubsets[subset.lower()]
        except AttributeError:
            raise TypeError('subset must be a string')
        except KeyError:
            raise ValueError('{0} is not a valid subset'
                             .format(repr(subset)))

    if chain is not None:
        if not isinstance(chain, str):
            raise TypeError('chain must be a string')
        elif len(chain) == 0:
            raise ValueError('chain must not be an empty string')
        title_suffix = '_' + chain + title_suffix
    
    ag = None
    if 'ag' in kwargs:
        ag = kwargs['ag']
        if not isinstance(ag, AtomGroup):
            raise TypeError('ag must be an AtomGroup instance')
        n_csets = ag.numCoordsets()
    elif model != 0:
        ag = AtomGroup(str(kwargs.get('title', 'Unknown')) + title_suffix)
        n_csets = 0

    biomol = kwargs.get('biomol', False)
    auto_secondary = None
    secondary = kwargs.get('secondary')
    if not secondary:
        auto_secondary = SETTINGS.get('auto_secondary')
        secondary = auto_secondary
    
    split = 0
    hd = None
    if model != 0:
        LOGGER.timeit()
        # Read items in category "atom_site"
        try:
            atom_dict = data[cif_category]
        except AttributeError as err:
            raise err
        if not len(atom_dict):
            raise ValueError('empty mmcif file or stream')
        if header or biomol or secondary:
            hd, split = get_pdb_info(data)
    
        _addAtomSiteInfo(ag, atom_dict, model, chain, subset, altloc)

        if ag.numAtoms() > 0:
            LOGGER.report(
                '{0} atoms and {1} coordinate set(s) were '
                'parsed in %.2fs.'.format(
                    ag.numAtoms(),
                    ag.numCoordsets() - n_csets))
        else:
            ag = None
            LOGGER.warn('Atomic data could not be parsed, please '
                            'check the input file.')
    elif header:
        hd, split = get_pdb_info(data)
    
    if ag is not None and isinstance(hd, dict):
        if secondary:
            if auto_secondary:
                try:
                    ag = assignSecstr(hd, ag)
                except ValueError:
                    pass
            else:
                ag = assignSecstr(hd, ag)
        if biomol:
            ag = buildBiomolecules(hd, ag)

            if isinstance(ag, list):
                LOGGER.info('Biomolecular transformations were applied, {0} '
                            'biomolecule(s) are returned.'.format(len(ag)))
            else:
                LOGGER.info(
                    'Biomolecular transformations were applied to the '
                    'coordinate data.')
    
    if model != 0:
        if header:
            return ag, hd
        else:
            return ag
    else:
        return hd


def _addAtomSiteInfo(atomgroup, atom_dict, model, chain, subset, altloc_torf, format='PDB'):
    """Return an AtomGroup.

    :arg atom_dict: a dictionary variable for category "atom_site" in mmcif."""
    format = format.upper()
    if format == 'PDB':
        isPDB = True
    else:
        isPDB = False
    
    if subset:
        if subset == 'ca':
            subset = set(('CA',))
        elif subset in 'bb':
            subset = flags.BACKBONE
            only_subset = True
            protein_resnames = flags.AMINOACIDS
    else:
        only_subset = False
    
    if chain is None:
        only_chains = False
    else:
        only_chains = True
    onlycoords = False
    n_atoms = atomgroup.numAtoms()
    if n_atoms > 0:
        asize = n_atoms
    else:
        asize = len(atom_dict['group_PDB'])
    addcoords = False
    if atomgroup.numCoordsets() > 0:
        addcoords = True

    alength = asize
    coordinates = np.zeros((asize, 3), dtype=float)
    atomnames = np.zeros(asize, dtype=ATOMIC_FIELDS['name'].dtype)
    resnames = np.zeros(asize, dtype=ATOMIC_FIELDS['resname'].dtype)
    resnums = np.zeros(asize, dtype=ATOMIC_FIELDS['resnum'].dtype)
    chainids = np.zeros(asize, dtype=ATOMIC_FIELDS['chain'].dtype)
    hetero = np.zeros(asize, dtype=bool)
    termini = np.zeros(asize, dtype=bool)
    altlocs = np.zeros(asize, dtype=ATOMIC_FIELDS['altloc'].dtype)
    icodes = np.zeros(asize, dtype=ATOMIC_FIELDS['icode'].dtype)
    serials = np.zeros(asize, dtype=ATOMIC_FIELDS['serial'].dtype)
    if isPDB:
        segnames = np.zeros(asize, dtype=ATOMIC_FIELDS['segment'].dtype)
        elements = np.zeros(asize, dtype=ATOMIC_FIELDS['element'].dtype)
        bfactors = np.zeros(asize, dtype=ATOMIC_FIELDS['beta'].dtype)
        occupancies = np.zeros(asize, dtype=ATOMIC_FIELDS['occupancy'].dtype)
        anisou = None
        siguij = None
    else:
        charges = np.zeros(asize, dtype=ATOMIC_FIELDS['charge'].dtype)
        radii = np.zeros(asize, dtype=ATOMIC_FIELDS['radius'].dtype)
    
    total_atoms = len(atom_dict['group_PDB']) 
    asize = 2000  # increase array length by this much when needed
    start = 0
    stop = total_atoms  #  len(atom_dict['group_PDB'])
    
    nmodel = 0
    # if a specific model is requested, skip lines until that one
    if isPDB and model is not None and model != 1:
        for i in range(total_atoms):
            if atom_dict['pdbx_PDB_model_num'][i] != nmodel:
                nmodel += 1
                if model == nmodel:
                    start = i + 1
                    stop = n_atoms
                    break
        if nmodel != model:
            raise PDBParseError('model {0} is not found'.format(model))
    if isinstance(altloc_torf, str):
        if altloc_torf.strip() != '.':
            LOGGER.info('Parsing alternate locations {0}.'
                        .format(altloc_torf))
            which_altlocs = ' ' + ''.join(altloc_torf.split())
        else:
            which_altlocs = ' .'
        altloc_torf = False
    else:
        which_altlocs = ' .'
        altloc_torf = True
    
    acount = 0
    altloc = defaultdict(list)
    i = start
    END = False
    while i < stop:
        startswith = atom_dict['group_PDB'][i]
        #if nmodel != atom_dict['pdbx_PDB_model_num']:
        #    startswith = "END"
        if startswith == 'ATOM' or startswith == 'HETATM':
            if only_subset:
                atomname = atom_dict['auth_atom_id'][i]
                resname = atom_dict['auth_comp_id'][i]
                if not (atomname in subset and resname in protein_resnames):
                    i += 1
                    continue
            else:
                atomname = atom_dict['auth_atom_id'][i]
                resname = atom_dict['auth_comp_id'][i]
            
            chid = atom_dict['auth_asym_id'][i]
            if only_chains:
                if chid not in chain:
                    i += 1
                    continue
            alt = atom_dict['label_alt_id'][i]
            if alt not in which_altlocs:
                altloc[alt].append(i)
                i += 1
                continue
            try:
                coordinates[acount, 0] = atom_dict['Cartn_x'][i]
                coordinates[acount, 1] = atom_dict['Cartn_y'][i]
                coordinates[acount, 2] = atom_dict['Cartn_z'][i]
            except:
                if acount >= n_atoms > 0:
                    if nmodel == 0:
                        raise ValueError(format + 'file and AtomGroup ag must '
                                         'have same number of atoms')
                    LOGGER.warn(
                        'Discarding model {0}, which contains more '
                        'atoms than first model does.'.format(nmodel + 1))
                    acount = 0
                    nmodel += 1
                    coordinates = np.zeros((n_atoms, 3), dtype=float)
                else:
                    raise ("Invalid or missing coordinate(s) at {}".format(i))
            if onlycoords:
                acount += 1
                i += 1
                continue
            try:
                serials[acount] = atom_dict['id'][i]
            except ValueError:
                LOGGER.warn('Failed to parse serial number {0}.'
                                .format(i))
                serials[acount] = serials[acount - 1] + 1
            
            altlocs[acount] = alt
            atomnames[acount] = atomname
            resnames[acount] = resname
            chainids[acount] = chid
            resnums[acount] = atom_dict['auth_seq_id'][i]  # .split()[0])
            icodes[acount] = atom_dict['pdbx_PDB_ins_code'][i]
            if isPDB:
                try:
                    occupancies[acount] = atom_dict['occupancy'][i]
                except:
                    LOGGER.warn('failed to parse occupancy at line {0}'
                                .format(i))
                try:
                    bfactors[acount] = atom_dict['B_iso_or_equiv'][i]
                except:
                    LOGGER.warn('failed to parse beta-factor at line {0}'
                                .format(i))
                hetero[acount] = startswith[0] == 'H'
                segnames[acount] = ""
                elements[acount] = atom_dict['type_symbol'][i]
            else:
                try:
                    charges[acount] = "" # not found corresponding part to atom_dict['']
                except:
                    LOGGER.warn('failed to parse charge at line {0}'
                                .format(i))
                try:
                    radii[acount] = "" #line[62:69]
                except:
                    LOGGER.warn('failed to parse radius at line {0}'
                                .format(i))
            acount += 1  
            #if n_atoms == 0 and acount >= alength:
            #    # if arrays are short extend them with zeros
            #    alength += asize
            #    coordinates = np.concatenate(
            #        (coordinates, np.zeros((asize, 3), float)))
            #    atomnames = np.concatenate((
            #        atomnames,
            #        np.zeros(asize, ATOMIC_FIELDS['name'].dtype)))
            #    resnames = np.concatenate((
            #        resnames,
            #        np.zeros(asize, ATOMIC_FIELDS['resname'].dtype)))
            #    resnums = np.concatenate((
            #        resnums,
            #        np.zeros(asize, ATOMIC_FIELDS['resnum'].dtype)))
            #    chainids = np.concatenate((
            #        chainids,
            #        np.zeros(asize, ATOMIC_FIELDS['chain'].dtype)))
            #    hetero = np.concatenate((hetero, np.zeros(asize, bool)))
            #    termini = np.concatenate((termini, np.zeros(asize, bool)))
            #    altlocs = np.concatenate((
            #        altlocs,
            #        np.zeros(asize, ATOMIC_FIELDS['altloc'].dtype)))
            #    icodes = np.concatenate((
            #        icodes,
            #        np.zeros(asize, ATOMIC_FIELDS['icode'].dtype)))
            #    serials = np.concatenate((
            #        serials,
            #        np.zeros(asize, ATOMIC_FIELDS['serial'].dtype)))
            #    if isPDB:
            #        bfactors = np.concatenate((
            #            bfactors,
            #            np.zeros(asize, ATOMIC_FIELDS['beta'].dtype)))
            #        occupancies = np.concatenate((
            #            occupancies,
            #            np.zeros(asize, ATOMIC_FIELDS['occupancy'].dtype)))
            #        segnames = np.concatenate((
            #            segnames,
            #            np.zeros(asize, ATOMIC_FIELDS['segment'].dtype)))
            #        elements = np.concatenate((
            #            elements,
            #            np.zeros(asize, ATOMIC_FIELDS['element'].dtype)))
            #        if anisou is not None:
            #            anisou = np.concatenate((
            #                anisou,
            #                np.zeros(
            #                    (asize, 6),
            #                    ATOMIC_FIELDS['anisou'].dtype)))
            #        if siguij is not None:
            #            siguij = np.concatenate((
            #                siguij,
            #                np.zeros(
            #                    (asize, 6),
            #                    ATOMIC_FIELDS['siguij'].dtype)))
            #    else:
            #        charges = np.concatenate((
            #            charges,
            #            np.zeros(asize, ATOMIC_FIELDS['charge'].dtype)))
            #        radii = np.concatenate((
            #            radii,
            #            np.zeros(asize, ATOMIC_FIELDS['radius'].dtype)))
        elif not onlycoords and nmodel != atom_dict['pdbx_PDB_model_num'][i]:  
            # (startswith == 'TER   ' or
            # startswith.strip() == 'TER'):
            termini[acount - 1] = True
        elif nmodel != atom_dict['pdbx_PDB_model_num'][i]:  # Start with new model
            # elif startswith == 'ENDMDL' or startswith[:3] == 'END':
            if acount == 0:
                # If there is no atom record between ENDMDL & END skip to next
                i += 1
                continue
            if model is not None:
                i += 1
                break
            diff = stop - i - 1
            if diff < acount:
                END = True
            if onlycoords:
                if acount < n_atoms:
                    LOGGER.warn('Discarding model {0}, which contains '
                                '{1} fewer atoms than the first model '
                                'does.'.format(nmodel+1, n_atoms-acount))
                else:
                    coordsets[nmodel] = coordinates
                    nmodel += 1
                acount = 0
                if not END:
                    coordinates = coordsets[nmodel]
            else:
                if acount != n_atoms > 0:
                    raise ValueError(
                        'PDB file and AtomGroup ag must have '
                        'same number of atoms')
                # this is where to decide if more coordsets should be expected
                if END:
                    coordinates.resize((acount, 3), refcheck=False)
                    if addcoords:
                        atomgroup.addCoordset(coordinates)
                    else:
                        atomgroup._setCoords(coordinates)
                else:
                    coordsets = np.zeros((diff/acount+1, acount, 3))
                    coordsets[0] = coordinates[:acount]
                    onlycoords = True
                atomnames.resize(acount, refcheck=False)
                resnames.resize(acount, refcheck=False)
                resnums.resize(acount, refcheck=False)
                chainids.resize(acount, refcheck=False)
                hetero.resize(acount, refcheck=False)
                termini.resize(acount, refcheck=False)
                altlocs.resize(acount, refcheck=False)
                icodes.resize(acount, refcheck=False)
                serials.resize(acount, refcheck=False)
                if not only_subset:
                    atomnames = np.char.strip(atomnames)
                    resnames = np.char.strip(resnames)
                atomgroup.setNames(atomnames)
                atomgroup.setResnames(resnames)
                atomgroup.setResnums(resnums)
                atomgroup.setChids(chainids)
                atomgroup.setFlags('hetatm', hetero)
                atomgroup.setFlags('pdbter', termini)
                atomgroup.setAltlocs(altlocs)
                atomgroup.setIcodes(np.char.strip(icodes))
                atomgroup.setSerials(serials)
                if isPDB:
                    bfactors.resize(acount, refcheck=False)
                    occupancies.resize(acount, refcheck=False)
                    segnames.resize(acount, refcheck=False)
                    elements.resize(acount, refcheck=False)
                    atomgroup.setBetas(bfactors)
                    atomgroup.setOccupancies(occupancies)
                    atomgroup.setSegnames(np.char.strip(segnames))
                    atomgroup.setElements(np.char.strip(elements))
                    if anisou is not None:
                        anisou.resize((acount, 6))
                        atomgroup.setAnisous(anisou / 10000)
                    if siguij is not None:
                        siguij.resize((acount, 6))
                        atomgroup.setAnistds(siguij / 10000)
                else:
                    charges.resize(acount)
                    radii.resize(acount)
                    atomgroup.setCharges(charges)
                    atomgroup.setRadii(radii)

                nmodel += 1
                n_atoms = acount
                acount = 0
                coordinates = np.zeros((n_atoms, 3), dtype=float)
                if altloc and altloc_torf:
                    _evalAltlocs(atomgroup, altloc, chainids, resnums,
                                 resnames, atomnames)
                    altloc = defaultdict(list)
                if END:
                    break
        elif isPDB and startswith == 'ANISOU':
            if anisou is None:
                anisou = True
                anisou = np.zeros(
                    (alength, 6),
                    dtype=ATOMIC_FIELDS['anisou'].dtype)
            try:
                index = acount - 1
                anisou[index, 0] = None  # line[28:35] not used now
                anisou[index, 1] = None  # line[35:42] not used now
                anisou[index, 2] = None  # line[43:49] not used now
                anisou[index, 3] = None  # line[49:56] not used now
                anisou[index, 4] = None  # line[56:63] not used now
                anisou[index, 5] = None  # line[63:70] not used now
            except:
                LOGGER.warn(
                    'failed to parse anisotropic temperature '
                    'factors at line {0}'.format(i))
        elif isPDB and startswith == 'SIGUIJ':
            if siguij is None:
                siguij = np.zeros(
                    (alength, 6),
                    dtype=ATOMIC_FIELDS['siguij'].dtype)
            try:
                index = acount - 1
                siguij[index, 0] = None # line[28:35] not used now
                siguij[index, 1] = None # line[35:42] not used now
                siguij[index, 2] = None # line[43:49] not used now
                siguij[index, 3] = None # line[49:56] not used now
                siguij[index, 4] = None # line[56:63] not used now
                siguij[index, 5] = None # line[63:70] not used now
            except:
                LOGGER.warn(
                    'failed to parse standard deviations of '
                    'anisotropic temperature factors at line {0}'.format(i))
        elif startswith == 'SIGATM':
            pass
        i += 1
    if onlycoords:
        if acount == atomgroup.numAtoms():
            coordsets[nmodel] = coordinates
            nmodel += 1
        del coordinates
        coordsets.resize((nmodel, atomgroup.numAtoms(), 3))
        if addcoords:
            atomgroup.addCoordset(coordsets)
        else:
            atomgroup._setCoords(coordsets)
    elif not END:
        # this means last line was an ATOM line, so atomgroup is not decorated
        if acount < total_atoms:
            coordinates.resize((acount, 3))
        if addcoords:
            atomgroup.addCoordset(coordinates)
        else:
            atomgroup._setCoords(coordinates)
        if acount < total_atoms:
            atomnames.resize(acount)
            resnames.resize(acount)
            resnums.resize(acount)
            chainids.resize(acount)
            hetero.resize(acount)
            termini.resize(acount)
            altlocs.resize(acount)
            icodes.resize(acount)
            serials.resize(acount)
        if not only_subset:
            atomnames = np.char.strip(atomnames)
            resnames = np.char.strip(resnames)
        atomgroup.setNames(atomnames)
        atomgroup.setResnames(resnames)
        atomgroup.setResnums(resnums)
        atomgroup.setChids(chainids)
        atomgroup.setFlags('hetatm', hetero)
        atomgroup.setFlags('pdbter', termini)
        atomgroup.setAltlocs(altlocs)
        atomgroup.setIcodes(np.char.strip(icodes))
        atomgroup.setSerials(serials)
        if isPDB:
            if anisou is not None:
                anisou.resize((acount, 6))
                atomgroup.setAnisous(anisou / 10000)
            if siguij is not None:
                siguij.resize((acount, 6))
                atomgroup.setAnistds(siguij / 10000)
            if acount < total_atoms:
                bfactors.resize(acount)
                occupancies.resize(acount)
                segnames.resize(acount)
                elements.resize(acount)
            atomgroup.setSegnames(np.char.strip(segnames))
            atomgroup.setElements(np.char.strip(elements))
            atomgroup.setBetas(bfactors)
            atomgroup.setOccupancies(occupancies)
        else:
            # problem?
            charges.resize(acount)
            radii.resize(acount)
            atomgroup.setCharges(charges)
            atomgroup.setRadii(radii)

    if altloc and altloc_torf:
        _evalAltlocs(atomgroup, altloc, chainids, resnums, resnames, atomnames)

    return atomgroup


def _evalAltlocs(atomgroup, altloc, chainids, resnums, resnames, atomnames):
    altloc_keys = list(altloc)
    altloc_keys.sort()
    indices = {}
    for key in altloc_keys:
        xyz = atomgroup.getCoords()
        success = 0
        lines = altloc[key]
        for line, i in lines:
            aan = line[12:16].strip()
            arn = line[17:21].strip()
            ach = line[21]
            ari = int(line[22:26].split()[0])
            rn, ids, ans = indices.get((ach, ari), (None, None, None))
            if ids is None:
                ids = indices.get(ach, None)
                if ids is None:
                    ids = (chainids == ach).nonzero()[0]
                    indices[ach] = ids
                ids = ids[resnums[ids] == ari]
                if len(ids) == 0:
                    LOGGER.warn("failed to parse altloc {0} at line {1}, "
                                "residue not present for altloc 'A'".format(
                                    repr(key), i+1))
                    continue
                rn = resnames[ids[0]]
                ans = atomnames[ids]
                indices[(ach, ari)] = (rn, ids, ans)
            if rn != arn:
                LOGGER.warn("failed to parse altloc {0} at line {1}, "
                            "residue name mismatch (expected {2}, "
                            "parsed {3})".format(repr(key), i+1, repr(rn),
                                                 repr(arn)))
                continue
            index = ids[(ans == aan).nonzero()[0]]
            if len(index) != 1:
                LOGGER.warn("failed to parse altloc {0} at line {1}, atom"
                            " {2} not found in the residue"
                            .format(repr(key), i+1, repr(aan)))
                continue
            try:
                xyz[index[0], 0] = float(line[30:38])
                xyz[index[0], 1] = float(line[38:46])
                xyz[index[0], 2] = float(line[46:54])
            except:
                LOGGER.warn('failed to parse altloc {0} at line {1}, could'
                            ' not read coordinates'.format(repr(key), i+1))
                continue
            success += 1
        LOGGER.info('{0} out of {1} altloc {2} lines were parsed.'
                    .format(success, len(lines), repr(key)))
        if success > 0:
            LOGGER.info('Altloc {0} is appended as a coordinate set to the '
                        'atom group.'.format(repr(key), atomgroup.getTitle()))
            atomgroup.addCoordset(xyz, label='altloc ' + key)

