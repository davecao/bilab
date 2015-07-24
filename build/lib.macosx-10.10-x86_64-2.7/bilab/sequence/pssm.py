# -*- coding: utf-8 -*-
"""
   This module defines the PSSM class
"""

import os.path
import sys
import pprint
import numpy as np

from types import *

from .sequence import Sequence
from .psiblastIO import PsiBlast_pssm_score

__all__ = ['PSSM']

class PSSM(object):

    """
        class for storing a PSSM matrix and alignment info
    """

    def __init__(self, stream, format = 'xml', rawscorefile=None, sequence = None):
        super(PSSM, self).__init__()
        if format == 'xml':
            self.__parse_xml(stream)
        else:
            print("Unknow format: only xml format is supported now")

        if rawscorefile is not None:
            self.__parse_raw_score(rawscorefile)

    def __dictElement(self, element, prefix=None):
        """ This function is borrowed from prody """
        dict_ = {}
        length = False
        if isinstance(prefix, str):
            length = len(prefix)
        for child in element:
            tag = child.tag
            if length and tag.startswith(prefix):
                tag = tag[length:]
            if len(child) == 0:
                dict_[tag] = child.text
            else:
                dict_[tag] = child
        return dict_

    def __get_header(self, hit_id, hit_accession, hit_def):
        """ Analyse the header """
        dict_ = {}
        dbname = self._querydb
        if dbname == 'pdb':
            for item in (hit['id'] + hit['def']).split('>gi'):
                #>gi|1633462|pdb|4AKE|A Chain A, Adenylate Kinase
                #                        __________TITLE__________
                head, title = item.split(None, 1)
                head = head.split('|')
                pdb_id = head[-2].lower()
                chain_id = head[-1][:1]

                dict_['pdb_id'] = pdb_id
                dict_['chain_id'] = chain_id
                dict_['title'] = (head[-1][1:] + title).strip()
        elif dbname == 'cdd_delta':
            dict_['pssm-id'] = hit_accession
            dict_['accession'] = hit_def[:hit_def.find(',')]
            dict_['title'] = hit_def[hit_def.find(',')+1:].strip()
        elif dbname == 'nr':
            print("Parser for nr db is not implemented yet.")
            return None
        dict_ = {'db':dbname}
        return dict_

    def __parse_xml(self, xml, prefix = None):
        """ parse xml file/string"""

        import xml.etree.cElementTree as ET

        if os.path.isfile(xml):
            xml = ET.parse(xml)
            root = xml.getroot()
        elif isinstance(xml, StringTypes):
            pprint.pprint(xml)
            root = ET.fromstring(xml)
        else:
            raise ValueError('xml is not a filename or does not look like'
                                 ' a valid XML string')
        root = self.__dictElement(root, 'BlastOutput_')
        self._param = self.__dictElement(root['param'][0], 'Parameters_')
        self._querydb = root['db']
        self._program = root['program']
        self._program_version = root['version']

        self._querylen = int(root['query-len'])

        # Save
        hits = []
        for iteration in root['iterations']:
            for hit in self.__dictElement(iteration, 'Iteration_')['hits']:
                hit = self.__dictElement(hit, 'Hit_')
                data = self.__dictElement(hit['hsps'][0], 'Hsp_')
                for key in ['align-len', 'gaps', 'hit-frame', 'hit-from',
                            'hit-to', 'identity', 'positive', 'query-frame',
                            'query-from', 'query-to']:
                    data[key] = int(data[key])

                data['query-len'] = self._querylen

                for key in ['evalue', 'bit-score', 'score']:
                    data[key] = float(data[key])

                p_identity = 100.0 * data['identity'] / (data['query-to'] -
                                                    data['query-from'] + 1)
                # Overlap:
                p_overlap = (100.0 * (data['align-len'] - data['gaps']) /
                              self._querylen)
                data['percent_identity'] = p_identity
                #data['percent_coverage'] = p_overlap
                data['percent_overlap'] = p_overlap
                pdbch = data.copy()
                d = self.__get_header(hit['id'], hit['accession'], hit['def'])
                #pprint.pprint(d)
                pdbch.update(d)
                hits.append((p_identity, p_overlap, pdbch))

        hits.sort(key=lambda hit: hit[0], reverse=True)
        self._hits = hits
        #pprint.pprint(root)
        #pprint.pprint(self._param)
        #pprint.pprint(self._hits[0])

    def __parse_raw_score(self, score_file):
        """
            Read pssm score file in ascii code
            ._pssm_data contains keys:
                'score': np matrix of sequence in rows and score in columns.

                'obs_percent': np matrix same to 'score'

                'seq': a dictionary, seqno as keys and a.a. as value

                'order': order of amino acids for columns used in 'score' and
                        'obs_percent'

                'info_per_pos': a float value store information per position

                'relative_weight_gapless_pseudocounts': relative weight of
                              gapless real matches to pseudocounts
                'params':  paramters for K and lambda
        """
        self._pssm_data = PsiBlast_pssm_score(score_file)
        #np.savetxt('score_stdout.dat', self._pssm_data['score'], fmt="%2d")
