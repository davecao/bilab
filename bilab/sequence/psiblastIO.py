# -*- coding: utf-8 -*-

from __future__ import print_function

import os.path
import re
import pprint
import numpy as np

__all__ = ['PsiBlast_pssm_score']

def PsiBlast_pssm_score(handle):
    """
      Parse pssm ascii output
        order='ARNDCQEGHILKMFPSTWYV' output from psiblast
    """
    order = 'ARNDCQEGHILKMFPSTWYV'
    pssm_scores = []
    pssm_obs_percent = []
    # information per position,
    info_per_pos = []
    # relative weight of gapless real matches to pseudocounts
    relative_weight_gapless_pseudocounts = []
    aa_in_pos = {}
    parameters = {}
    order_len = len(order)
    result = {}
    if os.path.isfile(handle):
            with open(handle, 'r') as fh:
                for line in fh:
                    line = line.strip("\n")

                    if line.startswith('Last') or (not line):
                        continue

                    if line.isalpha():
                        if re.search("Lambda", line) is None:
                            # order line
                            given_order = ''.join(line.split())
                            if given_order[:order_len] != order:
                                order = given_order[:order_len]
                                print("order is different")
                            else:
                                print("order is same")
                        continue
                    elif re.search("Lambda", line):
                        continue

                    content = line.split()

                    if content[0].isdigit():
                        # matrix line
                        aa_in_pos[content[0]] = content[1]
                        data = [int(s) for s in content[2:order_len + 2]]
                        pssm_scores.append(data)
                        data = [int(s) for s in \
                                   content[order_len +2: order_len*2 + 2]]
                        pssm_obs_percent.append(data)
                        data = float(content[-2])
                        info_per_pos.append(data)
                        data = float(content[-1])
                        relative_weight_gapless_pseudocounts.append(data)
                    else:
                        # parameters
                        parameters['K'+'.'+content[0]+'_'+ content[1]] = \
                                         content[2]
                        parameters['Lambda'+'.'+content[0]+'_'+ content[1]] = \
                                        content[3]
    else:
        raise ValueError('handle is not a filename or does not exist')
    result['score'] = np.asmatrix(pssm_scores)
    result['obs_percent'] = np.asmatrix(pssm_obs_percent)
    result['seq'] = aa_in_pos
    result['order'] = order
    result['info_per_pos'] = info_per_pos
    result['relative_weight_gapless_pseudocounts'] = \
                              relative_weight_gapless_pseudocounts
    result['params'] = parameters
    return result

