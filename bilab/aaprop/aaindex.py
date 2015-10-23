# -*- coding: utf-8 -*-

import os
import re
import math
import pprint

from types import *
from bisect import bisect_left
import bilab

__all__ = [ 'AAindex', 'AAindex1Parser', 'AAindex2Parser' ]


def find_name_in_aaindex(lst, i):
    """
        binary search return index
    """
    index = bisect_left(lst, i)
    if index == len(lst) or lst[index] != i:
        return False
    return index

def str2float(str):
    #for string NA, N/A
    NaN_str = "NA:N/A:None"
    found = re.search(str, NaN_str)
    if found:
        val = float('NaN')
    else:
        try:
            val = float(str)
        except ValueError:
            raise ValueError("Failed to convert to float")
    return val

def string2dict(str_list, state = 0):
    """ string to dictionary """
    if type(str_list) is not StringType:
        raise TypeError("Input should be a string")
    x = str_list.split()
    return dict(zip(x[0::2], x[1::2]))

def index2dict(str_list, state = 0):
    """ Convert aaindex values into dictionary """
    order = 'ARNDCQEGHILKMFPSTWYV'
    dict_ = {}

    if type(str_list) is not StringType:
        raise TypeError("Input should be a string")
    str = "".join(str_list.replace("/", '').split())

    if str.isalpha():
        return dict_

    for inx, val in enumerate(str_list.split()):
        if state == 2:
            dict_[order[inx+10]] = str2float(val)
        elif state == 1:
            dict_[order[inx]] = str2float(val)

    return dict_

def literature2dict(str_list, separator=":", state=0):
    dict_ = {}
    if type(str_list) is not StringType:
        raise TypeError("Input should be a string")
    str_list = re.sub(r'\s*:\s*',':', str_list)
    if not str_list:
        return dict_
    for liter in str_list.split():
        k, v = liter.split(separator)
        dict_[k] = v
    return dict_

def extract_field(line, field, separator=" "):
    return line.split(separator)[field]

# 'C': lambda x : dict(zip(x.split()[0::2], x.split()[1::2]))

AAINDEX_TAGS={
    'C': string2dict,
    'R': literature2dict,
    'I': index2dict
}

class AAindex1Parser(object):
    """
    Parse aaindex1 format

    """
    def __init__(self):
        """ initialize """
        self._data_location = bilab.data + os.sep + 'aaindex' + os.sep + 'aaindex1'
        self._parsername = 'AAindex1Parser'
        self._data = self.__dictRecord(self._data_location)
        self._scale_name_ind = [x['H'] for x in self._data]

    def __call__(self):
        """
            When directly call from obj.
            e.g., aa = AAindex1Parser(), aa()
        """
        print("aaindex1 infile:{} from __call__".format(self._data_location))
        #self.__parse_record()

    def __dictRecord(self, infile, delimitor="//"):
        """ Read into dict """
        data = []
        dict_ = {}
        old_tag = ''
        tag = ''
        state_continue = 0
        with open(self._data_location,'r') as fhandle:
            for line in fhandle:
                #print("line:{}".format(line))
                if line.startswith('//'):
                    #dict_['C'].update(dict(dict_['C'].split()))
                    data.append(dict_.copy())
                    state_continue = 0
                    # copy data
                    dict_.clear()
                    old_tag = tag
                elif line.startswith(" "):
                    # continue line
                    state_continue += 1
                    if tag in AAINDEX_TAGS:
                        list_dict = AAINDEX_TAGS[tag](line[1:].strip(),
                                     state = state_continue)
                        dict_[tag].update(list_dict)
                    else:
                        dict_[tag] = "{} {}".format(dict_[tag], line.strip())
                else:
                    tag = line.strip()[0]

                    if tag != old_tag:
                        # New item
                        old_tag = tag
                        state_continue = 0

                    if tag in AAINDEX_TAGS:
                        if tag in dict_:
                            list_dict = AAINDEX_TAGS[tag](line[1:].strip(),
                                state = state_continue)
                            dict_[tag].update(list_dict)
                        else:
                            dict_[tag] = AAINDEX_TAGS[tag](line[1:].strip(),
                                state = state_continue)
                    else:
                        # process as string
                        dict_[tag] = line[1:].strip()
        data.sort(key=lambda x: x['H'])
        return data

    def parse_record(self, infile, verbose=False):
        """ Parse line

        Args:
            infile(str): the name of a file
            verbose(bool): show detail info
        """
        if verbose:
            print("aaindex1 infile:{}".format(self._data_location))
        self._data = self.__dictRecord(infile)
        self._scale_name_ind = [ x['H'] for x in self._data ]

    def get_parser_name(self):
        return self._parsername

    def get_scale(self, scale_name):
        """ Return scale values in dict """
        inx = find_name_in_aaindex(self._scale_name_ind, scale_name.upper())
        #pprint.pprint(inx)
        #pprint.pprint(self._data)
        #pprint.pprint(self._data[inx])
        if inx:
            return self._data[inx]['I']
        print("Could not found the scale named {}\nCheck the scale_name please"
            .format(scale_name))
        return None

class AAindex2Parser(object):
    """ Parse aaindex2 format
        Not implemented yet

    """
    def __init__(self):
        """ initialize """
        self._data_location = bilab.data + os.sep + 'aaindex' + os.sep + 'aaindex2'
        self._parsername = 'AAindex2Parser'
        self.parse_record()

    def __call__(self):
        """
            When directly call from obj.
            e.g., aa = AAindex2Parser(), aa()
        """
        print("aaindex2 infile:{} from __call__".format(self._data_location))

    def parse_record(self, verbose=False):
        if verbose:
            print("aaindex2 infile:{}".format(self._data_location))

    def get_parser_name(self):
        return self._parsername


class AAindex(object):
    """Proxy class for AAindex parsers

    usage:

    1. create an object

    .. ipython:: python

        import bilab
        aa1prop = bilab.aaprop.AAindex(parser=bilab.aaprop.AAindex1Parser)

    data: stored in aa1prop._data  

    **Keys** are   

        'H': string,

        'D': string,    

        'R': dict(),    

        'A': string,    

        'T': string,    

        'J': string,    

        'C': dict()  correlated aaindex,  

        'I': dict()  index values associated with keys of 20 amino acids  
     

    2. get scale values

    .. ipython:: python

        aa1prop.get_scale('KYTJ820101')

    """

    def __init__(self, parser=AAindex1Parser, adapted_methods=None):
        """
           aaindex = 1,2,3 means aaindex1, aaindex2, aaindex3
        """
        super(AAindex, self).__init__()
        self._parser = parser()
        if adapted_methods is not None:
            self.__dict__.update(adapted_methods)

    def __getattr__(self, attr):
        """ get delegation to the object """
        return getattr(self._parser, attr)

