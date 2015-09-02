# -*- coding: utf-8 -*-

import os.path
import errno
import sys
import imp
import string
import numpy as np

from threading import Thread
try:
    from Queue import Queue, Empty
except ImportError:
    from queue import Queue, Empty

__all__ = [ 'AsynchronousFileReader',
            'prettify_xml',
            'translator',
            'chars_occ_in_str',
            'import_module',
            'create_dir',
            'exist_dir',
            'three2one',
            'get_loaded_modules',
            'search_loaded_module',
            'find_loaded_module']
"""
      `ASX`_ (B)  asparagine or aspartic acid
      `GLX`_ (Z)  glutamine or glutamic acid
      `CSO`_ (C)  S-hydroxycysteine
      `HIP`_ (H)  ND1-phosphohistidine
       HSD   (H)  prototropic tautomer of histidine, H on ND1 (CHARMM)
       HSE   (H)  prototropic tautomer of histidine, H on NE2 (CHARMM)
       HSP   (H)  protonated histidine
      `MSE`_      selenomethionine
      `SEC`_ (U)  selenocysteine
      `SEP`_ (S)  phosphoserine
      `TPO`_ (T)  phosphothreonine
      `PTR`_ (Y)  O-phosphotyrosine
       XLE   (J)  leucine or isoleucine
       XAA   (X)  unspecified or unknown
"""
code = {
    "GLY" : "G",
    "ALA" : "A",
    "LEU" : "L",
    "ILE" : "I",
    "ARG" : "R",
    "LYS" : "K",
    "MET" : "M",
    "CYS" : "C",
    "TYR" : "Y",
    "THR" : "T",
    "PRO" : "P",
    "SER" : "S",
    "TRP" : "W",
    "ASP" : "D",
    "GLU" : "E",
    "ASN" : "N",
    "GLN" : "Q",
    "PHE" : "F",
    "HIS" : "H",
    "VAL" : "V",
    
    "M3L" : "K",
    "MSE" : "X",
    "CAS" : "X",
    "ASX" : "B",
    "GLX" : "Z",
    "XJE" : "J",
    "HIP" : "X",
    "HSD" : "X",
    "HSE" : "X",
    "HSP" : "X",
    "CSO" : "X",
    "SEC" : "X",
    "SEP" : "X",
    "TPO" : "X",
    "PTR" : "X",
    "XAA" : "X"
}


def three2one(aa):
    """
        Convert three letters of amino acids into one letter
    """
    aa_upper = aa.upper()
    if aa_upper in code:
        return code[aa_upper]
    return None


def get_loaded_modules():
    """
        print out all loaded modules in sys
    """
    for key in sys.modules.keys():
        print("{}".format(key))

def search_loaded_module(str):
    """
        Search a module started with str within the sys

    """
    found = False
    for key in sys.modules.keys():
        if key.startswith(str):
            found = True
            print("{}".format(key))
    if not found:
        print("Could not found module startswith {}".format(str))
    return found

def find_loaded_module(mod):
    """
       find module with exact name in sys
    """
    #for key in sys.modules.keys():
    #    print key
    found = False
    if mod in sys.modules.keys():
        print("found {} in sys".format(mod))
        found = True
    else:
        print("Could not found {} in sys".format(mod))
    return found

def exist_dir(path, verbose=True):
    """ Check the given path is a directory or not """
    test = True
    if not os.path.isdir(path):
        if verbose:
            print("The specified path is not a directory:{0}".format(path))
        test = False
    else:
        if not os.path.exists(path):
            if verbose:
                print("The directory {0} is not existed.".format(path))
            test = False
    return test

def create_dir(path):
    if not exist_dir(path):
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno!= errno.EEXIST:
                raise

def import_module(modname, forall=True, path=".", verbose=False):
    """
        import package dynamically
    """
    try:
        return sys.modules[modname]
    except KeyError:
        pass

    sys.path.append(path)
    fp, pathname, description = imp.find_module(modname)
    sys.path.pop(-1)
    try:
        #print("Path:{}".format(path))
        # load the module
        mod = imp.load_module(modname, fp, pathname, description)
        if verbose:
            print("Loading module {} from {}".format(modname, pathname))
        if forall:
            # add it to global space
            globals().update(mod.__dict__)
            if verbose:
                print("Add loaded module {} to global name space"
                    .format(modname))
        return mod
    finally:
        if fp:
            fp.close()

def prettify_xml(elem):
    """
        Pretty print xml
        referred to pymotw.com
    """
    try:
        from xml.etree import ElementTree
        from xml.dom import minidom
    except ImportError:
        raise ImportError('Failed to import xml.etree and xml.dom.')

    rough_string = ElementTree.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")

def translator(frm='', to='', delete='', keep=None):
    # Python Cookbook Recipe 1.9
    # Chris Perkins, Raymond Hettinger
    # usage:
    #    >>std_aa_only = translator(keep="ACDEFGHILKMNPQRSTVWY")
    #    >>std_aa_only("ADDGSFGHXXX")
    #    ADDGSFGH  <--- XXX will be delete
    #    >>
    if len(to) == 1: to = to * len(frm)
    trans = string.maketrans(frm, to)
    if keep is not None:
        allchars = string.maketrans('', '')
        # delete is expanded to delete everything except
        # what is mentioned in set(keep)-set(delete)
        delete = allchars.translate(allchars, keep.translate(allchars, delete))
    def translate(s):
        return s.translate(trans, delete)
    return translate

def chars_occ_in_str(str):
    """
        Count the occurrences of characters in a string by numpy
        see: https://gist.github.com/zed/347000
    """
    if isinstance(str, bytes):
        str = str.decode("utf-8")
    enc = 'utf-16' + ('le' if sys.byteorder == 'little' else 'be')
    a = np.frombuffer(str.encode(enc), dtype = np.uint16)
    counts = np.bincount(a)

    counts = [(unichr(i), v) for i, v in enumerate(counts) if v]
    return dict(counts)

class AsynchronousFileReader(Thread):
    """
    Helper class to implement asynchronous reading of a file
    in a separate thread. Pushes read lines on a queue to
    be consumed in another thread.
    """

    def __init__(self, fd, queue):
        assert isinstance(queue, Queue)
        assert callable(fd.readline)
        Thread.__init__(self)
        self._fd = fd
        self._queue = queue

    def run(self):
        """
          The body of the tread: read lines and put them on the queue.
        """
        for line in iter(self._fd.readline, ''):
            self._queue.put(line)

    def eof(self):
        """Check whether there is no more content to expect."""
        return not self.is_alive() and self._queue.empty()

class ProgressBar(Thread):
    """ This class is borrowed from
        Corey Goldberg - 2010
       ascii command-line progress bar with percentage and elapsed time display
    """
    def __init__(self, threadId, name, interval=5, fill_char='#'):
        super(ProgressBar, self).__init__()
        self.threadId = threadId
        self.name = name
        self.prog_bar = '[]'
        self.fill_char = fill_char
        self.basic_char = '.'
        self.counter = 0
        self.width = 40
        self.start_time = time.time()
        self.interval = int(interval)
        #self.__update_amount(0)

    def reset(self):
        self.start_time = time.time()

    def elapsed_time(self):
        elapsed = int(time.time() - self.start_time)
        return elapsed

    def run(self):
        if sys.platform.lower().startswith('win'):
            print("{}\r".format(self))
            #print self, '\r',
        else:
            print("{}[A".format(chr(27)))
            #print self, chr(27) + '[A'
        self.update_time(self.elapsed_time())
        time.sleep(self.interval)
        self.counter += self.interval
        print("{}".format(self))

    def update_time(self, elapsed_secs):
        self.__update_amount(elapsed_secs)
        #self.prog_bar += '  %ds/%ss' % (elapsed_secs, self.duration)
        self.prog_bar += '  %ds' % (elapsed_secs)

    def __update_amount(self, new_amount):
        #print("start:{}".format(new_amount))
        percent_done = int(round((self.counter / 100.0) * 100.0))
        all_full = self.width - 2
        pct_place = self.counter
        pct_string = self.fill_char
        self.prog_bar = '[' + self.basic_char * all_full + ']'
        #num_hashes = int(round((percent_done / 100.0) * all_full))
        #pct_place = (len(self.prog_bar) / 2) - len(str(percent_done))
        #pct_string = '%d%%' % percent_done
        #print("pct:{}".format(pct_place))
        self.prog_bar = self.prog_bar[0:pct_place] + \
            (pct_string + self.prog_bar[pct_place + len(pct_string):])
        self.counter += 1

    def __reset(self):
        self.prog_bar = '[]'

    def __str__(self):
        return str(self.prog_bar)
