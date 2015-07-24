# -*- coding: utf-8 -*-

from __future__ import print_function

import os.path
import re
import pprint
import numpy as np

__all__ = ["PSIC"]


class PSIC(object):
    """
        class PSIC read sequence multiple alignment file
    """
    def __init__(self, aligned):
        super(PSIC, self).__init__()
        
