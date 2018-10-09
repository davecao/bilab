# -*- coding: utf-8 -*-
# sparse6.py

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

__all__ = []


def adj2sparse(adj_list):
    """
    Encode an adjacency list to sparse6 format
    
    Graph 1, order 7.
      0 : 1 2;
      1 : 2;
      2 : ;
      3 : ;
      4 : ;
      5 : 6;
      6 : ;
    
    Sparse6 format:
    ----------------
    >>sparse6<<
    :Fa@x^
    ----------------
    ':' indicates sparse6 format.
      Subtract 63 from the other bytes and write them in binary, 
      six bits each.

       000111 100010 000001 111001 011111

      The first byte is not 63, so it is n.  n=7
      n-1 needs 3 bits (k=3).  Write the other bits in groups
      of 1 and k:

        1 000  1 000  0 001  1 110   0 101  1 111
  
        This is the b/x sequence  1,0 1,0 0,1 1,6 0,5 1,7.
        The 1,7 at the end is just padding.
        The remaining parts give the edges 0-1 0-2 1-2 5-6.

    
    see the original description at the following URL
        https://users.cecs.anu.edu.au/~bdm/data/formats.txt
    
    
    Parameters
    ----------
    adj_list : dict
        The adjacency list of a graph
    
    Returns
    -------
    str
        
    """
    pass


def sparse2adj(sp, full=False):
    """
    Convert sparse6 format to a dict, representing an adjacency list of
    a graph.
    
    Parameters
    ----------
    sp : a list of str
        Convert to adjacency lists of a graph
    
    Returns
    -------
    A dict
         
    """
    pass