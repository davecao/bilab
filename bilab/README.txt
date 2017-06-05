bilab package
=============

Prerequisites
-------------

1. Numpy
2. distutil and setuptools
3. C/C++ compiler
4. Boost C++ Libaray

Download and Installation
-------------------------
1. Numpy

    git clone git://github.com/numpy/numpy.git numpy
    cd numpy
    python setup.py build install

2. Boost Library

3. bilab
    cd bilab
    python setup.py install install_data

Test
-----
    ./prInteract.py --pdbdir . --pdbid 9mht --outfmt txt --out 9mht --distance 5.0 --source protein --target nucleic -v

  the output file named 9mht.txt
----------------------------------------------------------------------------------------------------
#ATOM: source : target : the distance : sum of covalent bond length :is a covalent bond (True/False)
# format of source: resname.chid.Icode.resnum.atomname
# format of target: resname.chid.Icode.resnum.atomname
ATOM:TYR.A._.44.OH:DC.C._.402.P:4.612:1.730:False
ATOM:TYR.A._.44.CE2:DC.C._.402.P:4.829:1.830:False
ATOM:TYR.A._.44.OH:DC.C._.402.OP2:4.023:1.320:False
ATOM:TYR.A._.44.CZ:DC.C._.402.OP2:4.305:1.420:False
ATOM:TYR.A._.44.CE2:DC.C._.402.OP2:3.622:1.420:False



