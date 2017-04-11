bilab package
=============

Prerequisites
-------------

1. Numpy
2. Distutils2
3. C/C++ compiler

Download and Installation
-------------------------
1. Numpy

    git clone git://github.com/numpy/numpy.git numpy  
    cd numpy  
    python setup.py build install  

2. Distutil2 and setuptools  
    cd Distutils2-1.0a4  
    python setup.py build install  

3. Boost Library
   Please refer to http://www.boost.org

4. bilab
    cd bilab
    pysetup run build install_dist
or  specify the installation prefix
    pysetup run build install_dist --prefix=/path/to/install

Test
-----
    ./prInteract.py --pdbdir . --pdbid 9mht --outfmt txt --out 9mht --distance 5.0 --source protein --target nucleic -v  

Output
--------------
The result is saved to the file named 9mht.txt    

\#ATOM: source : target : the distance : sum of covalent bond length :is a covalent bond (True/False)  
\#format of source: resname.chid.Icode.resnum.atomname  
\#format of target: resname.chid.Icode.resnum.atomname  
ATOM:TYR.A.\_.44.OH:DC.C.\_.402.P:4.612:1.730:False    
ATOM:TYR.A.\_.44.CE2:DC.C.\_.402.P:4.829:1.830:False    
ATOM:TYR.A.\_.44.OH:DC.C.\_.402.OP2:4.023:1.320:False    
ATOM:TYR.A.\_.44.CZ:DC.C.\_.402.OP2:4.305:1.420:False    
ATOM:TYR.A.\_.44.CE2:DC.C.\_.402.OP2:3.622:1.420:False    
