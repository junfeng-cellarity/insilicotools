#!/usr/bin/env python
import sys
from openeye.oechem import *

if len(sys.argv)!=3:
    print ("Usage: %s input.smi output.sdf")
else:
    ifs = oemolistream()
    ifs.open(sys.argv[1])
    
    ofs = oemolostream()
    ofs.open(sys.argv[2])

    mol = OEGraphMol()
    while OEReadMolecule(ifs,mol):
        OEWriteMolecule(ofs,mol)
    ifs.close()
    ofs.close()