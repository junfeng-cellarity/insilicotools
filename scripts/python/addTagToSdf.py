#!/usr/bin/env python
from openeye.oechem import *
import sys
import os

if __name__=="__main__":
    if len(sys.argv)!=5:
        print ("Add a common tag to the sdf file")
        print ("Usage:%s input.sdf output.sdf tagName tagValue"%(os.path.basename(sys.argv[0])))
    else:
        ifs = oemolistream()
        ifs.open(sys.argv[1])

        ofs = oemolostream()
        ofs.open(sys.argv[2])

        tagName = sys.argv[3]
        tagValue = sys.argv[4]

        mol = OEGraphMol()
        while OEReadMolecule(ifs, mol):
            OESetSDData(mol,tagName,tagValue)
            OEWriteMolecule(ofs,mol)
        ifs.close()
        ofs.close()
