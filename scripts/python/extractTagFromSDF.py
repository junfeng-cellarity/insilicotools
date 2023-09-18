#!/usr/bin/env python
from openeye.oechem import *
import sys

if __name__ == "__main__":
    if len(sys.argv)!=3:
        print "Usage:%s tag input.sdf"%sys.argv[0]
    else:
        tag = sys.argv[1]
        sdf = sys.argv[2]

        ifs = oemolistream()
        ifs.open(sdf)
        mol = OEGraphMol()
        while OEReadMolecule(ifs,mol):
            print OEGetSDData(mol,tag)
        ifs.close()
