#!/usr/bin/env python
from openeye.oechem import *
import sys

if __name__ == "__main__":
    if len(sys.argv)!=4:
        print "Usage:%s mollist input.sdf output.sdf"%sys.argv[0]
    else:
        list = sys.argv[1]
        molList = open(list, "r").read().splitlines()

        sdf = sys.argv[2]

        ifs = oemolistream()
        ifs.open(sdf)

        output = sys.argv[3]
        ofs = oemolostream()
        ofs.open(output)

        mol = OEGraphMol()
        while OEReadMolecule(ifs,mol):
            if mol.GetTitle() in molList:
                OEWriteMolecule(ofs,mol)
        ofs.close()
        ifs.close()
