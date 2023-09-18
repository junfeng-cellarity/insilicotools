#!/usr/bin/env python
from openeye.oechem import *
import sys

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: %s input.sdf output.sdf"%sys.argv[0])
    else:
        ifs = oemolistream()
        ifs.open(sys.argv[1])

        ofs = oemolostream()
        ofs.open(sys.argv[2])

        dict = {}
        mol = OEGraphMol()
        while OEReadMolecule(ifs,mol):
            if mol.GetTitle() not in dict:
                dict[mol.GetTitle()] = 1
                OEWriteMolecule(ofs,mol)
        ifs.close()
        ofs.close()
