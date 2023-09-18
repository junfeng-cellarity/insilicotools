#!/usr/bin/env python
import sys
from openeye.oechem import *

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print "Usage: %s origin.sdf target.sdf"%sys.argv[0]
    else:
        dict = {}
        origin = []
        target = []
        ifs = oemolistream()
        mol = OEGraphMol()
        ifs.open(sys.argv[1])
        while OEReadMolecule(ifs,mol):
            origin.append(OEGraphMol(mol))
        ifs.close()

        ifs.open(sys.argv[2])
        while OEReadMolecule(ifs,mol):
            target.append(OEGraphMol(mol))
        ifs.close()

        duplicate = []
        w_o_duplicate = []

        for mol1 in origin:
            canSmi = OECreateCanSmiString(mol1)
            for mol2 in target:
                canSmi2 = OECreateCanSmiString(mol2)
                if canSmi == canSmi2:
                    print canSmi2,mol2.GetTitle(), mol1.GetTitle()
                    break


