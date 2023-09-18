#!/usr/bin/env python
from openeye.oechem import *
import sys
import math

if len(sys.argv)!=3:
    print("Usage:%s input.sdf numOfTrunks"%sys.argv[0])
else:
    ifs = oemolistream()
    ifs.open(sys.argv[1])

    numTrunks = int(sys.argv[2])

    moldb = OEMolDatabase(ifs)
    numMols = moldb.GetMaxMolIdx()
    idx = 0
    trunkSize = int(math.ceil(float(numMols)/numTrunks))
    for i in range(0,numTrunks):
        fname = "MolTrunk%d.sdf"%i
        ofs = oemolostream()
        ofs.open(fname)
        for idx in range(i*trunkSize,(i+1)*trunkSize):
            if idx >= numMols:
                break
            mol1 = OEGraphMol()
            if moldb.GetMolecule(mol1,idx):
                OEWriteMolecule(ofs,mol1)
        ofs.close()


