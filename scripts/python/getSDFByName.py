#!/usr/bin/env python
import sys
from openeye.oechem import *

if len(sys.argv)!=3:
    print("Usage:%s molName.list input.sdf > output.sdf"%sys.argv[0])
else:
    nameList = open(sys.argv[1],"r").read().split("\n")
    mol = OEGraphMol()
    ifs = oemolistream()
    ifs.open(sys.argv[2])
    ofs = oemolostream()
    ofs.SetFormat(OEFormat_SDF)
    ofs.openstring()
    while OEReadMolecule(ifs,mol):
        if mol.GetTitle() in nameList:
            OEWriteMolecule(ofs,mol)
    ifs.close()
    print(ofs.GetString())