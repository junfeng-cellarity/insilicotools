#!/usr/bin/env python
from openeye.oechem import *
import sys
if __name__ == "__main__":
    if len(sys.argv)!=5:
        print ("Usage:%s input.sdf tag_name tags.txt output.sdf"%sys.argv[0])
    else:
        ifs = oemolistream()
        ifs.open(sys.argv[1])
        ofs = oemolostream()
        ofs.open(sys.argv[4])
        mol = OEGraphMol()
        tag_name = sys.argv[2]
        tags = open(sys.argv[3],"r").read().splitlines()
        while OEReadMolecule(ifs,mol):
            name = OEGetSDData(mol,tag_name)
            if name in tags:
                OEWriteMolecule(ofs,mol)
        ifs.close()
        ofs.close()
