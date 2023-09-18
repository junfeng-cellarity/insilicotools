#!/usr/bin/env python
from openeye.oechem import *
import sys

if __name__ == "__main__":
    if len(sys.argv)!=5:
        print "Usage:%s tagName mollist input.sdf output.sdf"%sys.argv[0]
    else:
        tag = sys.argv[1]
        list = sys.argv[2]
        molList = open(list, "r").read().splitlines()

        sdf = sys.argv[3]

        ifs = oemolistream()
        ifs.open(sdf)

        output = sys.argv[4]
        ofs = oemolostream()
        ofs.open(output)

        dict = {}
        mol = OEGraphMol()
        while OEReadMolecule(ifs,mol):
            tagValue = OEGetSDData(mol,tag)
            if tagValue in molList:
                dict[tagValue] = OEGraphMol(mol)
        for molname in molList:
            if dict.has_key(molname):
                OEWriteMolecule(ofs,dict[molname])
        ofs.close()
        ifs.close()
